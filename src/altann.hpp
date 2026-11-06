// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#pragma once
// AltAnn's compatibility model uses seven coarse states. Metadata belongs to
// each retained path and does not create additional states.
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
using namespace std;
static const int NUM_STATES = 7;

static const array<string, NUM_STATES> state_name = {
    "N", "E0", "E1", "E2", "I0", "I1", "I2"
};

// This finite sentinel is part of the inherited scoring convention. Baseline
// decoding deliberately preserves its arithmetic for compatibility.
static const double NEG_INF = -1e9;
static const int MIN_INTRON = 40;
static const int MIN_EXON   = 3;
static const int MIN_INTER  = 30;
static const int MIN_SINGLE = 100;
struct PathMetadata {
    int intron_len = 0;
    int exon_len = 0;
    int inter_len = 0;
    int exon_from = 0;
    // The predecessor of the cell carrying this metadata.  The legacy start
    // rule checks this value to ensure the intergenic run itself came from N.
    int predecessor_state = -1;
};

struct GenomicInterval {
    int start_0based = 0;
    int end_0based_exclusive = 0;
    bool operator==(const GenomicInterval &other) const {
        return start_0based == other.start_0based &&
               end_0based_exclusive == other.end_0based_exclusive;
    }
    bool operator<(const GenomicInterval &other) const {
        return tie(start_0based, end_0based_exclusive) <
               tie(other.start_0based, other.end_0based_exclusive);
    }
};

struct SpliceJunction {
    // Both coordinates are the first intronic bases of their motifs.  The
    // interval is [donor_0based, acceptor_0based + 2).
    int donor_0based = -1;
    int acceptor_0based = -1;
    char strand = '+';
    bool operator==(const SpliceJunction &other) const {
        return donor_0based == other.donor_0based &&
               acceptor_0based == other.acceptor_0based &&
               strand == other.strand;
    }
    bool operator<(const SpliceJunction &other) const {
        return tie(donor_0based, acceptor_0based, strand) <
               tie(other.donor_0based, other.acceptor_0based, other.strand);
    }
};

struct TranscriptAnnotation {
    string sequence_id;
    char strand = '+';
    int genomic_start_0based = -1;
    int genomic_end_0based_exclusive = -1;
    int path_start = -1;
    int path_end = -1;
    vector<GenomicInterval> exons;
    vector<GenomicInterval> cds_intervals;
    vector<SpliceJunction> junctions;
    double score = NEG_INF;
};

inline bool is_exon(int s) {
    return (s == 1 || s == 2 || s == 3);
}

inline bool is_intron(int s) {
    return (s == 4 || s == 5 || s == 6);
}

// Producers round emissions to decimal text; retain UniAnn's float storage.
// Accumulated path scores and position-specific transition scores use double.
template <class T>
struct PositionWindow {
    const T *values = nullptr;
    int first = 0;
    int count = 0;
    const T &operator[](int position) const {
        assert(position >= first && position - first < count);
        return values[position - first];
    }
};

// Scores are only readable inside one loaded window, addressed by sequence
// position. Sequence and transition tables stay available for the whole record.
struct ModelInputs {
    PositionWindow<array<float, NUM_STATES>> emit;
    PositionWindow<double> gt_score;
    PositionWindow<double> ag_score;
    PositionWindow<double> atg_score;
    PositionWindow<double> stop_score;
    const vector<char> &seq;
    const vector<vector<double>> &trans;
};

// One loaded interval of the score tables, with the same layout as ModelInputs.
class ScoreWindow {
public:
    ModelInputs inputs() const;
private:
    friend class ScoreSource;
    const vector<char> *seq_ = nullptr;
    const vector<vector<double>> *trans_ = nullptr;
    int first_ = 0;
    vector<array<float, NUM_STATES>> emit_;
    array<vector<double>, 4> sparse_;  // GT, AG, ATG, stop
};

// Validates the score tables once, then serves any interval on demand. Emissions
// stay in the input file behind a position index; only the motif tables, which
// are sparse, are held in memory.
class ScoreSource {
public:
    struct ByteRange { uint64_t begin = numeric_limits<uint64_t>::max(); uint64_t end = 0; };
    ScoreSource(const vector<char> &seq, const vector<vector<double>> &trans,
                const string &emissions, const array<string, 4> &sparse);
    int length() const { return static_cast<int>(seq_.size()); }
    // Interval [first, last], clipped to the sequence. Thread-safe.
    ScoreWindow window(int first, int last) const;
private:
    struct SparseScores { vector<int> positions; vector<double> values; };
    const vector<char> &seq_;
    const vector<vector<double>> &trans_;
    string emissions_path_;
    vector<ByteRange> emission_ranges_;
    array<SparseScores, 4> sparse_;
};

// Forward walk over the sequence, holding one block of scores at a time. Every
// emission row is parsed and validated once by a complete walk.
class ScoreCursor {
public:
    explicit ScoreCursor(const ScoreSource &source) : source_(source) {}
    // Scores for a transition into position, including its motif offsets.
    ModelInputs at(int position);
private:
    const ScoreSource &source_;
    ScoreWindow window_;
    int loaded_last_ = -1;
};

// Decoding a transition reads the destination, its successor (TAG/AG emission
// rule) and up to two earlier bases (motif starts).
inline ScoreWindow transition_window(const ScoreSource &source, int first_position, int last_position) {
    return source.window(first_position - 2, last_position + 1);
}

struct TransitionResult {
    bool allowed = false;
    double transition_score = NEG_INF;
    double emission_score = NEG_INF;
    PathMetadata next_metadata;
};

struct Options {
    int k = 10;
    int flank = 1000;
    int threads = 0;
    string output_filename;
    string report_filename;
    string viterbi_log_filename;
    string models_filename;
    bool raw = false;
    // Shared by all concurrent loci; tracebacks that do not fit spill to a
    // temporary file beside the output.
    size_t traceback_memory_bytes = size_t(1024) * 1024 * 1024;
};
struct Fasta { string id; vector<char> sequence; };
// Run-length encoded state path over a whole sequence.
class StatePath {
public:
    StatePath() = default;
    // Runs in ascending start order; the first run starts at position 0.
    StatePath(vector<int> run_starts, vector<int8_t> run_states, int length)
        : run_starts_(std::move(run_starts)), run_states_(std::move(run_states)), length_(length) {}
    int size() const { return length_; }
    int operator[](int position) const {
        const auto next = upper_bound(run_starts_.begin(), run_starts_.end(), position);
        return run_states_[next - run_starts_.begin() - 1];
    }
    const vector<int> &run_starts() const { return run_starts_; }
    const vector<int8_t> &run_states() const { return run_states_; }
private:
    vector<int> run_starts_;
    vector<int8_t> run_states_;
    int length_ = 0;
};

// Metadata of the winning N cell at each base, needed for capped windows. An N
// cell only carries inter_len and its predecessor, and N->N steps just increment
// inter_len, so only the positions breaking that pattern are stored.
class NoncodingPrefix {
public:
    void append(int position, const PathMetadata &metadata) {
        if (!anchors_.empty() && metadata.predecessor_state == 0 &&
            metadata.inter_len == anchors_.back().inter_len + position - anchors_.back().position)
            return;
        anchors_.push_back({position, metadata.inter_len, metadata.predecessor_state});
    }
    PathMetadata operator[](int position) const {
        const auto next = upper_bound(anchors_.begin(), anchors_.end(), position,
            [](int value, const Anchor &anchor) { return value < anchor.position; });
        const Anchor &anchor = *(next - 1);
        PathMetadata metadata;
        metadata.inter_len = anchor.inter_len + position - anchor.position;
        metadata.predecessor_state = position == anchor.position ? anchor.predecessor_state : 0;
        return metadata;
    }
private:
    struct Anchor { int position; int inter_len; int predecessor_state; };
    vector<Anchor> anchors_;
};

struct Baseline {
    StatePath states;
    NoncodingPrefix noncoding_prefix;
};
Fasta read_fasta(const string &path);
vector<vector<double>> init_transitions();
int splice_motif_start_0based(int);
int codon_start_0based(int);
TransitionResult evaluate_transition(int, int, int, const PathMetadata &, const ModelInputs &);
Baseline decode_baseline(const ScoreSource &);
Baseline load_baseline_log(const ScoreSource &, const string &);
vector<TranscriptAnnotation> build_transcript_annotations(const vector<int> &, const vector<int> &, const string &);
vector<TranscriptAnnotation> build_transcript_annotations(const StatePath &, const string &);
string intron_chain_key(const TranscriptAnnotation &);
string full_transcript_key(const TranscriptAnnotation &);
string classify_splice_difference(const TranscriptAnnotation &, const TranscriptAnnotation &);
double splice_site_jaccard(const TranscriptAnnotation &, const TranscriptAnnotation &);
string tsv_escape_chain(const TranscriptAnnotation &);
void write_hierarchical_transcript(ostream &, const TranscriptAnnotation &, const string &, const string &);
void run_model_scoring(const string &, const vector<TranscriptAnnotation> &, const Baseline &, const ScoreSource &, const Options &);
bool run_local_k_best(const string &, const vector<TranscriptAnnotation> &, const Baseline &, const ScoreSource &, const Options &);

// Stage the complete output bundle; preserve prior files if publication fails.
void write_output_bundle(const vector<pair<string, string>> &files);
