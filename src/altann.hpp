// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#pragma once
// AltAnn's compatibility model uses seven coarse states. Metadata belongs to
// each retained path and does not create additional states.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
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
struct ModelInputs {
    const vector<array<float, NUM_STATES>> &emit;
    const vector<double> &gt_score;
    const vector<double> &ag_score;
    const vector<double> &atg_score;
    const vector<double> &stop_score;
    const vector<char> &seq;
    const vector<vector<double>> &trans;
};

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
};
struct Fasta { string id; vector<char> sequence; };
struct Baseline {
    vector<int> states;
    // Metadata of the winning N prefix at each base, needed for capped windows.
    vector<PathMetadata> noncoding_prefix;
};
Fasta read_fasta(const string &path);
vector<array<float, NUM_STATES>> load_emissions(const string &, int);
vector<double> load_sparse_scores(const string &, int);
vector<vector<double>> init_transitions();
int splice_motif_start_0based(int);
int codon_start_0based(int);
TransitionResult evaluate_transition(int, int, int, const PathMetadata &, const ModelInputs &);
Baseline decode_baseline(const ModelInputs &);
Baseline load_baseline_log(const ModelInputs &, const string &);
vector<TranscriptAnnotation> build_transcript_annotations(const vector<int> &, const vector<int> &, const string &);
vector<TranscriptAnnotation> build_transcript_annotations(const vector<int> &, const string &);
string intron_chain_key(const TranscriptAnnotation &);
string full_transcript_key(const TranscriptAnnotation &);
string classify_splice_difference(const TranscriptAnnotation &, const TranscriptAnnotation &);
string tsv_escape_chain(const TranscriptAnnotation &);
void write_hierarchical_transcript(ostream &, const TranscriptAnnotation &, const string &, const string &);
bool run_local_k_best(const string &, const vector<TranscriptAnnotation> &, const vector<int> &, const vector<PathMetadata> &, const ModelInputs &, const Options &);

// Stage the complete output bundle; preserve prior files if publication fails.
void write_output_bundle(const vector<pair<string, string>> &files);
