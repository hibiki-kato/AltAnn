// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"
#include <omp.h>
#include <unistd.h>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <atomic>
#include <climits>
#include <future>
#include <chrono>

// Local decoding retains K paths per coarse state. Search, annotation ranking,
// and per-locus output are separate steps; the public runner only orchestrates
// parallel work and publishes the ordered output bundle.
namespace {
// K counts retained paths per coarse state, before annotation deduplication.
struct PathCandidate {
    double score;
    int from_state;
    int from_k;
    TransitionResult evaluated;

    bool operator<(const PathCandidate &other) const {
        if (score != other.score) return score > other.score;
        if (from_state != other.from_state) return from_state < other.from_state;
        return from_k < other.from_k;
    }
};

struct LocalKBestCandidate {
    int rank;
    double score;        // score of the whole local path this gene belongs to
    double score_delta;
    double actual_start_score;
    int genes_in_path = 1;
    int gene_in_path = 1;
    string classification;
    double splice_jaccard = 1.0;
    bool identical_to_reference;
    TranscriptAnnotation transcript;
};

struct ActiveCell {
    double dp = NEG_INF;
    PathMetadata metadata;
};

// K is limited to 32767 at the CLI, so signed 16-bit predecessor ranks
// can represent every retained path and the -1 unreachable sentinel.
struct CompactBP {
    int16_t bt_state;
    int16_t bt_k;
};

struct LocusOutput {
    string gff, reference, report, error;
    ExactSearchStats exact;
    double search_seconds = 0;
    int window_length = 0;
};

// Backpointer rows are written forward and read backward. A traceback that does
// not fit the shared budget moves through an unlinked temporary file in blocks,
// with each write and read-ahead running in the background while decoding
// continues. A spilling locus holds at most two blocks in memory.
class BackpointerStore {
public:
    BackpointerStore(size_t rows, size_t row_cells, TracebackBudget &budget,
                     const string &spill_directory,
                     size_t block_bytes = size_t(16) * 1024 * 1024)
        : budget_(budget), rows_(rows), row_cells_(row_cells) {
        const size_t row_bytes = row_cells * sizeof(CompactBP);
        spilling_ = !budget_.try_reserve(rows * row_bytes);
        reserved_ = spilling_ ? 0 : rows * row_bytes;
        rows_per_block_ = spilling_ ? min(rows, max<size_t>(1, block_bytes / row_bytes)) : rows;
        front_.assign(rows_per_block_ * row_cells, CompactBP{-1, -1});
        if (!spilling_) return;
        back_.resize(front_.size());
        string path = (spill_directory.empty() ? string(".") : spill_directory) +
                      "/.altann-traceback-XXXXXX";
        descriptor_ = mkstemp(path.data());
        if (descriptor_ < 0) throw runtime_error("Cannot create traceback file in " + spill_directory +
                                                 ": " + strerror(errno));
        unlink(path.c_str());
    }
    BackpointerStore(const BackpointerStore &) = delete;
    BackpointerStore &operator=(const BackpointerStore &) = delete;
    ~BackpointerStore() {
        if (pending_.valid()) pending_.wait();
        if (descriptor_ >= 0) close(descriptor_);
        budget_.release(reserved_);
    }

    bool spilling() const { return spilling_; }

    // Rows must be requested in ascending order. Unwritten cells hold {-1, -1}.
    CompactBP *write_row(size_t row) {
        const size_t block = row / rows_per_block_;
        if (block != current_block_) {
            write_block_in_background(current_block_);
            current_block_ = block;
            fill(front_.begin(), front_.end(), CompactBP{-1, -1});
        }
        return front_.data() + (row - block * rows_per_block_) * row_cells_;
    }

    // Call once after the last write_row; rows are then read in descending order.
    void finish_writes() {
        if (!spilling_) return;
        write_block_in_background(current_block_);
        pending_.get();
        current_block_ = NO_BLOCK;
    }

    const CompactBP *read_row(size_t row) {
        const size_t block = row / rows_per_block_;
        if (block != current_block_) {
            if (pending_.valid() && ahead_block_ == block) {
                pending_.get();
                front_.swap(back_);
            } else {
                if (pending_.valid()) pending_.get();
                read_block(block, front_);
            }
            current_block_ = block;
            if (block > 0) {
                ahead_block_ = block - 1;
                pending_ = async(launch::async, [this, ahead = block - 1] { read_block(ahead, back_); });
            }
        }
        return front_.data() + (row - block * rows_per_block_) * row_cells_;
    }

private:
    static constexpr size_t NO_BLOCK = numeric_limits<size_t>::max();

    size_t block_cells(size_t block) const {
        return min(rows_per_block_, rows_ - block * rows_per_block_) * row_cells_;
    }
    off_t block_offset(size_t block) const {
        return static_cast<off_t>(block * rows_per_block_ * row_cells_ * sizeof(CompactBP));
    }
    void write_block_in_background(size_t block) {
        if (pending_.valid()) pending_.get();
        front_.swap(back_);
        pending_ = async(launch::async, [this, block] {
            const char *cursor = reinterpret_cast<const char *>(back_.data());
            size_t remaining = block_cells(block) * sizeof(CompactBP);
            off_t offset = block_offset(block);
            while (remaining) {
                const ssize_t written = pwrite(descriptor_, cursor, remaining, offset);
                if (written < 0 && errno == EINTR) continue;
                if (written <= 0) throw runtime_error(string("Cannot write traceback file: ") + strerror(errno));
                cursor += written;
                remaining -= written;
                offset += written;
            }
        });
    }
    void read_block(size_t block, vector<CompactBP> &buffer) const {
        char *cursor = reinterpret_cast<char *>(buffer.data());
        size_t remaining = block_cells(block) * sizeof(CompactBP);
        off_t offset = block_offset(block);
        while (remaining) {
            const ssize_t count = pread(descriptor_, cursor, remaining, offset);
            if (count < 0 && errno == EINTR) continue;
            if (count <= 0) throw runtime_error(string("Cannot read traceback file: ") + strerror(errno));
            cursor += count;
            remaining -= count;
            offset += count;
        }
    }

    TracebackBudget &budget_;
    size_t reserved_ = 0;
    size_t rows_, row_cells_, rows_per_block_ = 0;
    bool spilling_ = false;
    int descriptor_ = -1;
    size_t current_block_ = 0, ahead_block_ = NO_BLOCK;
    vector<CompactBP> front_, back_;
    // Declared after its buffers: the destructor joins before they are freed.
    future<void> pending_;
};

// Decode the bounded interval with both endpoint states pinned to N. This
// function owns all local DP memory and returns paths in descending score order.
// Cells whose state is absent from allowed_states stay unreachable.
vector<LocalPathResult> search_paths_impl(
    int L_bound, int R_bound, const PathMetadata &initial_metadata,
    const ModelInputs &inputs, const Options &options, TracebackBudget &budget,
    const vector<uint8_t> *allowed_states = nullptr
) {
    const int K = options.k;
    const int local_len = R_bound - L_bound + 1;
    // The local decoder solves exactly the UniAnn problem on [L_bound, R_bound]
    // with both endpoints pinned to N: same states, same transitions, same
    // structural constraints. The difference from baseline decoding is that
    // each cell keeps the top K partial paths instead of the single best.
    vector<vector<ActiveCell>> prev_dp(NUM_STATES, vector<ActiveCell>(K));
    vector<vector<ActiveCell>> curr_dp(NUM_STATES, vector<ActiveCell>(K));
    const size_t cells_per_position = static_cast<size_t>(NUM_STATES) * K;
    BackpointerStore backpointers(local_len, cells_per_position, budget, spill_directory(options));

    prev_dp[0][0].dp = 0.0 + inputs.emit[L_bound][0];
    prev_dp[0][0].metadata = initial_metadata;

    for (int i = 1; i < local_len; i++) {
        int p = L_bound + i;
        CompactBP *row = backpointers.write_row(i);

        for (int to = 0; to < NUM_STATES; to++) {
            if (allowed_states && !((*allowed_states)[i] >> to & 1)) continue;
            vector<PathCandidate> candidates;

            for (int from = 0; from < NUM_STATES; from++) {
                // Start codons span p-2..p: a first-step N->E could
                // otherwise annotate one base before the capped window.
                if (options.flank >= 0 && from == 0 && is_exon(to) &&
                    codon_start_0based(p) < L_bound) continue;
                for (int k = 0; k < K; k++) {
                    if (prev_dp[from][k].dp <= NEG_INF) continue;

                    TransitionResult evaluated = evaluate_transition(
                        p, from, to, prev_dp[from][k].metadata, inputs);

                    if (evaluated.allowed) {
                        double score = prev_dp[from][k].dp +
                            evaluated.transition_score + evaluated.emission_score;
                        candidates.push_back({score, from, k, evaluated});
                    }
                }
            }

            sort(candidates.begin(), candidates.end());
            for (int k = 0; k < min(K, static_cast<int>(candidates.size())); k++) {
                curr_dp[to][k].dp = candidates[k].score;
                curr_dp[to][k].metadata = candidates[k].evaluated.next_metadata;
                row[to * K + k] =
                    CompactBP{static_cast<int16_t>(candidates[k].from_state),
                              static_cast<int16_t>(candidates[k].from_k)};
            }
        }
        prev_dp.swap(curr_dp);
        for (int to = 0; to < NUM_STATES; to++) {
            for (int k = 0; k < K; k++) {
                curr_dp[to][k].dp = NEG_INF;
            }
        }
    }

    backpointers.finish_writes();

    // All paths step back together, so one backward pass reads each row once.
    struct Walker { LocalPathResult result; int state; int k; vector<int> starts; vector<int8_t> states; };
    vector<Walker> walkers;
    for (int k = 0; k < K; k++) {
        if (prev_dp[0][k].dp <= NEG_INF) continue;
        Walker walker{};
        walker.result.score = prev_dp[0][k].dp;
        walker.state = 0;
        walker.k = k;
        walkers.push_back(std::move(walker));
    }
    for (int i = local_len - 1; i > 0 && !walkers.empty(); i--) {
        const int p = L_bound + i;
        const CompactBP *row = backpointers.read_row(i);
        size_t kept = 0;
        for (size_t w = 0; w < walkers.size(); ++w) {
            Walker &walker = walkers[w];
            const CompactBP bp = row[walker.state * K + walker.k];
            if (bp.bt_state < 0) continue;
            if (bp.bt_state == 0 && is_exon(walker.state))
                walker.result.start_scores.push_back(inputs.atg_score[codon_start_0based(p)]);
            if (bp.bt_state != walker.state) {
                walker.starts.push_back(p);
                walker.states.push_back(static_cast<int8_t>(walker.state));
            }
            walker.state = bp.bt_state;
            walker.k = bp.bt_k;
            if (kept != w) walkers[kept] = std::move(walker);
            ++kept;
        }
        walkers.resize(kept);
    }

    vector<LocalPathResult> paths;
    for (auto &walker : walkers) {
        walker.starts.push_back(L_bound);
        walker.states.push_back(static_cast<int8_t>(walker.state));
        reverse(walker.starts.begin(), walker.starts.end());
        reverse(walker.states.begin(), walker.states.end());
        reverse(walker.result.start_scores.begin(), walker.result.start_scores.end());
        walker.result.path = StatePath(std::move(walker.starts), std::move(walker.states), local_len);
        paths.push_back(std::move(walker.result));
    }

    return paths;
}

// Coding intervals are ordered by genomic coordinate. Require a shared coding
// base, rather than overlapping transcript spans: a separate gene can sit
// entirely inside this reference's intron or within its search flank.
bool shares_cds_base(const TranscriptAnnotation &reference,
                     const TranscriptAnnotation &candidate) {
    if (reference.sequence_id != candidate.sequence_id ||
        reference.strand != candidate.strand) return false;
    size_t i = 0, j = 0;
    while (i < reference.cds_intervals.size() &&
           j < candidate.cds_intervals.size()) {
        const auto &a = reference.cds_intervals[i];
        const auto &b = candidate.cds_intervals[j];
        if (max(a.start_0based, b.start_0based) <
            min(a.end_0based_exclusive, b.end_0based_exclusive)) return true;
        if (a.end_0based_exclusive <= b.end_0based_exclusive) ++i;
        else ++j;
    }
    return false;
}

// Distinct state paths may describe the same retained annotation. Filter each
// path to genes sharing CDS bases with its own global-best reference, then
// collapse duplicates while retaining the original whole-path score order.
// Raw mode skips both steps: every gene of every path is returned, and rank is
// the decoder's path index, so a path without genes still consumes a rank.
vector<LocalKBestCandidate> rank_annotations(
    const vector<LocalPathResult> &paths, const TranscriptAnnotation &ref,
    const string &seqid, double ref_fixed_score, int K, bool raw = false
) {
    vector<LocalKBestCandidate> unique_candidates;
    set<string> seen;

    // A local path is one candidate. Unrelated genes do not reach the output or
    // consume ranks. Multiple overlapping genes still share one path rank.
    int rank_count = 0;
    for (const auto &path : paths) {
        if (raw) ++rank_count;
        auto anns = build_transcript_annotations(path.path, seqid);
        vector<size_t> retained;
        size_t overlapping = 0;
        for (size_t j = 0; j < anns.size(); ++j) {
            const bool overlaps = shares_cds_base(ref, anns[j]);
            overlapping += overlaps;
            if (raw || overlaps) retained.push_back(j);
        }
        if (retained.empty()) continue;

        if (!raw) {
            string key;
            for (size_t j : retained) key += full_transcript_key(anns[j]) + "|";
            if (!seen.insert(key).second) continue;
            ++rank_count;
        }

        const string ref_key = full_transcript_key(ref);
        for (size_t j : retained) {
            TranscriptAnnotation ann = anns[j];
            const bool overlaps = shares_cds_base(ref, ann);
            ann.score = path.score;

            LocalKBestCandidate cand;
            cand.rank = rank_count;
            cand.transcript = ann;
            cand.score = path.score;
            cand.score_delta = path.score - ref_fixed_score;
            cand.actual_start_score = j < path.start_scores.size()
                ? path.start_scores[j] : NEG_INF;
            // Preserve original path membership and start-score indexing even
            // when unrelated genes were omitted from the exported annotation.
            cand.genes_in_path = static_cast<int>(anns.size());
            cand.gene_in_path = static_cast<int>(j) + 1;
            // Two disjoint intronless genes would otherwise compare as identical.
            cand.classification = overlaps ? classify_splice_difference(ref, ann)
                                           : "no_cds_overlap";
            cand.splice_jaccard = overlaps ? splice_site_jaccard(ref, ann) : 0.0;
            cand.identical_to_reference =
                (overlapping == 1 && overlaps && full_transcript_key(ann) == ref_key);
            unique_candidates.push_back(cand);
        }
        if (rank_count >= K) break;
    }

    return unique_candidates;
}

// Establish the reference's window and prefix history, then score its best
// path and alternatives on precisely the same interval and initial condition.
LocusOutput decode_locus(
    size_t ref_idx, const string &seqid,
    const vector<TranscriptAnnotation> &references,
    const Baseline &baseline,
    const ScoreSource &source, const Options &options, TracebackBudget &budget
) {
    const StatePath &global_path_states = baseline.states;
    const int L = source.length();
    ostringstream gff_out, reference_out, tsv_out;
    tsv_out << setprecision(17);
    const auto &ref = references[ref_idx];

    int requested_left = (ref_idx == 0) ? 0 : references[ref_idx - 1].path_end;
    int requested_right = (ref_idx == references.size() - 1)
        ? L - 1 : references[ref_idx + 1].path_start - 1;
    if (options.flank >= 0) {
        int flank = options.flank;
        if (options.flank_cds_cap) {
            long long coding = 0;
            for (const auto &cds : ref.cds_intervals)
                coding += cds.end_0based_exclusive - cds.start_0based;
            flank = static_cast<int>(min<long long>(flank, coding));
        }
        int left_flank = flank, right_flank = flank;
        if (options.flank_junctions) {
            const int start = ref.genomic_start_0based, end = ref.genomic_end_0based_exclusive;
            const auto &junctions = *options.flank_junctions;
            // Relative support: a fraction of the median reads of this reference's
            // own CDS introns, or of the whole file when none of them is present.
            double min_reads = 0;
            if (options.flank_junction_min_fraction > 0) {
                vector<double> own;
                for (size_t k = 0; k + 1 < ref.cds_intervals.size(); ++k) {
                    const double n = junctions.reads(ref.cds_intervals[k].end_0based_exclusive,
                                                     ref.cds_intervals[k + 1].start_0based);
                    if (n >= 0) own.push_back(n);
                }
                double median = junctions.median_reads;
                if (!own.empty()) {
                    sort(own.begin(), own.end());
                    const size_t n = own.size();
                    median = n % 2 ? own[n / 2] : (own[n / 2 - 1] + own[n / 2]) / 2;
                }
                min_reads = options.flank_junction_min_fraction * median;
            }
            const auto supported = [&](const pair<int, int> &j) {
                return min_reads <= 0 || junctions.reads(j.first, j.second) >= min_reads;
            };
            // Left: the intron ends (acceptor-side exon starts) inside the span and
            // starts before it. by_end is sorted by end.
            // The donor-side exon ends at base it->first - 1, before the span.
            int farthest = INT_MAX;
            for (auto it = lower_bound(junctions.by_end.begin(), junctions.by_end.end(), make_pair(start, INT_MIN),
                     [](const pair<int, int> &j, const pair<int, int> &key) { return j.second < key.first; });
                 it != junctions.by_end.end() && it->second < end; ++it)
                if (it->first <= start && supported(*it)) farthest = min(farthest, it->first);
            if (farthest != INT_MAX)
                left_flank = max(left_flank, min(options.flank_junction_cap,
                                                 start + 1 - farthest + options.flank_junction_pad));
            // Right: the intron starts (donor-side exon ends) inside the span and
            // ends after it.
            // The acceptor-side exon starts at base it->second, after the span.
            farthest = INT_MIN;
            for (auto it = lower_bound(junctions.by_start.begin(), junctions.by_start.end(), make_pair(start + 1, INT_MIN));
                 it != junctions.by_start.end() && it->first <= end; ++it)
                if (it->second >= end && supported(*it)) farthest = max(farthest, it->second);
            if (farthest != INT_MIN)
                right_flank = max(right_flank, min(options.flank_junction_cap,
                                                   farthest + 1 - end + options.flank_junction_pad));
        }
        requested_left = max(requested_left, ref.genomic_start_0based - left_flank);
        requested_right = static_cast<int>(min<long long>(requested_right,
            static_cast<long long>(ref.genomic_end_0based_exclusive) - 1 + right_flank));
    }

    // Complete annotations omit terminal partial genes and invalid structures,
    // but those coding segments still exist in the global state path. Extend
    // each flank only through the contiguous N run adjoining this reference;
    // otherwise its requested endpoint could fall inside an omitted neighbor.
    int L_bound = ref.path_start - 1;
    int R_bound = ref.path_end;
    if (L_bound < 0 || R_bound >= L || L_bound >= R_bound ||
        global_path_states[L_bound] != 0 || global_path_states[R_bound] != 0) {
        throw runtime_error("Reference gene " + to_string(ref_idx + 1) +
                            " has invalid intergenic boundary anchors");
    }
    while (L_bound > requested_left && global_path_states[L_bound - 1] == 0)
        --L_bound;
    while (R_bound < requested_right && global_path_states[R_bound + 1] == 0)
        ++R_bound;

    const auto window = transition_window(source, L_bound, R_bound);
    const auto inputs = window.inputs();
    double ref_fixed_score = inputs.emit[L_bound][0];
    PathMetadata initial_metadata = options.flank >= 0
        ? baseline.noncoding_prefix[L_bound] : PathMetadata{};
    if (L_bound == 0 || options.flank < 0) {
        initial_metadata.predecessor_state = 0;
        initial_metadata.inter_len = 1;
    }
    PathMetadata ref_meta = initial_metadata;
    int ref_cur_state = 0;
    for (int p = L_bound + 1; p <= R_bound; p++) {
        int next_state = global_path_states[p];
        TransitionResult eval = evaluate_transition(p, ref_cur_state, next_state, ref_meta, inputs);
        if (!eval.allowed) {
            throw runtime_error("Forbidden reference transition for gene " + to_string(ref_idx + 1) + " at position " + to_string(p + 1));
        }
        ref_fixed_score += eval.transition_score + eval.emission_score;
        ref_meta = eval.next_metadata;
        ref_cur_state = next_state;
    }

    {
        auto reference = ref;
        reference.score = ref_fixed_score;
        const string locus = "locus" + to_string(ref_idx + 1);
        const string gene = locus + ".reference.g1";
        reference_out << ref.sequence_id << "\tUniAnn\tlocus\t" << L_bound + 1
                      << '\t' << R_bound + 1 << "\t.\t+\t.\tID=" << locus << '\n';
        reference_out << ref.sequence_id << "\tUniAnn\tgene\t" << ref.genomic_start_0based + 1
                      << '\t' << ref.genomic_end_0based_exclusive
                      << "\t.\t+\t.\tID=" << gene << ";Parent=" << locus << '\n';
        write_hierarchical_transcript(reference_out, reference, gene, gene + ".t1");
    }

    ExactSearchStats exact_stats;
    const auto search_begin = chrono::steady_clock::now();
    const auto paths = options.metadata_states
        ? search_paths_exact(L_bound, R_bound, initial_metadata, inputs, options, budget, &exact_stats)
        : search_paths_impl(L_bound, R_bound, initial_metadata, inputs, options, budget);
    const double search_seconds = chrono::duration<double>(chrono::steady_clock::now() - search_begin).count();

    const auto unique_candidates = rank_annotations(paths, ref, seqid, ref_fixed_score, options.k, options.raw);

    // The search interval is reported as a locus; each gene of each ranked
    // path becomes its own gene record inside it.
    string locus_id = "locus" + to_string(ref_idx + 1);
    if (!unique_candidates.empty()) {
        gff_out << ref.sequence_id << "\tUniAnn\tlocus\t"
            << L_bound + 1 << '\t' << R_bound + 1
            << "\t.\t+\t.\tID=" << locus_id << "\n";
    }

    for (const auto &cand : unique_candidates) {
        string g_id = locus_id + ".k" + to_string(cand.rank) +
                      ".g" + to_string(cand.gene_in_path);
        string t_id = g_id + ".t1";

        gff_out << cand.transcript.sequence_id << "\tUniAnn\tgene\t"
            << cand.transcript.genomic_start_0based + 1 << '\t'
            << cand.transcript.genomic_end_0based_exclusive
            << "\t.\t+\t.\tID=" << g_id
            << ";Parent=" << locus_id << "\n";
        write_hierarchical_transcript(gff_out, cand.transcript, g_id, t_id);

        tsv_out << (ref_idx + 1) << '\t'
                << cand.rank << '\t'
                << cand.genes_in_path << '\t'
                << cand.gene_in_path << '\t'
                << (ref.genomic_start_0based + 1) << '\t'
                << ref.genomic_end_0based_exclusive << '\t'
                << (cand.transcript.genomic_start_0based + 1) << '\t'
                << cand.transcript.genomic_end_0based_exclusive << '\t'
                << cand.score << '\t'
                << ref_fixed_score << '\t'
                << cand.score_delta << '\t'
                << cand.actual_start_score << '\t'
                << tsv_escape_chain(cand.transcript) << '\t'
                << tsv_escape_chain(ref) << '\t'
                << (cand.identical_to_reference ? "1" : "0") << '\t'
                << cand.classification << '\t' << L_bound + 1 << '\t'
                << R_bound + 1 << '\t' << initial_metadata.inter_len << '\t'
                << cand.splice_jaccard;
        if (!options.adaptive_k_spec.empty()) tsv_out << '\t' << options.k;
        tsv_out << '\n';
    }

    LocusOutput output;
    output.gff = gff_out.str();
    output.reference = reference_out.str();
    output.report = tsv_out.str();
    output.exact = exact_stats;
    output.search_seconds = search_seconds;
    output.window_length = R_bound - L_bound + 1;
    return output;
}
} // namespace

vector<LocalPathResult> search_paths(
    int left, int right, const PathMetadata &initial, const ModelInputs &inputs,
    const Options &options, const vector<uint8_t> *allowed_states
) {
    TracebackBudget budget(options.traceback_memory_bytes);
    return search_paths_impl(left, right, initial, inputs, options, budget, allowed_states);
}

bool run_local_k_best(
    const string &seqid,
    const vector<TranscriptAnnotation> &references,
    const Baseline &baseline,
    const ScoreSource &source,
    const Options &options
) {
    ostringstream gff_out, reference_out, tsv_out;
    gff_out << "##gff-version 3\n";
    reference_out << "##gff-version 3\n";
    tsv_out << "reference_gene_index\trank\tgenes_in_path\tgene_in_path\t"
            << "reference_start_1based\t"
            << "reference_end_1based\tcandidate_start_1based\t"
            << "candidate_end_1based\tfixed_interval_score\t"
            << "reference_fixed_interval_score\tscore_delta\t"
            << "actual_start_score\tintron_chain\t"
            << "reference_intron_chain\t"
            << "identical_to_reference\tclassification\twindow_start_1based\twindow_end_1based\tinitial_inter_len\t"
            << "splice_site_jaccard";
    if (!options.adaptive_k_spec.empty()) tsv_out << "\tsearch_k";
    tsv_out << '\n';

    // Each locus writes to independent buffers. Flush by reference index after
    // the parallel region, so scheduling never changes IDs or record order.
    vector<LocusOutput> outputs(references.size());
    TracebackBudget budget(options.traceback_memory_bytes);
    const int requested_threads = options.threads > 0 ? options.threads : omp_get_num_procs();
    const int thread_count = max(1, min(requested_threads, static_cast<int>(references.size())));
    // There cannot be more loci than sequence bases (an int-sized input).
    // A signed canonical loop is also accepted by older OpenMP front ends.
    const int locus_count = static_cast<int>(references.size());
    const auto k_rule = parse_adaptive_k_rule(options.adaptive_k_spec, options.k);
    vector<int> locus_k;
    locus_k.reserve(references.size());
    for (const auto &reference : references)
        locus_k.push_back(k_rule.budget(reference.cds_intervals.size()));
    #pragma omp parallel for schedule(dynamic) num_threads(thread_count)
    for (int ref_idx = 0; ref_idx < locus_count; ++ref_idx) {
        try {
            Options local_options = options;
            local_options.k = locus_k[ref_idx];
            outputs[ref_idx] = decode_locus(ref_idx, seqid, references,
                baseline, source, local_options, budget);
        }
        catch (const exception &error) { outputs[ref_idx].error = error.what(); }
        catch (...) { outputs[ref_idx].error = "Unknown decoder error"; }
    }
    // Exceptions cannot cross an OpenMP boundary. Report them on the caller.
    ExactSearchStats total;
    double search_seconds = 0, max_search_seconds = 0;
    size_t max_trace = 0, positions = 0;
    int longest_window = 0;
    for (const auto &result : outputs) {
        if (!result.error.empty()) throw runtime_error(result.error);
        total.entries += result.exact.entries;
        total.max_entries_row = max(total.max_entries_row, result.exact.max_entries_row);
        total.max_slots_row = max(total.max_slots_row, result.exact.max_slots_row);
        total.spilled_bytes += result.exact.spilled_bytes;
        max_trace = max(max_trace, result.exact.trace_bytes);
        search_seconds += result.search_seconds;
        max_search_seconds = max(max_search_seconds, result.search_seconds);
        positions += result.window_length;
        longest_window = max(longest_window, result.window_length);
        gff_out << result.gff;
        reference_out << result.reference;
        tsv_out << result.report;
    }
    cerr << setprecision(6) << "search_stats\tloci=" << references.size()
         << "\twindow_bases=" << positions << "\tlongest_window=" << longest_window
         << "\tsum_search_seconds=" << search_seconds << "\tmax_locus_seconds=" << max_search_seconds;
    if (options.metadata_states)
        cerr << "\tmean_entries_per_base=" << (positions ? double(total.entries) / positions : 0.0)
             << "\tmax_entries_row=" << total.max_entries_row << "\tmax_keys_row=" << total.max_slots_row
             << "\tmax_trace_MiB=" << max_trace / 1048576.0
             << "\tspilled_trace_MiB=" << total.spilled_bytes / 1048576.0;
    if (!options.adaptive_k_spec.empty()) {
        long double weighted = 0;
        long long sum_k = 0;
        for (size_t i = 0; i < locus_k.size(); ++i) {
            sum_k += locus_k[i];
            weighted += static_cast<long double>(locus_k[i]) * outputs[i].window_length;
        }
        cerr << "\tadaptive_min_k=" << (locus_k.empty() ? options.k : *min_element(locus_k.begin(), locus_k.end()))
             << "\tadaptive_max_k=" << (locus_k.empty() ? options.k : *max_element(locus_k.begin(), locus_k.end()))
             << "\tadaptive_mean_k=" << (locus_k.empty() ? 0.0 : double(sum_k) / locus_k.size())
             << "\tadaptive_window_weighted_k=" << (positions ? double(weighted / positions) : 0.0);
    }
    cerr << '\n';
    write_output_bundle({{options.output_filename, gff_out.str()},
                         {options.output_filename + ".reference.gff", reference_out.str()},
                         {options.report_filename, tsv_out.str()}});
    return true;
}
