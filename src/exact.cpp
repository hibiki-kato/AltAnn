// Exact local K-best decoding over (coarse state, canonical metadata) keys.
//
// The default search keeps K paths per coarse state, but transitions also read
// PathMetadata (exon/intron/intergenic lengths, the exon's origin, the
// unfinished spliced codon). Two paths in the same coarse state can therefore
// have different futures, and a pruned path could have led to a better
// complete path. Here each hidden state is a coarse state plus every metadata
// distinction a future transition can observe, so keeping the top K partial
// paths per key is the exact K-best list Viterbi for this model.
#include "altann.hpp"
#include <unistd.h>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <unordered_map>

namespace {

// Spliced-codon suffix classes. Only equality with TAA/TAG/TGA is ever tested,
// so a pending suffix is equivalent to any other with the same stop prefix.
int suffix_class(const PathMetadata &m) {
    if (m.coding_suffix_length == 0) return 0;
    if (m.coding_suffix_length == 1) return m.coding_suffix[0] == 'T' ? 1 : 2;
    if (m.coding_suffix[0] == 'T' && m.coding_suffix[1] == 'A') return 3;
    if (m.coding_suffix[0] == 'T' && m.coding_suffix[1] == 'G') return 4;
    return 5;
}

void set_suffix_class(PathMetadata &m, int cls) {
    static const char first[6]  = {0, 'T', 'C', 'T', 'T', 'C'};
    static const char second[6] = {0, 0, 0, 'A', 'G', 'C'};
    static const int length[6]  = {0, 1, 1, 2, 2, 2};
    m.coding_suffix = {first[cls], second[cls]};
    m.coding_suffix_length = length[cls];
}

constexpr int N_INTER_CAP = MIN_INTER - 2;          // inter_len + 2 >= MIN_INTER
constexpr int FIRST_EXON_CAP = MIN_SINGLE + 3;      // exon_len - 2 > MIN_SINGLE
constexpr int LATER_EXON_CAP = MIN_EXON + 2;        // exon_len - 2 >= MIN_EXON
constexpr int INTRON_THRESHOLD = MIN_INTRON - 2;    // intron_len + 2 >= MIN_INTRON

// Canonical metadata: the representative of every raw metadata that yields the
// same allowed transitions, scores and canonical successors.
PathMetadata canonicalize(int state, const PathMetadata &m) {
    PathMetadata key;
    if (state == 0) {
        key.inter_len = min(m.inter_len, N_INTER_CAP);
        key.predecessor_state = m.predecessor_state == 0 ? 0 : 1;
        return key;
    }
    key.coding_started = m.coding_started;
    set_suffix_class(key, suffix_class(m));
    if (is_exon(state)) {
        key.exon_from = m.exon_from == 0 ? 0 : 4;
        key.exon_len = min(m.exon_len, key.exon_from == 0 ? FIRST_EXON_CAP : LATER_EXON_CAP);
    } else {
        key.intron_len = m.intron_len < INTRON_THRESHOLD
            ? m.intron_len : INTRON_THRESHOLD + (m.intron_len - INTRON_THRESHOLD) % 3;
    }
    return key;
}

// Dense key ids. Coding fields (suffix class 0..5, started flag) multiply the
// exon and intron blocks.
constexpr int N_WIDTH = (N_INTER_CAP + 1) * 2;
constexpr int EXON_LEN_WIDTH = (FIRST_EXON_CAP + 1) + (LATER_EXON_CAP + 1);
constexpr int CODING_WIDTH = 6 * 2;
constexpr int EXON_WIDTH = EXON_LEN_WIDTH * CODING_WIDTH;
constexpr int INTRON_LEN_WIDTH = INTRON_THRESHOLD + 3;
constexpr int INTRON_WIDTH = INTRON_LEN_WIDTH * CODING_WIDTH;
constexpr int EXON_OFFSET = N_WIDTH;
constexpr int INTRON_OFFSET = EXON_OFFSET + 3 * EXON_WIDTH;
constexpr int KEY_COUNT = INTRON_OFFSET + 3 * INTRON_WIDTH;

int key_id(int state, const PathMetadata &key) {
    if (state == 0) {
        const int inter = max(0, key.inter_len);
        return inter * 2 + key.predecessor_state;
    }
    const int coding = suffix_class(key) * 2 + (key.coding_started ? 1 : 0);
    if (is_exon(state)) {
        const int length = key.exon_from == 0 ? key.exon_len : FIRST_EXON_CAP + 1 + key.exon_len;
        return EXON_OFFSET + (state - 1) * EXON_WIDTH + length * CODING_WIDTH + coding;
    }
    return INTRON_OFFSET + (state - 4) * INTRON_WIDTH + key.intron_len * CODING_WIDTH + coding;
}

struct Slot {
    int state = -1;
    PathMetadata metadata;
    // Indices into the row's entry array, in descending score order.
    int first = 0;
    int count = 0;
};

struct Candidate {
    double score;
    int source_slot;   // position in the previous row's active list
    int source_rank;
    bool operator<(const Candidate &o) const {
        if (score != o.score) return score > o.score;
        if (source_slot != o.source_slot) return source_slot < o.source_slot;
        return source_rank < o.source_rank;
    }
};

// Backpointer entry: its parent entry in the previous row and its coarse state.
struct TraceEntry {
    // Seven coarse states need three bits; the remaining 29 bits address the
    // previous row. Keep four-byte entries while supporting K=100 and beyond.
    uint32_t packed;
    TraceEntry(uint32_t parent, int state) : packed((parent << 3) | uint32_t(state)) {}
    explicit TraceEntry(uint32_t raw) : packed(raw) {}
    uint32_t parent() const { return packed >> 3; }
    int state() const { return packed & 7u; }
};
static_assert(sizeof(TraceEntry) == 4, "Exact traceback entries must stay compact");

// Exact rows hold a data-dependent number of entries, so the traceback is one
// word stream: each row's entries, then its entry count. Walking it backward
// needs no per-row index. Full blocks stay in memory while the shared budget
// lasts; later blocks go to an unlinked temporary file and are read back into a
// two-block cache. The open block is not counted against the budget.
class TraceStream {
public:
    TraceStream(TracebackBudget &budget, string directory, size_t block_words = size_t(1) << 20)
        : budget_(budget), directory_(std::move(directory)), block_words_(block_words) {}
    TraceStream(const TraceStream &) = delete;
    TraceStream &operator=(const TraceStream &) = delete;
    ~TraceStream() {
        if (descriptor_ >= 0) close(descriptor_);
        budget_.release(reserved_);
    }
    void push(uint32_t word) {
        tail_.push_back(word);
        if (tail_.size() == block_words_) seal();
    }
    size_t size() const { return blocks_.size() * block_words_ + tail_.size(); }
    size_t spilled_bytes() const { return spilled_ * block_words_ * sizeof(uint32_t); }
    uint32_t operator[](size_t index) {
        const size_t block = index / block_words_, offset = index % block_words_;
        if (block == blocks_.size()) return tail_[offset];
        if (!blocks_[block].empty()) return blocks_[block][offset];
        for (int c = 0; c < 2; ++c)
            if (cached_[c] == block) { recent_ = c; return cache_[c][offset]; }
        const int slot = 1 - recent_;
        read_block(block, cache_[slot]);
        cached_[slot] = block;
        recent_ = slot;
        return cache_[slot][offset];
    }
private:
    static constexpr size_t NONE = numeric_limits<size_t>::max();
    void seal() {
        const size_t bytes = block_words_ * sizeof(uint32_t);
        if (budget_.try_reserve(bytes)) {
            reserved_ += bytes;
            blocks_.push_back(std::move(tail_));
            tail_ = vector<uint32_t>();
            return;
        }
        if (descriptor_ < 0) {
            string path = (directory_.empty() ? string(".") : directory_) + "/.altann-exact-trace-XXXXXX";
            descriptor_ = mkstemp(path.data());
            if (descriptor_ < 0)
                throw runtime_error("Cannot create traceback file in " + directory_ + ": " + strerror(errno));
            unlink(path.c_str());
        }
        const char *cursor = reinterpret_cast<const char *>(tail_.data());
        size_t remaining = bytes;
        off_t offset = static_cast<off_t>(blocks_.size() * bytes);
        while (remaining) {
            const ssize_t written = pwrite(descriptor_, cursor, remaining, offset);
            if (written < 0 && errno == EINTR) continue;
            if (written <= 0) throw runtime_error(string("Cannot write traceback file: ") + strerror(errno));
            cursor += written;
            remaining -= written;
            offset += written;
        }
        blocks_.emplace_back();  // empty: the block is on disk
        ++spilled_;
        tail_.clear();
    }
    void read_block(size_t block, vector<uint32_t> &buffer) const {
        buffer.resize(block_words_);
        char *cursor = reinterpret_cast<char *>(buffer.data());
        size_t remaining = block_words_ * sizeof(uint32_t);
        off_t offset = static_cast<off_t>(block * remaining);
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
    string directory_;
    size_t block_words_;
    size_t reserved_ = 0, spilled_ = 0;
    int descriptor_ = -1;
    vector<vector<uint32_t>> blocks_;
    vector<uint32_t> tail_;
    array<vector<uint32_t>, 2> cache_;
    array<size_t, 2> cached_{NONE, NONE};
    int recent_ = 0;
};

// Exact K-best on [L_bound, R_bound] with both ends in N. input_at(p) returns
// ModelInputs valid for a transition into p; positions are visited in
// ascending order except for start scores, read only when want_start_scores.
template <class InputAt>
vector<LocalPathResult> exact_search(
    int L_bound, int R_bound, const PathMetadata &initial_metadata, InputAt &&input_at,
    const Options &options, TracebackBudget &budget, ExactSearchStats *stats, bool want_start_scores
) {
    const int K = options.k;
    const int local_len = R_bound - L_bound + 1;
    if (static_cast<long long>(KEY_COUNT) * K >= (1LL << 29))
        throw runtime_error("K too large for exact metadata-state search");

    vector<Slot> prev_slots, curr_slots;
    vector<double> prev_scores, curr_scores;
    vector<int> slot_of_key(KEY_COUNT, -1);
    TraceStream trace(budget, spill_directory(options));

    {
        Slot s;
        s.state = 0;
        // Transitions are invariant to everything the canonical form drops.
        s.metadata = canonicalize(0, initial_metadata);
        s.first = 0;
        s.count = 1;
        prev_slots.push_back(s);
        prev_scores.push_back(0.0 + input_at(L_bound).emit[L_bound][0]);
        trace.push(TraceEntry(0, 0).packed);
        trace.push(1);
    }

    struct Pending { int state; PathMetadata metadata; vector<Candidate> candidates; };
    vector<Pending> pending;   // pool reused across positions; first pending_used are live
    size_t pending_used = 0;
    vector<int> touched, order;
    size_t max_entries_row = 1, sum_entries = 1;
    size_t max_slots_row = 1;

    for (int i = 1; i < local_len; i++) {
        const int p = L_bound + i;
        const ModelInputs inputs = input_at(p);
        pending_used = 0;
        touched.clear();
        for (int s = 0; s < static_cast<int>(prev_slots.size()); ++s) {
            const Slot &source = prev_slots[s];
            const int from = source.state;
            int targets[4];
            int target_count = 0;
            targets[target_count++] = from;
            if (from == 0) { targets[target_count++] = 1; targets[target_count++] = 2; targets[target_count++] = 3; }
            else if (is_exon(from)) { targets[target_count++] = from + 3; targets[target_count++] = 0; }
            else { targets[target_count++] = 1; targets[target_count++] = 2; targets[target_count++] = 3; }
            for (int t = 0; t < target_count; ++t) {
                const int to = targets[t];
                if (options.flank >= 0 && from == 0 && is_exon(to) &&
                    codon_start_0based(p) < L_bound) continue;
                const TransitionResult edge = evaluate_transition(p, from, to, source.metadata, inputs);
                if (!edge.allowed) continue;
                const PathMetadata key = canonicalize(to, edge.next_metadata);
                const int id = key_id(to, key);
                int &index = slot_of_key[id];
                if (index < 0) {
                    index = static_cast<int>(pending_used++);
                    if (pending.size() < pending_used) pending.emplace_back();
                    pending[index].state = to;
                    pending[index].metadata = key;
                    pending[index].candidates.clear();
                    touched.push_back(id);
                }
                auto &cands = pending[index].candidates;
                // Same summation order as the coarse search.
                for (int r = 0; r < source.count; ++r)
                    cands.push_back({prev_scores[source.first + r] + edge.transition_score +
                                     edge.emission_score, s, r});
            }
        }
        // Deterministic slot order: by key id.
        order.resize(touched.size());
        for (size_t j = 0; j < order.size(); ++j) order[j] = static_cast<int>(j);
        sort(order.begin(), order.end(), [&](int a, int b) { return touched[a] < touched[b]; });

        curr_slots.clear();
        curr_scores.clear();
        for (int j : order) {
            Pending &pend = pending[slot_of_key[touched[j]]];
            auto &cands = pend.candidates;
            const int keep = min(K, static_cast<int>(cands.size()));
            partial_sort(cands.begin(), cands.begin() + keep, cands.end());
            Slot slot;
            slot.state = pend.state;
            slot.metadata = pend.metadata;
            slot.first = static_cast<int>(curr_scores.size());
            slot.count = keep;
            for (int r = 0; r < keep; ++r) {
                curr_scores.push_back(cands[r].score);
                const int parent = prev_slots[cands[r].source_slot].first + cands[r].source_rank;
                trace.push(TraceEntry(static_cast<uint32_t>(parent), pend.state).packed);
            }
            curr_slots.push_back(slot);
        }
        trace.push(static_cast<uint32_t>(curr_scores.size()));
        for (int id : touched) slot_of_key[id] = -1;
        max_entries_row = max(max_entries_row, curr_scores.size());
        max_slots_row = max(max_slots_row, curr_slots.size());
        sum_entries += curr_scores.size();
        prev_slots.swap(curr_slots);
        prev_scores.swap(curr_scores);
    }

    if (stats) {
        stats->positions += local_len;
        stats->entries += sum_entries;
        stats->max_entries_row = max(stats->max_entries_row, max_entries_row);
        stats->max_slots_row = max(stats->max_slots_row, max_slots_row);
        stats->trace_bytes = trace.size() * sizeof(uint32_t);
        stats->spilled_bytes += trace.spilled_bytes();
    }

    // Final paths end in N; merge every N key's list.
    vector<Candidate> finals;
    for (int s = 0; s < static_cast<int>(prev_slots.size()); ++s) {
        if (prev_slots[s].state != 0) continue;
        for (int r = 0; r < prev_slots[s].count; ++r)
            finals.push_back({prev_scores[prev_slots[s].first + r], s, r});
    }
    sort(finals.begin(), finals.end());
    if (static_cast<int>(finals.size()) > K) finals.resize(K);

    // One backward pass over the stream traces every final path.
    const size_t n = finals.size();
    vector<uint32_t> entry(n);
    vector<int> state(n, 0);
    vector<vector<int>> starts(n), start_positions(n);
    vector<vector<int8_t>> states(n);
    for (size_t f = 0; f < n; ++f)
        entry[f] = static_cast<uint32_t>(prev_slots[finals[f].source_slot].first + finals[f].source_rank);
    size_t row_start = n ? trace.size() - 1 - trace[trace.size() - 1] : 0;
    for (int i = local_len - 1; n && i > 0; --i) {
        const size_t prev_start = row_start - 1 - trace[row_start - 1];
        const int p = L_bound + i;
        for (size_t f = 0; f < n; ++f) {
            const TraceEntry te(trace[row_start + entry[f]]);
            const int s = te.state();
            const TraceEntry pe(trace[prev_start + te.parent()]);
            if (pe.state() == 0 && is_exon(s)) start_positions[f].push_back(p);
            if (pe.state() != s) {
                starts[f].push_back(p);
                states[f].push_back(static_cast<int8_t>(s));
            }
            entry[f] = te.parent();
            state[f] = pe.state();
        }
        row_start = prev_start;
    }

    vector<LocalPathResult> paths;
    for (size_t f = 0; f < n; ++f) {
        LocalPathResult result;
        result.score = finals[f].score;
        starts[f].push_back(L_bound);
        states[f].push_back(static_cast<int8_t>(state[f]));
        reverse(starts[f].begin(), starts[f].end());
        reverse(states[f].begin(), states[f].end());
        if (want_start_scores)
            for (auto it = start_positions[f].rbegin(); it != start_positions[f].rend(); ++it)
                result.start_scores.push_back(input_at(*it).atg_score[codon_start_0based(*it)]);
        result.path = StatePath(std::move(starts[f]), std::move(states[f]), local_len);

        if (options.verify_exact) {
            // Replay with raw metadata: the score must match and every step be allowed.
            PathMetadata meta = initial_metadata;
            double score = input_at(L_bound).emit[L_bound][0];
            int cur = 0;
            for (int i = 1; i < local_len; ++i) {
                const int p = L_bound + i;
                const int next = result.path[p];
                const auto ev = evaluate_transition(p, cur, next, meta, input_at(p));
                if (!ev.allowed) throw runtime_error("Exact path replay hit a forbidden transition at " + to_string(p + 1));
                score += ev.transition_score + ev.emission_score;
                meta = ev.next_metadata;
                cur = next;
            }
            if (fabs(score - result.score) > 1e-6 * max(1.0, fabs(score)))
                throw runtime_error("Exact path replay score mismatch at window " + to_string(L_bound + 1) + ": " + to_string(score) + " vs " + to_string(result.score));
        }
        paths.push_back(std::move(result));
    }
    return paths;
}

} // namespace

string spill_directory(const Options &options) {
    if (!options.spill_dir.empty()) return options.spill_dir;
    return filesystem::path(options.output_filename).parent_path().string();
}

vector<LocalPathResult> search_paths_exact(
    int L_bound, int R_bound, const PathMetadata &initial_metadata,
    const ModelInputs &inputs, const Options &options, TracebackBudget &budget, ExactSearchStats *stats
) {
    return exact_search(L_bound, R_bound, initial_metadata,
                        [&](int) -> const ModelInputs & { return inputs; }, options, budget, stats, true);
}

// Whole-sequence exact best path: K=1 over (coarse state, canonical metadata)
// keys, the exact counterpart of decode_baseline. It starts from decode_baseline's
// position-0 N cell and, like local search, ends in N. Scores are streamed through
// a ScoreCursor and the traceback spills past --traceback-memory, so memory does
// not grow with the sequence. The noncoding prefix holds this path's own N
// metadata, which is where every local window starts.
Baseline decode_baseline_exact(const ScoreSource &source, const Options &options) {
    const int length = source.length();
    PathMetadata initial;
    initial.inter_len = 1;
    Options single = options;
    single.k = 1;
    single.flank = max(single.flank, 0);  // keep the start-codon guard at position 0
    single.verify_exact = false;           // replayed below in any case
    TracebackBudget budget(options.traceback_memory_bytes);
    ExactSearchStats stats;
    vector<LocalPathResult> paths;
    {
        ScoreCursor cursor(source);
        paths = exact_search(0, length - 1, initial, [&](int p) { return cursor.at(p); },
                             single, budget, &stats, false);
    }
    if (paths.empty()) throw runtime_error("Exact baseline found no path ending in N");
    cerr << "exact_baseline\tpositions=" << length << "\ttrace_MiB=" << stats.trace_bytes / 1048576.0
         << "\tspilled_MiB=" << stats.spilled_bytes / 1048576.0 << '\n';
    Baseline result;
    result.states = std::move(paths[0].path);
    ScoreCursor cursor(source);
    PathMetadata meta = initial;
    double score = cursor.at(0).emit[0][0];
    int current = 0;
    result.noncoding_prefix.append(0, meta);
    for (int p = 1; p < length; ++p) {
        const int next = result.states[p];
        const auto edge = evaluate_transition(p, current, next, meta, cursor.at(p));
        if (!edge.allowed)
            throw runtime_error("Exact baseline replay hit a forbidden transition at " + to_string(p + 1));
        score += edge.transition_score + edge.emission_score;
        meta = edge.next_metadata;
        current = next;
        if (current == 0) result.noncoding_prefix.append(p, meta);
    }
    if (fabs(score - paths[0].score) > 1e-6 * max(1.0, fabs(score)))
        throw runtime_error("Exact baseline replay score mismatch: " + to_string(score) + " vs " +
                            to_string(paths[0].score));
    result.final_scores.fill(NEG_INF);
    result.final_scores[0] = score;
    return result;
}
