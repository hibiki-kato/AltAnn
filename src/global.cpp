// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"

Baseline decode_baseline(const ModelInputs &inputs) {
    const int length = inputs.seq.size();
    struct Cell { double score; PathMetadata metadata; };
    array<Cell, NUM_STATES> previous{}, current{};
    vector<array<int8_t, NUM_STATES>> traceback(length);
    Baseline result;
    result.noncoding_prefix.resize(length);
    result.states.resize(length);
    for (int state = 0; state < NUM_STATES; ++state) {
        previous[state].score = (state == 0 ? 0.0 : NEG_INF) + inputs.emit[0][state];
        previous[state].metadata = {is_intron(state) ? 1 : 0,
                                   is_exon(state) ? 1 : 0, state == 0 ? 1 : 0, 0, -1};
        traceback[0][state] = -1;
    }
    result.noncoding_prefix[0] = previous[0].metadata;
    // Only seven scores are needed at once. One signed byte per state retains
    // the traceback; only N prefix metadata is needed by local decoding.
    for (int position = 1; position < length; ++position) {
        for (int to = 0; to < NUM_STATES; ++to) {
            double best = -1e18;
            int best_from = -1;
            PathMetadata metadata;
            for (int from = 0; from < NUM_STATES; ++from) {
                const auto transition = evaluate_transition(position, from, to,
                                                           previous[from].metadata, inputs);
                // Preserve finite unreachable-cell arithmetic and strict tie
                // order from UniAnn. Local decoding checks allowed separately.
                const double score = previous[from].score + transition.transition_score
                                     + transition.emission_score;
                if (score > best) {
                    best = score;
                    best_from = from;
                    metadata = transition.next_metadata;
                }
            }
            current[to] = {best, metadata};
            traceback[position][to] = static_cast<int8_t>(best_from);
        }
        previous.swap(current);
        result.noncoding_prefix[position] = previous[0].metadata;
    }
    int state = 0;
    for (int candidate = 1; candidate < NUM_STATES; ++candidate)
        if (previous[candidate].score > previous[state].score) state = candidate;
    for (int position = length - 1; position >= 0; --position) {
        if (state < 0 || state >= NUM_STATES) throw runtime_error("Invalid baseline traceback");
        result.states[position] = state;
        state = traceback[position][state];
    }
    return result;
}
