// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"

Baseline load_baseline_log(const ModelInputs &inputs, const string &path) {
    ifstream log(path);
    if (!log) throw runtime_error("Cannot open Viterbi log: " + path);
    const int length = inputs.seq.size();
    struct Cell { double score; PathMetadata metadata; };
    array<Cell, NUM_STATES> previous{}, current{};
    array<int, NUM_STATES> logged_scores{};
    vector<array<int8_t, NUM_STATES>> traceback(length);
    Baseline result;
    result.noncoding_prefix.resize(length);
    result.states.resize(length);
    int expected_position = 0;
    bool expect_traceback = false;
    size_t line_number = 0;
    const auto fail = [&](const string &reason) {
        throw runtime_error("Invalid Viterbi log " + path + " at line " +
                            to_string(line_number) + ": " + reason);
    };
    const auto parse_integer = [&](const string &token) {
        size_t end = 0;
        int value = 0;
        try { value = stoi(token, &end); }
        catch (const exception &) { fail("expected an integer, got '" + token + "'"); }
        if (end != token.size()) fail("expected an integer, got '" + token + "'");
        return value;
    };
    string line;
    while (getline(log, line)) {
        ++line_number;
        istringstream row(line);
        string position_token, kind;
        if (!(row >> position_token >> kind) || (kind != "dp" && kind != "bt"))
            continue; // UniAnn also writes warnings and progress to stderr.
        const int position = parse_integer(position_token);
        if (position != expected_position || position >= length)
            fail("expected position " + to_string(expected_position) +
                 ", got " + to_string(position));
        if (kind != (expect_traceback ? "bt" : "dp"))
            fail("expected " + string(expect_traceback ? "bt" : "dp") +
                 " row at position " + to_string(position));
        array<int, NUM_STATES> values{};
        for (int state = 0; state < NUM_STATES; ++state) {
            string token;
            if (!(row >> token)) fail("each dp/bt row must contain seven values");
            values[state] = parse_integer(token);
        }
        string extra;
        if (row >> extra) fail("each dp/bt row must contain exactly seven values");
        if (kind == "dp") {
            logged_scores = values;
            expect_traceback = true;
            continue;
        }
        for (int state = 0; state < NUM_STATES; ++state) {
            const int from = values[state];
            if (position == 0) {
                if (from != -1) fail("initial traceback values must be -1");
                current[state].score = (state == 0 ? 0.0 : NEG_INF) + inputs.emit[0][state];
                current[state].metadata = {is_intron(state) ? 1 : 0,
                    is_exon(state) ? 1 : 0, state == 0 ? 1 : 0, 0, -1};
            } else {
                if (from < 0 || from >= NUM_STATES)
                    fail("traceback predecessor must be between 0 and 6");
                const auto transition = evaluate_transition(position, from, state,
                                                           previous[from].metadata, inputs);
                // Replay the supplied predecessor, without searching all seven
                // alternatives. Finite unreachable-cell arithmetic must match
                // UniAnn; transition.allowed is only used by local decoding.
                current[state].score = previous[from].score + transition.transition_score
                                       + transition.emission_score;
                current[state].metadata = transition.next_metadata;
            }
            // Upstream prints int(dp), truncating toward zero. Recover exact
            // scores from the model, then compare their printed representation.
            // Rounded log values must never become accumulated path scores.
            const double printed = trunc(current[state].score);
            if (!isfinite(printed) || printed < numeric_limits<int>::min() ||
                printed > numeric_limits<int>::max())
                fail("reconstructed score cannot be represented by UniAnn's integer log");
            if (static_cast<int>(printed) != logged_scores[state])
                fail("score mismatch at position " + to_string(position) +
                     ", state " + state_name[state] +
                     "; check the sequence, scores, and UniAnn model version");
            traceback[position][state] = static_cast<int8_t>(from);
        }
        previous.swap(current);
        result.noncoding_prefix[position] = previous[0].metadata;
        ++expected_position;
        expect_traceback = false;
    }
    if (log.bad()) fail("could not finish reading the file");
    if (expected_position != length || expect_traceback)
        fail("incomplete dp/bt rows; expected " + to_string(length) + " positions");
    int state = 0;
    for (int candidate = 1; candidate < NUM_STATES; ++candidate)
        if (previous[candidate].score > previous[state].score) state = candidate;
    for (int position = length - 1; position >= 0; --position) {
        result.states[position] = state;
        state = traceback[position][state];
    }
    return result;
}

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
