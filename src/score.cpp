#include "altann.hpp"
#include <omp.h>

// Score externally supplied coding models under the same seven-state model.
// Each model replaces the global-best path over its own span, extended to
// whole global coding runs it touches, and both paths are replayed over one
// N-bounded interval. The difference is comparable to the K-best delta.
namespace {
struct ModelRecord { string id; vector<GenomicInterval> cds; };

struct Replay {
    bool allowed = true;
    double score = NEG_INF;
    int failed_position = -1;
    int failed_from = -1;
    int failed_to = -1;
};

const char *transition_kind(int from, int to) {
    if (from == 0 && is_exon(to)) return "start";
    if (is_exon(from) && is_intron(to)) return "donor";
    if (is_intron(from) && is_exon(to)) return "acceptor";
    if (is_exon(from) && to == 0) return "stop";
    if (is_exon(from) && is_exon(to)) return "internal_stop";
    return "state";
}

vector<ModelRecord> read_models(const string &path) {
    ifstream input(path);
    if (!input) throw runtime_error("Cannot open model table: " + path);
    vector<ModelRecord> models;
    set<string> ids;
    string line;
    while (getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto tab = line.find('\t');
        if (tab == string::npos || line.find('\t', tab + 1) != string::npos)
            throw runtime_error("Model rows must contain an ID and CDS intervals: " + line);
        ModelRecord model{line.substr(0, tab), {}};
        if (model.id.empty() || !ids.insert(model.id).second)
            throw runtime_error("Missing or duplicate model ID: " + model.id);
        istringstream intervals(line.substr(tab + 1));
        string item;
        while (getline(intervals, item, ',')) {
            const auto dash = item.find('-');
            size_t used_start = 0, used_end = 0;
            int start = -1, end = -1;
            try {
                if (dash == string::npos) throw invalid_argument("");
                start = stoi(item.substr(0, dash), &used_start);
                end = stoi(item.substr(dash + 1), &used_end);
            } catch (const exception &) {
                throw runtime_error("Invalid CDS interval for " + model.id + ": " + item);
            }
            if (used_start != dash || used_end != item.size() - dash - 1 || start < 0 || end <= start)
                throw runtime_error("Invalid CDS interval for " + model.id + ": " + item);
            model.cds.push_back({start, end});
        }
        if (model.cds.empty()) throw runtime_error("Model has no CDS: " + model.id);
        sort(model.cds.begin(), model.cds.end());
        models.push_back(model);
    }
    return models;
}

// Inverse of build_transcript_annotations: the N->E switch is two bases after
// the CDS start, a donor's G is the last exon state, an acceptor's G is the
// first exon state, and the stop codon's last base is N. Exon frames are the
// absolute codon phase; each intron keeps the preceding exon's frame.
string model_states(const vector<GenomicInterval> &cds, int length,
                    int &first, vector<int> &states) {
    for (size_t i = 0; i + 1 < cds.size(); ++i)
        if (cds[i + 1].start_0based - cds[i].end_0based_exclusive < 3)
            return "overlapping_or_adjacent_cds";
    first = cds.front().start_0based + 2;
    const int last = cds.back().end_0based_exclusive - 2;
    if (first > last) return "cds_too_short";
    if (first < 1 || last + 1 >= length) return "touches_sequence_end";
    states.assign(last - first + 1, 0);
    int frame = cds.front().start_0based % 3;
    for (size_t i = 0; i < cds.size(); ++i) {
        const int a = i == 0 ? first : cds[i].start_0based - 1;
        const int b = i + 1 == cds.size() ? last : cds[i].end_0based_exclusive;
        for (int p = a; p <= b; ++p) states[p - first] = 1 + frame;
        if (i + 1 < cds.size()) {
            for (int p = b + 1; p <= cds[i + 1].start_0based - 2; ++p)
                states[p - first] = 4 + frame;
            frame = (frame + cds[i + 1].start_0based - cds[i].end_0based_exclusive) % 3;
        }
    }
    return "";
}

struct Window {
    bool valid = false;
    int left = -1, right = -1, replace_left = -1, replace_right = -1;
    PathMetadata initial;
};

// Replace whole global coding runs touching [first, last], so no half-gene of
// the baseline remains, then extend through adjoining N runs as decode_locus
// does for a reference spanning [span_start, span_end).
Window model_window(int first, int last, int span_start, int span_end,
                    const Baseline &baseline, const Options &options) {
    const auto &global = baseline.states;
    const int L = global.size();
    Window w;
    w.replace_left = first;
    w.replace_right = last;
    while (w.replace_left > 1 && global[w.replace_left - 1] != 0) --w.replace_left;
    while (w.replace_right + 2 < L && global[w.replace_right + 1] != 0) ++w.replace_right;
    if (global[w.replace_left - 1] != 0 || global[w.replace_right + 1] != 0) return w;
    w.left = w.replace_left - 1;
    w.right = w.replace_right + 1;
    const long long requested_left = static_cast<long long>(span_start) - options.flank;
    const long long requested_right = static_cast<long long>(span_end) - 1 + options.flank;
    while (w.left > 0 && (options.flank < 0 || w.left > requested_left) && global[w.left - 1] == 0) --w.left;
    while (w.right + 1 < L && (options.flank < 0 || w.right < requested_right) && global[w.right + 1] == 0) ++w.right;
    w.initial = options.flank >= 0 ? baseline.noncoding_prefix[w.left] : PathMetadata{};
    if (w.left == 0 || options.flank < 0) {
        w.initial.predecessor_state = 0;
        w.initial.inter_len = 1;
    }
    w.valid = true;
    return w;
}

template <class StateAt>
Replay replay(int left, int right, PathMetadata metadata, StateAt state_at,
              const ModelInputs &inputs) {
    Replay result;
    result.score = inputs.emit[left][0];
    int current = state_at(left);
    for (int p = left + 1; p <= right; ++p) {
        const int next = state_at(p);
        const auto eval = evaluate_transition(p, current, next, metadata, inputs);
        if (!eval.allowed) {
            result = Replay{false, NEG_INF, p, current, next};
            return result;
        }
        result.score += eval.transition_score + eval.emission_score;
        metadata = eval.next_metadata;
        current = next;
    }
    return result;
}

string score_model(const ModelRecord &model, const string &seqid,
                   const vector<TranscriptAnnotation> &references,
                   const Baseline &baseline, const ScoreSource &source,
                   const Options &options) {
    const auto &global = baseline.states;
    const int L = global.size();
    ostringstream row;
    row << setprecision(17) << model.id << '\t';
    int first = -1;
    vector<int> states;
    const string error = model_states(model.cds, L, first, states);
    if (!error.empty()) {
        row << "unscorable\t" << error << "\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n";
        return row.str();
    }
    const int last = first + static_cast<int>(states.size()) - 1;
    const auto window = model_window(first, last, model.cds.front().start_0based,
                                     model.cds.back().end_0based_exclusive, baseline, options);
    if (!window.valid) {
        row << "unscorable\ttouches_sequence_end\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n";
        return row.str();
    }
    const int left = window.left, right = window.right;
    const int replace_left = window.replace_left, replace_right = window.replace_right;
    const PathMetadata &initial = window.initial;
    const auto scores = transition_window(source, left, right);
    const auto inputs = scores.inputs();
    const auto reference = replay(left, right, initial,
                                  [&](int p) { return global[p]; }, inputs);
    if (!reference.allowed)
        throw runtime_error("Forbidden baseline transition at position " +
                            to_string(reference.failed_position + 1));
    const auto candidate = replay(left, right, initial, [&](int p) {
        if (p < replace_left || p > replace_right) return global[p];
        return (p < first || p > last) ? 0 : states[p - first];
    }, inputs);

    // Relate the model to the global-best genes it replaces.
    TranscriptAnnotation annotation;
    {
        vector<int> positions, path;
        for (int p = left; p <= right; ++p) {
            positions.push_back(p);
            path.push_back((p < replace_left || p > replace_right) ? global[p]
                           : (p < first || p > last) ? 0 : states[p - first]);
        }
        for (const auto &t : build_transcript_annotations(positions, path, seqid))
            if (t.genomic_start_0based == model.cds.front().start_0based) annotation = t;
    }
    string indices;
    bool identical = false;
    double best_jaccard = 0.0;
    for (size_t i = 0; i < references.size(); ++i) {
        const auto &ref = references[i];
        if (ref.path_end < replace_left || ref.path_start > replace_right) continue;
        if (!indices.empty()) indices += ',';
        indices += to_string(i + 1);
        identical = identical || ref.cds_intervals == model.cds;
        if (!annotation.cds_intervals.empty())
            best_jaccard = max(best_jaccard, splice_site_jaccard(ref, annotation));
    }
    row << (candidate.allowed ? "ok\t.\t." : "forbidden\t");
    if (!candidate.allowed)
        row << transition_kind(candidate.failed_from, candidate.failed_to) << '\t'
            << candidate.failed_position + 1;
    row << '\t';
    if (candidate.allowed) row << candidate.score << '\t';
    else row << ".\t";
    row << reference.score << '\t';
    if (candidate.allowed) row << candidate.score - reference.score << '\t';
    else row << ".\t";
    row << left + 1 << '\t' << right + 1 << '\t' << replace_left + 1 << '\t'
        << replace_right + 1 << '\t' << (indices.empty() ? "." : indices) << '\t'
        << (identical ? 1 : 0) << '\t';
    if (indices.empty()) row << ".\n";
    else row << best_jaccard << '\n';
    return row.str();
}
// Decode the best path forced through every intron of a row, then score the
// gene carrying those introns exactly as a supplied model. Introns are
// zero-based half-open [GT, AG + 2), so the donor and acceptor G stay exon.
string force_row(const ModelRecord &row, const string &seqid,
                 const vector<TranscriptAnnotation> &references,
                 const Baseline &baseline, const ScoreSource &source,
                 const Options &options) {
    const auto &introns = row.cds;
    const int L = baseline.states.size();
    const auto fail = [&](const string &reason) {
        return row.id + "\t.\tunscorable\t" + reason + "\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n";
    };
    for (size_t i = 0; i < introns.size(); ++i) {
        if (introns[i].end_0based_exclusive - introns[i].start_0based < 4)
            return fail("intron_too_short");
        if (i + 1 < introns.size() && introns[i + 1].start_0based < introns[i].end_0based_exclusive + 3)
            return fail("overlapping_introns");
    }
    const int first = introns.front().start_0based;
    const int last = introns.back().end_0based_exclusive - 1;
    if (first < 2 || last + 2 >= L) return fail("touches_sequence_end");
    // Flank the replaced genes, not the introns, so their own starts and
    // stops stay reachable.
    const auto replaced = model_window(first, last, first, last + 1, baseline, options);
    if (!replaced.valid) return fail("touches_sequence_end");
    const auto window = model_window(first, last, replaced.replace_left,
                                     replaced.replace_right + 1, baseline, options);
    const auto scores = transition_window(source, window.left, window.right);
    const auto inputs = scores.inputs();
    const uint8_t exon = 0x0e, intron = 0x70;
    vector<uint8_t> allowed(window.right - window.left + 1, 0x7f);
    for (const auto &i : introns) {
        allowed[i.start_0based - window.left] = exon;
        allowed[i.end_0based_exclusive - 1 - window.left] = exon;
        for (int p = i.start_0based + 1; p <= i.end_0based_exclusive - 2; ++p)
            allowed[p - window.left] = intron;
    }
    Options single = options;
    single.k = 1;
    const auto paths = search_paths(window.left, window.right, window.initial, inputs, single, &allowed);
    if (paths.empty()) return fail("no_path");
    for (const auto &t : build_transcript_annotations(paths[0].path, seqid)) {
        const bool carries = all_of(introns.begin(), introns.end(), [&](const GenomicInterval &i) {
            return find(t.junctions.begin(), t.junctions.end(),
                        SpliceJunction{i.start_0based, i.end_0based_exclusive - 2, '+'}) != t.junctions.end();
        });
        if (!carries) continue;
        string cds;
        for (const auto &c : t.cds_intervals)
            cds += (cds.empty() ? "" : ",") + to_string(c.start_0based) + "-" + to_string(c.end_0based_exclusive);
        const string scored = score_model({row.id, t.cds_intervals}, seqid, references, baseline, source, options);
        return row.id + "\t" + cds + scored.substr(row.id.size());
    }
    return fail("no_complete_gene");
}
} // namespace

void run_junction_forcing(const string &seqid, const vector<TranscriptAnnotation> &references,
                          const Baseline &baseline, const ScoreSource &source,
                          const Options &options) {
    const auto rows = read_models(options.junctions_filename);
    vector<string> out_rows(rows.size()), errors(rows.size());
    const int requested = options.threads > 0 ? options.threads : omp_get_num_procs();
    const int threads = max(1, min(requested, static_cast<int>(max<size_t>(1, rows.size()))));
    // Each worker owns a constrained search budget; bound their combined total.
    Options worker_options = options;
    worker_options.traceback_memory_bytes /= threads;
    const int count = static_cast<int>(rows.size());
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < count; ++i) {
        try { out_rows[i] = force_row(rows[i], seqid, references, baseline, source, worker_options); }
        catch (const exception &error) { errors[i] = error.what(); }
        catch (...) { errors[i] = "Unknown forcing error"; }
    }
    ostringstream out;
    out << "model_id\tcds\tstatus\treason\tfailed_position_1based\tmodel_score\t"
        << "reference_score\tdelta\twindow_start_1based\twindow_end_1based\t"
        << "replaced_start_1based\treplaced_end_1based\treference_gene_indices\t"
        << "identical_to_reference\tsplice_site_jaccard\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!errors[i].empty()) throw runtime_error(rows[i].id + ": " + errors[i]);
        out << out_rows[i];
    }
    write_output_bundle({{options.output_filename, out.str()}});
}

void run_model_scoring(const string &seqid, const vector<TranscriptAnnotation> &references,
                       const Baseline &baseline, const ScoreSource &source,
                       const Options &options) {
    const auto models = read_models(options.models_filename);
    vector<string> rows(models.size()), errors(models.size());
    const int requested = options.threads > 0 ? options.threads : omp_get_num_procs();
    const int threads = max(1, min(requested, static_cast<int>(max<size_t>(1, models.size()))));
    const int count = static_cast<int>(models.size());
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < count; ++i) {
        try { rows[i] = score_model(models[i], seqid, references, baseline, source, options); }
        catch (const exception &error) { errors[i] = error.what(); }
        catch (...) { errors[i] = "Unknown scoring error"; }
    }
    ostringstream out;
    out << "model_id\tstatus\treason\tfailed_position_1based\tmodel_score\t"
        << "reference_score\tdelta\twindow_start_1based\twindow_end_1based\t"
        << "replaced_start_1based\treplaced_end_1based\treference_gene_indices\t"
        << "identical_to_reference\tsplice_site_jaccard\n";
    for (size_t i = 0; i < models.size(); ++i) {
        if (!errors[i].empty()) throw runtime_error(models[i].id + ": " + errors[i]);
        out << rows[i];
    }
    write_output_bundle({{options.output_filename, out.str()}});
}
