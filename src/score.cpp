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
                   const Baseline &baseline, const ModelInputs &inputs,
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
    // Replace whole global coding runs, so no half-gene of the baseline remains.
    int replace_left = first, replace_right = last;
    while (replace_left > 1 && global[replace_left - 1] != 0) --replace_left;
    while (replace_right + 2 < L && global[replace_right + 1] != 0) ++replace_right;
    if (global[replace_left - 1] != 0 || global[replace_right + 1] != 0) {
        row << "unscorable\ttouches_sequence_end\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n";
        return row.str();
    }
    // Same window rule as decode_locus for a reference with this span.
    int left = replace_left - 1, right = replace_right + 1;
    const long long requested_left = static_cast<long long>(model.cds.front().start_0based) - options.flank;
    const long long requested_right = static_cast<long long>(model.cds.back().end_0based_exclusive) - 1 + options.flank;
    while (left > 0 && (options.flank < 0 || left > requested_left) && global[left - 1] == 0) --left;
    while (right + 1 < L && (options.flank < 0 || right < requested_right) && global[right + 1] == 0) ++right;

    PathMetadata initial = options.flank >= 0 ? baseline.noncoding_prefix[left] : PathMetadata{};
    if (left == 0 || options.flank < 0) {
        initial.predecessor_state = 0;
        initial.inter_len = 1;
    }
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
} // namespace

void run_model_scoring(const string &seqid, const vector<TranscriptAnnotation> &references,
                       const Baseline &baseline, const ModelInputs &inputs,
                       const Options &options) {
    const auto models = read_models(options.models_filename);
    vector<string> rows(models.size()), errors(models.size());
    const int requested = options.threads > 0 ? options.threads : omp_get_num_procs();
    const int threads = max(1, min(requested, static_cast<int>(max<size_t>(1, models.size()))));
    const int count = static_cast<int>(models.size());
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < count; ++i) {
        try { rows[i] = score_model(models[i], seqid, references, baseline, inputs, options); }
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
