// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"
#include <omp.h>

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

struct LocalPathResult {
    vector<int> path_positions;
    vector<int> path_states;
    double score;
    // One ATG score per gene on the path, in genomic order. A local path may
    // contain zero, one or several genes.
    vector<double> start_scores;
};


struct LocalKBestCandidate {
    int rank;
    double score;        // score of the whole local path this gene belongs to
    double score_delta;
    double actual_start_score;
    int genes_in_path = 1;
    int gene_in_path = 1;
    string classification;
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

struct LocusOutput { string gff, reference, report, error; };

// Decode the bounded interval with both endpoint states pinned to N. This
// function owns all local DP memory and returns paths in descending score order.
vector<LocalPathResult> search_paths(
    int L_bound, int R_bound, const PathMetadata &initial_metadata,
    const ModelInputs &inputs, const Options &options
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
    if (static_cast<size_t>(local_len) > vector<CompactBP>().max_size() / cells_per_position)
        throw runtime_error("Local traceback exceeds the supported allocation size");
    vector<CompactBP> backpointers(static_cast<size_t>(local_len) * cells_per_position,
                                   CompactBP{-1, -1});

    prev_dp[0][0].dp = 0.0 + inputs.emit[L_bound][0];
    prev_dp[0][0].metadata = initial_metadata;

    for (int i = 1; i < local_len; i++) {
        int p = L_bound + i;

        for (int to = 0; to < NUM_STATES; to++) {
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
                backpointers[static_cast<size_t>(i) * NUM_STATES * K + to * K + k] =
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

    vector<LocalPathResult> paths;
    for (int k = 0; k < K; k++) {
        if (prev_dp[0][k].dp <= NEG_INF) continue;

        LocalPathResult res;
        res.score = prev_dp[0][k].dp;

        int cur_state = 0;
        int cur_k = k;
        bool valid = true;

        for (int i = local_len - 1; i >= 0; i--) {
            int p = L_bound + i;
            res.path_positions.push_back(p);
            res.path_states.push_back(cur_state);

            if (i > 0) {
                CompactBP bp = backpointers[static_cast<size_t>(i) * NUM_STATES * K +
                                            cur_state * K + cur_k];
                int next_state = bp.bt_state;
                int next_k = bp.bt_k;
                if (next_state < 0) { valid = false; break; }

                if (next_state == 0 && is_exon(cur_state)) {
                    res.start_scores.push_back(inputs.atg_score[codon_start_0based(p)]);
                }

                cur_state = next_state;
                cur_k = next_k;
            }
        }

        if (valid) {
            reverse(res.path_positions.begin(), res.path_positions.end());
            reverse(res.path_states.begin(), res.path_states.end());
            reverse(res.start_scores.begin(), res.start_scores.end());
            paths.push_back(res);
        }
    }

    return paths;
}

// Distinct state paths may describe the same complete annotation. Collapse
// those duplicates while retaining score order and every gene on each path.
vector<LocalKBestCandidate> rank_annotations(
    const vector<LocalPathResult> &paths, const TranscriptAnnotation &ref,
    const string &seqid, double ref_fixed_score, int K
) {
    vector<LocalKBestCandidate> unique_candidates;
    set<string> seen;

    // A local path is one candidate. Paths that annotate nothing (all N) are
    // dropped; paths carrying several genes are kept and reported as several
    // gene records under the same rank.
    int rank_count = 0;
    for (const auto &path : paths) {
        auto anns = build_transcript_annotations(path.path_positions, path.path_states, seqid);
        if (anns.empty()) continue;

        string key;
        for (const auto &ann : anns) key += full_transcript_key(ann) + "|";
        if (!seen.insert(key).second) continue;

        ++rank_count;
        const string ref_key = full_transcript_key(ref);
        for (size_t j = 0; j < anns.size(); ++j) {
            TranscriptAnnotation ann = anns[j];
            ann.score = path.score;

            LocalKBestCandidate cand;
            cand.rank = rank_count;
            cand.transcript = ann;
            cand.score = path.score;
            cand.score_delta = path.score - ref_fixed_score;
            cand.actual_start_score = j < path.start_scores.size()
                ? path.start_scores[j] : NEG_INF;
            cand.genes_in_path = static_cast<int>(anns.size());
            cand.gene_in_path = static_cast<int>(j) + 1;
            cand.classification = classify_splice_difference(ref, ann);
            cand.identical_to_reference =
                (anns.size() == 1 && full_transcript_key(ann) == ref_key);
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
    const vector<int> &global_path_states,
    const vector<PathMetadata> &noncoding_prefix,
    const ModelInputs &inputs, const Options &options
) {
    const int L = inputs.seq.size();
    ostringstream gff_out, reference_out, tsv_out;
    tsv_out << setprecision(17);
    const auto &ref = references[ref_idx];

    int requested_left = (ref_idx == 0) ? 0 : references[ref_idx - 1].path_end;
    int requested_right = (ref_idx == references.size() - 1)
        ? L - 1 : references[ref_idx + 1].path_start - 1;
    if (options.flank >= 0) {
        requested_left = max(requested_left, ref.genomic_start_0based - options.flank);
        requested_right = static_cast<int>(min<long long>(requested_right,
            static_cast<long long>(ref.genomic_end_0based_exclusive) - 1 + options.flank));
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

    double ref_fixed_score = inputs.emit[L_bound][0];
    PathMetadata initial_metadata = options.flank >= 0
        ? noncoding_prefix[L_bound] : PathMetadata{};
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

    const auto paths = search_paths(L_bound, R_bound, initial_metadata, inputs, options);

    const auto unique_candidates = rank_annotations(paths, ref, seqid, ref_fixed_score, options.k);

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
                << R_bound + 1 << '\t' << initial_metadata.inter_len << '\n';
    }

    return {gff_out.str(), reference_out.str(), tsv_out.str(), {}};
}
} // namespace

bool run_local_k_best(
    const string &seqid,
    const vector<TranscriptAnnotation> &references,
    const vector<int> &global_path_states,
    const vector<PathMetadata> &noncoding_prefix,
    const ModelInputs &inputs,
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
            << "identical_to_reference\tclassification\twindow_start_1based\twindow_end_1based\tinitial_inter_len\n";

    // Each locus writes to independent buffers. Flush by reference index after
    // the parallel region, so scheduling never changes IDs or record order.
    vector<LocusOutput> outputs(references.size());
    const int requested_threads = options.threads > 0 ? options.threads : omp_get_num_procs();
    const int thread_count = max(1, min(requested_threads, static_cast<int>(references.size())));
    // There cannot be more loci than sequence bases (an int-sized input).
    // A signed canonical loop is also accepted by older OpenMP front ends.
    const int locus_count = static_cast<int>(references.size());
    #pragma omp parallel for schedule(dynamic) num_threads(thread_count)
    for (int ref_idx = 0; ref_idx < locus_count; ++ref_idx) {
        try {
            outputs[ref_idx] = decode_locus(ref_idx, seqid, references,
                global_path_states, noncoding_prefix, inputs, options);
        }
        catch (const exception &error) { outputs[ref_idx].error = error.what(); }
        catch (...) { outputs[ref_idx].error = "Unknown decoder error"; }
    }
    // Exceptions cannot cross an OpenMP boundary. Report them on the caller.
    for (const auto &result : outputs) {
        if (!result.error.empty()) throw runtime_error(result.error);
        gff_out << result.gff;
        reference_out << result.reference;
        tsv_out << result.report;
    }
    write_output_bundle({{options.output_filename, gff_out.str()},
                         {options.output_filename + ".reference.gff", reference_out.str()},
                         {options.report_filename, tsv_out.str()}});
    return true;
}
