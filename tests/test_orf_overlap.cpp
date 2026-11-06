// Exercise output selection independently of DP search/scoring. Including the
// implementation keeps its path/ranking types private to the production unit.
#include "../src/local.cpp"

namespace {
void require(bool condition, const string &message) {
    if (!condition) throw runtime_error(message);
}

TranscriptAnnotation annotation(vector<GenomicInterval> cds,
                                string seqid = "fixture", char strand = '+') {
    TranscriptAnnotation result;
    result.sequence_id = seqid;
    result.strand = strand;
    result.cds_intervals = cds;
    result.exons = cds;
    if (!cds.empty()) {
        result.genomic_start_0based = cds.front().start_0based;
        result.genomic_end_0based_exclusive = cds.back().end_0based_exclusive;
    }
    return result;
}

// The state switches use the same codon-coordinate convention as the native
// annotation builder: start+2 is N->E; end-1 is E->N.
StatePath run_length(const vector<int> &states) {
    vector<int> starts;
    vector<int8_t> runs;
    for (size_t p = 0; p < states.size(); ++p)
        if (p == 0 || states[p] != states[p - 1]) {
            starts.push_back(p);
            runs.push_back(states[p]);
        }
    return StatePath(starts, runs, states.size());
}

LocalPathResult path(vector<GenomicInterval> genes, double score,
                     vector<double> start_scores) {
    LocalPathResult result;
    result.score = score;
    result.start_scores = start_scores;
    const int length = genes.empty() ? 600 : genes.back().end_0based_exclusive + 2;
    vector<int> states(length, 0);
    for (const auto &gene : genes)
        for (int p = gene.start_0based + 2; p < gene.end_0based_exclusive - 1; ++p)
            states[p] = 1;
    result.path = run_length(states);
    return result;
}

void test_backpointer_spill_round_trip() {
    const size_t cells = 7 * 3;
    const size_t row_bytes = cells * sizeof(CompactBP);
    // Limits and block sizes that keep everything in memory, split rows across
    // blocks with a partial last block, and hold one row per block.
    for (const auto &[limit, block] : vector<pair<size_t, size_t>>{
             {size_t(1) << 30, row_bytes}, {0, 4 * row_bytes}, {0, row_bytes}, {0, 1}}) {
        for (size_t rows : {size_t(1), size_t(4), size_t(5), size_t(41)}) {
            TracebackBudget budget(limit);
            BackpointerStore store(rows, cells, budget, ".", block);
            require(store.spilling() == (limit == 0), "the budget decides where a traceback lives");
            for (size_t row = 0; row < rows; ++row) {
                CompactBP *cell = store.write_row(row);
                require(cell[cells - 1].bt_state == -1, "a new row starts unset");
                cell[0] = CompactBP{static_cast<int16_t>(row % 7), static_cast<int16_t>(row)};
                cell[cells - 1] = CompactBP{static_cast<int16_t>(row % 5), static_cast<int16_t>(rows - row)};
            }
            store.finish_writes();
            for (size_t row = rows; row-- > 0;) {
                const CompactBP *cell = store.read_row(row);
                require(cell[0].bt_state == static_cast<int16_t>(row % 7) &&
                        cell[0].bt_k == static_cast<int16_t>(row), "first cell survives the round trip");
                require(cell[cells - 1].bt_state == static_cast<int16_t>(row % 5) &&
                        cell[cells - 1].bt_k == static_cast<int16_t>(rows - row),
                        "last cell survives the round trip");
            }
        }
    }
    bool rejected = false;
    TracebackBudget empty(0);
    try { BackpointerStore store(4, cells, empty, "./missing-spill-directory"); }
    catch (const runtime_error &) { rejected = true; }
    require(rejected, "an unusable spill directory is reported");
}

void test_traceback_budget_is_shared_across_loci() {
    const size_t cells = 7;
    const size_t row_bytes = cells * sizeof(CompactBP);
    TracebackBudget budget(10 * row_bytes);
    {
        BackpointerStore first(6, cells, budget, ".");
        BackpointerStore second(6, cells, budget, ".");
        BackpointerStore third(4, cells, budget, ".");
        require(!first.spilling(), "the first locus fits");
        require(second.spilling(), "a concurrent locus over the remaining budget spills");
        require(!third.spilling(), "a smaller locus still fits the remainder");
        require(budget.available() == 0, "reservations stay within the budget");
    }
    require(budget.available() == 10 * row_bytes, "finished loci release their memory");

    // Concurrent reservations can never exceed the budget.
    TracebackBudget shared(100);
    atomic<int> held{0}, peak{0};
    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < 20000; ++i) {
        if (!shared.try_reserve(10)) continue;
        const int now = ++held;
        int seen = peak.load();
        while (now > seen && !peak.compare_exchange_weak(seen, now)) {}
        --held;
        shared.release(10);
    }
    require(peak.load() <= 10, "concurrent reservations stay within the budget");
    require(shared.available() == 100, "concurrent releases restore the budget");
}

void test_cds_overlap() {
    const auto ref = annotation({{100, 200}, {400, 500}});
    require(!shares_cds_base(ref, annotation({{50, 100}})), "left edge touch is disjoint");
    require(!shares_cds_base(ref, annotation({{500, 550}})), "right edge touch is disjoint");
    require(shares_cds_base(ref, annotation({{99, 101}})), "one shared CDS base is sufficient");
    require(shares_cds_base(ref, annotation({{450, 470}})), "internal CDS overlap is retained");
    require(shares_cds_base(ref, annotation({{50, 550}})), "a spanning CDS is retained");
    require(!shares_cds_base(ref, annotation({{200, 400}})), "intron-only overlap is disjoint");
    require(!shares_cds_base(ref, annotation({{0, 100}, {200, 400}, {500, 600}})),
            "overlapping spans with disjoint CDS are rejected");
    require(!shares_cds_base(ref, annotation({{100, 200}}, "other")), "seqid must match");
    require(!shares_cds_base(ref, annotation({{100, 200}}, "fixture", '-')), "strand must match");
    require(shares_cds_base(annotation({{100, 200}}, "fixture", '-'),
                            annotation({{199, 250}}, "fixture", '-')),
            "same-strand reverse intervals can overlap");
    require(!shares_cds_base(ref, annotation({})), "no candidate CDS is rejected");
    require(!shares_cds_base(annotation({}), ref), "no reference CDS is rejected");
}

void test_filter_before_ranking_and_deduplication() {
    const auto ref = annotation({{100, 200}});
    const vector<LocalPathResult> paths = {
        path({}, 110, {}),
        path({{10, 30}}, 100, {9}),
        path({{10, 30}, {100, 200}, {300, 350}}, 90, {11, 22, 33}),
        path({{100, 200}}, 80, {44}),
        path({{150, 220}}, 70, {55}),
        path({{160, 230}}, 60, {66}),
    };
    const auto selected = rank_annotations(paths, ref, "fixture", 85, 2);
    require(selected.size() == 2, "only distinct retained annotations consume ranks");
    const auto &first = selected[0];
    require(first.rank == 1, "empty/disjoint paths must not consume ranks");
    require(first.transcript.cds_intervals == ref.cds_intervals, "only overlapping gene is exported");
    require(first.genes_in_path == 3 && first.gene_in_path == 2,
            "original whole-path gene membership is preserved");
    require(first.actual_start_score == 22, "start score uses original gene index");
    require(first.score == 90 && first.transcript.score == 90 && first.score_delta == 5,
            "filtering preserves exact whole-path score and reference delta");
    require(first.identical_to_reference, "reference equality describes retained annotation");
    require(selected[1].rank == 2 && selected[1].score == 70,
            "deduplication ignores discarded genes and retains the highest scoring path");
    require(selected[1].transcript.genomic_start_0based == 150,
            "partial overlap is retained after duplicate removal");
    require(!selected[1].identical_to_reference, "changed CDS is not reference-identical");
}

void test_raw_keeps_every_path_and_gene() {
    const auto ref = annotation({{100, 200}});
    const vector<LocalPathResult> paths = {
        path({}, 110, {}),
        path({{10, 30}}, 100, {9}),
        path({{10, 30}, {100, 200}}, 90, {11, 22}),
        path({{100, 200}}, 80, {44}),
        path({{100, 200}}, 70, {55}),
    };
    const auto all = rank_annotations(paths, ref, "fixture", 85, 10, true);
    require(all.size() == 5, "every gene of every path is returned");
    const vector<int> ranks = {2, 3, 3, 4, 5};
    for (size_t i = 0; i < all.size(); ++i)
        require(all[i].rank == ranks[i], "rank is the decoder path index");
    require(all[0].classification == "no_cds_overlap" && all[0].splice_jaccard == 0,
            "a disjoint intronless gene is not reported as identical");
    require(all[1].classification == "no_cds_overlap" && !all[1].identical_to_reference,
            "the disjoint gene of a shared path is marked");
    require(all[2].gene_in_path == 2 && all[2].actual_start_score == 22 &&
            all[2].classification == "identical" && all[2].identical_to_reference,
            "a disjoint neighbor does not change the overlapping gene's report");
    require(all[3].identical_to_reference && all[4].identical_to_reference,
            "repeated annotations are kept");
}

void test_intron_only_and_multiple_overlapping_genes() {
    const auto ref = annotation({{100, 200}, {400, 500}});
    const vector<LocalPathResult> paths = {
        path({{250, 350}}, 100, {9}),
        path({{10, 30}, {150, 180}, {250, 350}, {420, 460}}, 90, {11, 22, 33, 44}),
    };
    const auto selected = rank_annotations(paths, ref, "fixture", 85, 10);
    require(selected.size() == 2, "both overlapping genes survive, intronic/flanking genes do not");
    require(selected[0].rank == 1 && selected[1].rank == 1,
            "multiple retained genes share the original path rank");
    require(selected[0].gene_in_path == 2 && selected[1].gene_in_path == 4,
            "original indices survive removal of intervening intronic genes");
    require(selected[0].genes_in_path == 4 && selected[1].genes_in_path == 4,
            "path gene count includes discarded genes");
    require(selected[0].actual_start_score == 22 && selected[1].actual_start_score == 44,
            "all retained starts match their own original scores");
    require(selected[0].score == 90 && selected[1].score == 90,
            "retained genes continue to report the same complete path score");
    require(!selected[0].identical_to_reference && !selected[1].identical_to_reference,
            "split genes do not become reference-identical");
}
TranscriptAnnotation spliced(vector<GenomicInterval> cds) {
    auto result = annotation(cds);
    for (size_t i = 0; i + 1 < cds.size(); ++i)
        result.junctions.push_back({cds[i].end_0based_exclusive,
                                    cds[i + 1].start_0based - 2, '+'});
    return result;
}

void test_splice_site_jaccard() {
    const auto ref = spliced({{100, 200}, {300, 400}, {500, 600}});
    auto near = [](double a, double b) { return fabs(a - b) < 1e-12; };
    require(near(splice_site_jaccard(ref, ref), 1.0), "identical chains score 1");
    require(near(splice_site_jaccard(ref, spliced({{90, 200}, {300, 400}, {500, 650}})), 1.0),
            "start/stop changes do not affect splice sites");
    require(near(splice_site_jaccard(ref, spliced({{100, 210}, {300, 400}, {500, 600}})), 3.0 / 5.0),
            "one shifted donor breaks only that donor");
    require(near(splice_site_jaccard(ref, spliced({{100, 200}, {500, 600}})), 2.0 / 4.0),
            "a skipped exon removes its acceptor and donor");
    require(near(splice_site_jaccard(ref, annotation({{100, 600}})), 0.0),
            "an intronless candidate shares no splice sites with a spliced reference");
    require(near(splice_site_jaccard(annotation({{100, 200}}), annotation({{150, 300}})), 1.0),
            "two intronless models are identical");
    auto other = ref;
    other.strand = '-';
    require(near(splice_site_jaccard(ref, other), 0.0), "strand must match");
}

// Exons are written as E states and the gaps between them as I states; the
// annotation builder then derives motif coordinates with its usual offsets.
LocalPathResult spliced_path(vector<GenomicInterval> exons, double score) {
    auto result = path({{exons.front().start_0based, exons.back().end_0based_exclusive}},
                       score, {1});
    vector<int> states(result.path.size());
    for (size_t p = 0; p < states.size(); ++p) states[p] = result.path[p];
    for (size_t i = 0; i + 1 < exons.size(); ++i)
        for (int p = exons[i].end_0based_exclusive; p < exons[i + 1].start_0based; ++p)
            states[p] = 4;
    result.path = run_length(states);
    return result;
}

void test_ranked_candidates_report_jaccard() {
    const auto ref_path = spliced_path({{100, 200}, {300, 400}, {500, 600}}, 85);
    const auto ref = build_transcript_annotations(ref_path.path, "fixture").at(0);
    require(ref.junctions.size() == 2, "fixture reference has two introns");
    const vector<LocalPathResult> paths = {
        spliced_path({{100, 200}, {320, 420}, {520, 600}}, 100),  // 1 of 7 sites
        spliced_path({{100, 210}, {300, 400}, {500, 600}}, 90),   // 3 of 5 sites
    };
    const auto all = rank_annotations(paths, ref, "fixture", 85, 10);
    require(all.size() == 2, "every overlapping candidate is kept");
    for (const auto &cand : all)
        require(cand.splice_jaccard == splice_site_jaccard(ref, cand.transcript),
                "candidates report their splice-site Jaccard");
    require(fabs(all[0].splice_jaccard - 1.0 / 7.0) < 1e-12 &&
            fabs(all[1].splice_jaccard - 3.0 / 5.0) < 1e-12, "fixture Jaccard values");
}
} // namespace

int main() {
    try {
        test_cds_overlap();
        test_backpointer_spill_round_trip();
        test_traceback_budget_is_shared_across_loci();
        test_filter_before_ranking_and_deduplication();
        test_raw_keeps_every_path_and_gene();
        test_intron_only_and_multiple_overlapping_genes();
        test_splice_site_jaccard();
        test_ranked_candidates_report_jaccard();
    } catch (const exception &error) {
        cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
