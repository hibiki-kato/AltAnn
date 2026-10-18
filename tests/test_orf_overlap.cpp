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
LocalPathResult path(vector<GenomicInterval> genes, double score,
                     vector<double> start_scores) {
    LocalPathResult result;
    result.score = score;
    result.start_scores = start_scores;
    const int length = genes.empty() ? 600 : genes.back().end_0based_exclusive + 2;
    result.path_states.assign(length, 0);
    for (int p = 0; p < length; ++p) result.path_positions.push_back(p);
    for (const auto &gene : genes)
        for (int p = gene.start_0based + 2; p < gene.end_0based_exclusive - 1; ++p)
            result.path_states[p] = 1;
    return result;
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
} // namespace

int main() {
    try {
        test_cds_overlap();
        test_filter_before_ranking_and_deduplication();
        test_intron_only_and_multiple_overlapping_genes();
    } catch (const exception &error) {
        cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
