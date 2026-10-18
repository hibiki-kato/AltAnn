"""Regression checks for the existing strict-CDS and candidate-pool conventions."""

import importlib.util
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/evaluate_spom_predictions.py"
SPEC = importlib.util.spec_from_file_location("spom_evaluation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
TRMAP = shutil.which("trmap")
if TRMAP is None:
    local = Path("/home/hibiki/miniforge3/envs/eviann/bin/trmap")
    TRMAP = str(local) if local.is_file() else None


def model(identifier, cds, strand="+", **attrs):
    return MODULE.Model("chrI", strand, identifier, attrs, list(cds), list(cds))


class EvaluationTests(unittest.TestCase):
    def test_trmap_query_and_reference_ids_are_distinguished(self):
        pairs = MODULE.parse_trmap(
            ">query_1 chrI:100-300 + 100-150,200-300\n"
            "=\tchrI\t+\t100\t300\treference_1\t100-150,200-300\n"
            "=\tchrI\t+\t100\t300\treference_2\t100-150,200-300\n"
        )
        self.assertEqual(pairs, {("query_1", "reference_1"), ("query_1", "reference_2")})

    def test_candidates_share_loci_but_exact_cds_remains_nonredundant(self):
        reference = {"r": model("r", [(100, 150), (200, 300)], Parent="gene_r")}
        query = {
            "q0": model("q0", [(100, 150), (200, 300)], kbest_locus="locus1", kbest_origin="uniann_best"),
            "q1": model("q1", [(100, 150), (200, 300)], kbest_locus="locus1", kbest_origin="local_kbest"),
            "q2": model("q2", [(110, 150), (200, 300)], kbest_locus="locus1", kbest_origin="local_kbest"),
        }
        result = MODULE.evaluate(reference, query, {("q0", "r"), ("q1", "r")}, candidate=True)
        self.assertEqual(result["query_gene_loci"], 1)
        self.assertEqual(result["gene_precision_percent"], 100)
        self.assertEqual(result["query_unique_cds"], 2)
        self.assertEqual(result["cds_precision_percent"], 50)
        self.assertEqual(result["exact_match_query_transcripts"], 2)
        self.assertEqual(result["alternative_isoform_transcripts"], 2)
        self.assertEqual(result["precision_scope"], "candidate_pool")

    def test_same_cds_on_opposite_strand_is_not_an_exact_match(self):
        result = MODULE.evaluate(
            {"r": model("r", [(100, 300)])},
            {"q": model("q", [(100, 300)], strand="-")}, set(),
        )
        self.assertEqual(result["exact_match_unique_cds"], 0)

    def test_transcript_only_query_can_detect_a_reference_gene(self):
        reference = {"r": model("r", [(100, 300)], Parent="gene_r")}
        query = {"q": MODULE.Model("chrI", "+", "q", {}, [(90, 310)], [])}
        result = MODULE.evaluate(reference, query, {("q", "r")})
        self.assertEqual(result["gene_sensitivity_percent"], 100)
        self.assertEqual(result["query_transcripts"], 1)
        self.assertEqual(result["query_coding_transcripts"], 0)
        self.assertEqual(result["exact_match_unique_cds"], 0)
        self.assertIsNone(result["cds_precision_percent"])

    @unittest.skipUnless(TRMAP, "trmap required for the established evaluation convention")
    def test_real_trmap_strict_boundaries_and_strands(self):
        reference = {"r": model("r", [(100, 150), (200, 300)])}
        query = {
            "exact": model("exact", [(100, 150), (200, 300)]),
            "shifted": model("shifted", [(110, 150), (200, 300)]),
            "minus": model("minus", [(100, 150), (200, 300)], strand="-"),
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            MODULE.write_models(reference, root / "reference.gff3", cds_only=True)
            MODULE.write_models(query, root / "query.gff3", cds_only=True)
            transcript = MODULE.run_trmap(TRMAP, root / "reference.gff3", root / "query.gff3")
            strict = MODULE.run_trmap(TRMAP, root / "reference.gff3", root / "query.gff3", strict=True)
        self.assertEqual(transcript, {("model_0", "model_0"), ("model_1", "model_0")})
        self.assertEqual(strict, {("model_0", "model_0")})


if __name__ == "__main__":
    unittest.main()
