"""Spliced biology and exact serialized-score validation on original coordinates."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_prediction_cds.py"
SPEC = importlib.util.spec_from_file_location("prediction_cds_validation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def model(strand, first_length, origin="uniann_best", scores=None):
    coding = "ATGAAATAA"
    oriented = coding[:first_length] + "CCCC" + coding[first_length:]
    length = len(oriented)
    sequence = oriented if strand == "+" else oriented.translate(MODULE.COMPLEMENT)[::-1]
    intervals = [(1, first_length, 0), (first_length + 5, length, (-first_length) % 3)]
    if strand == "-":
        intervals = [(length - end + 1, length - start + 1, phase) for start, end, phase in intervals]
    attributes = "ID=t;kbest_origin=" + origin
    if scores is not None:
        attributes += ";" + ";".join(key + "=" + value for key, value in zip(
            ("kbest_score", "kbest_reference_score", "kbest_delta"), scores))
    rows = [f"chr\tUniAnn\ttranscript\t1\t{length}\t.\t{strand}\t.\t{attributes}"]
    for index, (start, end, phase) in enumerate(intervals, 1):
        rows.append(f"chr\tUniAnn\texon\t{start}\t{end}\t.\t{strand}\t.\tID=e{index};Parent=t")
        rows.append(f"chr\tUniAnn\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tID=c{index};Parent=t")
    return sequence, "##gff-version 3\n" + "\n".join(rows) + "\n"


class PredictionCDSValidationTests(unittest.TestCase):
    def validate(self, sequence, gff):
        with tempfile.TemporaryDirectory() as directory:
            genome, annotation = Path(directory) / "genome.fa", Path(directory) / "models.gff3"
            genome.write_text(">chr\n" + sequence + "\n")
            annotation.write_text(gff)
            return MODULE.validate(genome, annotation)

    def test_split_start_internal_and_stop_codons_on_both_strands(self):
        # Length 1 splits ATG; 4 splits the second codon; 8 splits terminal TAA.
        # Removing each segment's phase bases would corrupt all these CDSs.
        for strand in ("+", "-"):
            for split in (1, 4, 8):
                with self.subTest(strand=strand, split=split):
                    result = self.validate(*model(strand, split, scores=("0.2", "0.2", "0")))
                    self.assertTrue(result["all_checks_pass"], result)
                    self.assertEqual(result["by_origin_and_strand"]["baseline:" + strand]["CDS_bases_checked"], 9)

    def test_phase_discontinuity_is_diagnosed_in_transcript_direction(self):
        for strand in ("+", "-"):
            sequence, gff = model(strand, 4)
            gff = gff.replace("\t2\tID=c2", "\t0\tID=c2")
            result = self.validate(sequence, gff)
            self.assertIn("CDS_phase_discontinuity", result["failed_models"][0]["errors"])

    def test_real_spombe_native_start_donor_overlap_is_diagnosed(self):
        # Exact real DNA; both native UniAnn and AltAnn reproduce these exon
        # bounds when diagnostic emissions force the observed candidate path.
        fixture = json.loads((ROOT / "tests" / "data" / "spom_prediction_start_overlap.json").read_text())
        sequence = fixture["oriented_sequence"]
        rows = [f"chr\tUniAnn\ttranscript\t31\t507\t.\t+\t.\tID=t"]
        pieces = []
        for index, (start, end, phase) in enumerate(fixture["local_CDS_1based"], 1):
            pieces.append(sequence[start - 1:end])
            rows.append(f"chr\tUniAnn\tCDS\t{start}\t{end}\t.\t+\t{phase}\tID=c{index};Parent=t")
        self.assertEqual(sequence[30:34], "ATGT")
        self.assertEqual("".join(pieces)[:3], "ATC")
        result = self.validate(sequence, "\n".join(rows) + "\n")
        self.assertFalse(result["structural_diagnostics"])
        self.assertEqual(result["failed_models"][0]["errors"], ["non_ATG_start"])

    def test_internal_stop_spanning_exon_boundary_is_detected(self):
        sequence, gff = model("+", 4)
        # Same interval lengths, spliced CDS ATG TAA TAA: first TAA spans intron.
        sequence = "ATGT" + "CCCC" + "AATAA"
        result = self.validate(sequence, gff)
        self.assertIn("internal_stop", result["failed_models"][0]["errors"])
        self.assertEqual(result["failed_models"][0]["internal_stops"][0]["kind"], "split_across_junction")
        self.assertEqual(result["failed_models"][0]["internal_stops"][0]["genomic_positions_1based"], [4, 9, 10])

    def test_real_spombe_prediction_junction_stops_on_both_strands(self):
        # Actual model outputs and exact original-genome DNA slices, preserving
        # observed phases. Records are shifted only to make the fixture compact.
        fixture = json.loads((ROOT / "tests" / "data" / "spom_prediction_internal_stop.json").read_text())
        for case in fixture["cases"]:
            strand = case["strand"]
            rows = [f"chr\tUniAnn\ttranscript\t1\t{len(case['original_genome_sequence'])}\t.\t{strand}\t.\t"
                    + "ID=t;kbest_origin=" + ("uniann_best" if case["origin"] == "baseline" else "local_kbest")]
            for index, (start, end, phase) in enumerate(case["local_CDS_1based"], 1):
                rows.append(f"chr\tUniAnn\texon\t{start}\t{end}\t.\t{strand}\t.\tID=e{index};Parent=t")
                rows.append(f"chr\tUniAnn\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tID=c{index};Parent=t")
            with self.subTest(transcript=case["transcript"], strand=strand):
                result = self.validate(case["original_genome_sequence"], "\n".join(rows) + "\n")
                self.assertFalse(result["structural_diagnostics"])
                self.assertEqual(result["failed_models"][0]["errors"], ["internal_stop"])
                stop = result["failed_models"][0]["internal_stops"][0]
                source_stop = case["internal_stops"][0]
                self.assertEqual(stop["codon"], source_stop["codon"])
                self.assertEqual(stop["spliced_CDS_0based_offset"], source_stop["spliced_CDS_0based_offset"])
                self.assertEqual([position + case["genome_slice_1based"][0] - 1 for position in stop["genomic_positions_1based"]],
                                 source_stop["genomic_positions_1based"])
                self.assertTrue(all(stop["kind"] == "split_across_junction" for stop in case["internal_stops"]))

    def test_exact_decimal_delta_and_candidate_origin(self):
        sequence, gff = model("-", 4, "local_kbest", scores=(
            "100000000000000000000.3", "100000000000000000000.2", "0.1"))
        result = self.validate(sequence, gff)
        self.assertTrue(result["all_checks_pass"], result)
        self.assertEqual(result["by_origin_and_strand"]["candidates:-"]["score_triples_checked"], 1)
        gff = gff.replace("kbest_delta=0.1", "kbest_delta=0.100000000000000000001")
        result = self.validate(sequence, gff)
        self.assertIn("serialized_score_delta_mismatch", result["failed_models"][0]["errors"])

    def test_parent_and_original_coordinate_errors(self):
        sequence, gff = model("+", 4)
        gff = gff.replace("ID=c1;Parent=t", "ID=c1;Parent=missing")
        gff += "chr\tUniAnn\tCDS\t14\t16\t.\t-\t0\tID=bad;Parent=t\n"
        result = self.validate(sequence, gff)
        errors = {entry["error"] for entry in result["structural_diagnostics"]}
        self.assertTrue({"missing_parent_feature", "coordinates_outside_original_genome",
                         "parent_chromosome_or_strand_mismatch", "child_outside_parent"} <= errors)

    def test_gff_parent_cycle_and_invalid_score_are_diagnosed(self):
        sequence, gff = model("+", 4)
        gff = gff.replace("ID=t;", "ID=t;Parent=t;").replace("transcript\t1\t13\t.", "transcript\t1\t13\tNaN")
        result = self.validate(sequence, gff)
        errors = {entry["error"] for entry in result["structural_diagnostics"]}
        self.assertTrue({"cyclic_Parent_relationship", "invalid_gene_parent_type", "invalid_GFF_score"} <= errors)

    def test_cli_writes_diagnostics_and_returns_failure(self):
        sequence, gff = model("+", 4)
        with tempfile.TemporaryDirectory() as directory:
            genome, annotation, output = [Path(directory) / name for name in ("genome.fa", "models.gff3", "audit.json")]
            genome.write_text(">chr\n" + sequence + "\n")
            annotation.write_text(gff.replace("\t2\tID=c2", "\t0\tID=c2"))
            result = subprocess.run([sys.executable, str(SCRIPT), "--genome", str(genome),
                                     "--gff", str(annotation), "--output", str(output)], capture_output=True)
            self.assertEqual(result.returncode, 1)
            self.assertEqual(json.loads(output.read_text())["counts"]["failed_transcripts"], 1)
            self.assertEqual(genome.read_text(), ">chr\n" + sequence + "\n")


if __name__ == "__main__":
    unittest.main()
