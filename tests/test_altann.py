"""Integration checks for the native decoder and final GFF coordinate mapping.

Run after building: python3 -m unittest discover -s tests -v
Set ALTANN_REFERENCE_BINARY to a compiled Dmel snapshot to additionally check
numerical compatibility with the original seven-state implementation.
"""

import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from fixtures import ordered_inputs, write_fixture


ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get("ALTANN_CORE", ROOT / "build" / "altann-core")).resolve()


def read_report(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_gff(path):
    """Return feature columns and parsed attributes without altering coordinates."""
    records = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        attributes = dict(item.split("=", 1) for item in columns[8].split(";") if "=" in item)
        records.append((columns, attributes))
    return records


@unittest.skipUnless(CORE.is_file(), "Build altann-core before running integration tests")
class DecoderTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="altann-test-")
        self.addCleanup(self.directory.cleanup)
        self.work = Path(self.directory.name)
        self.inputs = write_fixture(self.work / "inputs")

    def run_core(self, name="core", k=8, threads=1, check=True):
        output = self.work / f"{name}.gff"
        report = self.work / f"{name}.tsv"
        result = subprocess.run(
            [str(CORE), *ordered_inputs(self.inputs), "--k", str(k),
             "--flank", "50", "--threads", str(threads),
             "--output", str(output), "--report", str(report)],
            capture_output=True, text=True, check=check, cwd=self.work,
        )
        return result, output, report

    def test_multiple_loci_and_alternative_splices(self):
        _, output, report = self.run_core()
        rows = read_report(report)
        self.assertEqual({row["reference_gene_index"] for row in rows}, {"1", "2"})
        for locus in ("1", "2"):
            candidates = [row for row in rows if row["reference_gene_index"] == locus]
            self.assertGreaterEqual(len({row["intron_chain"] for row in candidates}), 2)
            self.assertTrue(any(row["identical_to_reference"] == "1" for row in candidates))
            self.assertTrue(any(row["intron_chain"] == "" for row in candidates),
                            "The uninterrupted coding path should retain the intron")
            scores = [float(row["fixed_interval_score"]) for row in candidates]
            self.assertEqual(scores, sorted(scores, reverse=True))
            for row in candidates:
                self.assertAlmostEqual(
                    float(row["fixed_interval_score"]) - float(row["reference_fixed_interval_score"]),
                    float(row["score_delta"]), places=10,
                )
                self.assertLessEqual(int(row["window_start_1based"]), int(row["candidate_start_1based"]))
                self.assertLessEqual(int(row["candidate_end_1based"]), int(row["window_end_1based"]))
        features = read_gff(output)
        self.assertTrue(any(columns[2] == "CDS" for columns, _ in features))
        for columns, _ in features:
            self.assertGreaterEqual(int(columns[3]), 1)
            self.assertLessEqual(int(columns[4]), 1200)
            self.assertLessEqual(int(columns[3]), int(columns[4]))
        reference = read_gff(Path(str(output) + ".reference.gff"))
        coding_intervals = [(int(c[3]), int(c[4]), c[7]) for c, _ in reference if c[2] == "CDS"]
        # Start/stop codons are included, and the GT...AG intron is excluded.
        self.assertEqual(coding_intervals, [(61, 180, "0"), (241, 453, "0"),
                                            (661, 780, "0"), (841, 1053, "0")])

    def test_parallel_output_is_byte_identical(self):
        _, serial, serial_report = self.run_core("serial", threads=1)
        _, parallel, parallel_report = self.run_core("parallel", threads=2)
        self.assertEqual(serial.read_bytes(), parallel.read_bytes())
        self.assertEqual(serial_report.read_bytes(), parallel_report.read_bytes())
        self.assertEqual(Path(str(serial) + ".reference.gff").read_bytes(),
                         Path(str(parallel) + ".reference.gff").read_bytes())

    def test_k_limits_path_ranks(self):
        _, _, report = self.run_core(k=1)
        rows = read_report(report)
        self.assertTrue(rows)
        self.assertEqual({int(row["rank"]) for row in rows}, {1})

    def test_out_of_range_score_index_is_rejected(self):
        with self.inputs["ag"].open("a") as stream:
            stream.write("1200\t2\n")
        result, _, _ = self.run_core(check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(result.stderr.strip())

    def test_empty_fasta_is_rejected(self):
        self.inputs["fasta"].write_text(">empty\n")
        result, _, _ = self.run_core(check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(result.stderr.strip())

    def test_missing_emission_position_is_rejected(self):
        lines = self.inputs["emissions"].read_text().splitlines()
        self.inputs["emissions"].write_text("\n".join(lines[:500] + lines[501:]) + "\n")
        result, _, _ = self.run_core(check=False)
        self.assertNotEqual(result.returncode, 0)

    def test_missing_terminal_emission_positions_are_rejected(self):
        rows = self.inputs["emissions"].read_text().splitlines()
        self.inputs["emissions"].write_text("\n".join(rows[:-2]) + "\n")
        result, _, _ = self.run_core(check=False)
        self.assertNotEqual(result.returncode, 0)

    def test_optional_nucleotide_column_matches_numeric_only_input(self):
        _, baseline, baseline_report = self.run_core("numeric")
        sequence = self.inputs["fasta"].read_text().splitlines()[1]
        rows = self.inputs["emissions"].read_text().splitlines()
        self.inputs["emissions"].write_text("".join(
            f"{row}\t{sequence[index]}\n" for index, row in enumerate(rows)))
        _, annotated, annotated_report = self.run_core("nucleotides")
        self.assertEqual(baseline.read_bytes(), annotated.read_bytes())
        self.assertEqual(baseline_report.read_bytes(), annotated_report.read_bytes())

    def test_input_paths_with_shell_metacharacters(self):
        expected = self.run_cli("plus")
        # These characters are literal filename data, never shell syntax.
        self.inputs = write_fixture(self.work / "inputs % ' $(literal)")
        self.assertEqual(expected, self.run_cli("plus"))

    def test_output_cannot_overwrite_input(self):
        original = self.inputs["fasta"].read_bytes()
        result = subprocess.run(
            [str(CORE), *ordered_inputs(self.inputs), "--k", "2", "--threads", "1",
             "--output", str(self.inputs["fasta"]), "--report", str(self.work / "report.tsv")],
            capture_output=True, text=True, cwd=self.work,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.inputs["fasta"].read_bytes(), original)

    def test_all_intergenic_input_has_empty_annotation(self):
        self.inputs["emissions"].write_text("".join(
            f"{position}\t2\t-8\t-8\t-8\t-8\n" for position in range(1200)))
        _, output, report = self.run_core()
        self.assertEqual(read_gff(output), [])
        self.assertEqual(read_report(report), [])
        self.assertEqual(self.run_cli("plus"), [])

    @unittest.skipUnless(os.environ.get("ALTANN_REFERENCE_BINARY"), "Optional original snapshot binary not supplied")
    def test_original_snapshot_numerical_compatibility(self):
        _, output, report = self.run_core()
        baseline = self.work / "baseline.gff"
        baseline_report = self.work / "baseline.tsv"
        subprocess.run(
            [os.environ["ALTANN_REFERENCE_BINARY"], *ordered_inputs(self.inputs),
             "--local-k-best", "8", "--local-k-flank", "50",
             "--local-k-output", str(baseline), "--local-k-report", str(baseline_report),
             "--no-dp-dump"], cwd=self.work, capture_output=True, text=True, check=True,
        )
        self.assertEqual(read_report(report), read_report(baseline_report))
        self.assertEqual(output.read_text(), baseline.read_text())

    def run_cli(self, strand, extra=(), check=True):
        output = self.work / f"{strand}.gff3"
        args = [sys.executable, "-m", "altann", "decode"]
        for key, path in self.inputs.items():
            args.extend(["--" + key, str(path)])
        args.extend(["--strand", strand, "--output", str(output), "--core", str(CORE),
                     "--k", "8", "--flank", "50", "--threads", "2"])
        args.extend(extra)
        result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, check=check)
        if not check:
            return result
        return read_gff(output)

    def test_final_gff_reference_scores_and_minus_mapping(self):
        # Move both splice boundaries one base so the second CDS has nonzero
        # phase. This catches code that resets all reverse-strand phases to 0.
        self.inputs = write_fixture(self.work / "inputs", splice_shift=1)
        plus = self.run_cli("plus")
        minus = self.run_cli("minus")
        transcripts = [(columns, attrs) for columns, attrs in plus
                       if "kbest_rank" in attrs]
        self.assertTrue(transcripts)
        references = [attrs for _, attrs in transcripts if attrs["kbest_rank"] == "0"]
        self.assertEqual(len(references), 2)
        for columns, attrs in transcripts:
            self.assertEqual(columns[5], ".")
            self.assertAlmostEqual(float(attrs["kbest_score"]) - float(attrs["kbest_reference_score"]),
                                   float(attrs["kbest_delta"]), places=10)
        # Reflection preserves the phase of each corresponding CDS. Sorting by
        # biological feature data avoids depending on generated identifiers.
        reflected = sorted((c[2], 1201 - int(c[4]), 1201 - int(c[3]), c[7]) for c, _ in plus)
        actual = sorted((c[2], int(c[3]), int(c[4]), c[7]) for c, _ in minus)
        self.assertEqual(reflected, actual)
        self.assertTrue(all(c[6] == "-" for c, _ in minus))
        self.assertTrue(any(c[2] == "CDS" and c[7] in ("1", "2") for c, _ in minus))

    def test_best_gff_mismatch_is_rejected(self):
        best = self.work / "wrong-best.gff3"
        best.write_text("##gff-version 3\n"
                        "fixture\tUniAnn\ttranscript\t1\t120\t.\t+\t.\tID=wrong\n"
                        "fixture\tUniAnn\texon\t1\t120\t.\t+\t.\tParent=wrong\n")
        result = self.run_cli("plus", ["--best-gff", str(best)], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(result.stderr.strip())

    def test_empty_supplied_best_gff_is_valid(self):
        expected = self.run_cli("plus")
        best = self.work / "filtered-best.gff3"
        best.write_text("##gff-version 3\n")
        self.assertEqual(expected, self.run_cli("plus", ["--best-gff", str(best)]))

    def test_explicit_segment_offset_and_ownership(self):
        for strand in ("plus", "minus"):
            with self.subTest(strand=strand):
                baseline = self.run_cli(strand)
                shifted = self.run_cli(strand, ["--offset", "1000", "--sequence-length", "3000",
                                                "--seqid", "chr%20name", "--segment-overlap", "0",
                                                "--segment-margin", "0"])
                expected = sorted((c[2], int(c[3]) + 1000, int(c[4]) + 1000, c[6], c[7])
                                  for c, _ in baseline)
                actual = sorted((c[2], int(c[3]), int(c[4]), c[6], c[7]) for c, _ in shifted)
                self.assertEqual(expected, actual)
                self.assertTrue(all(c[0] == "chr%2520name" for c, _ in shifted))
        # A margin larger than the segment removes every model. This checks
        # that explicit coordinates activate ownership even without a split header.
        self.assertEqual(self.run_cli("plus", ["--offset", "1000", "--sequence-length", "3000",
                                               "--segment-overlap", "0", "--segment-margin", "1200"]), [])

    def test_manifest_both_strands_is_sorted_independently_of_job_order(self):
        manifest = self.work / "jobs.tsv"
        columns = ["strand", *self.inputs]
        contents = []
        for iteration, strands in enumerate((("minus", "plus"), ("plus", "minus"))):
            with manifest.open("w") as stream:
                writer = csv.DictWriter(stream, columns, delimiter="\t")
                writer.writeheader()
                for strand in strands:
                    writer.writerow(dict(strand=strand, **{
                        key: str(path.relative_to(self.work)) for key, path in self.inputs.items()}))
            output = self.work / f"manifest-{iteration}.gff3"
            subprocess.run([sys.executable, "-m", "altann", "decode", "--manifest", str(manifest),
                            "--core", str(CORE), "--output", str(output), "--k", "8", "--flank", "50",
                            "--threads", "2"], cwd=ROOT, check=True, capture_output=True, text=True)
            contents.append(output.read_bytes())
            features = read_gff(output)
            self.assertEqual({c[6] for c, _ in features}, {"+", "-"})
            keys = [(c[0], int(c[3]), int(c[4]), c[6], attrs["ID"])
                    for c, attrs in features if c[2] in ("transcript", "mRNA")]
            self.assertEqual(keys, sorted(keys))
        self.assertEqual(contents[0], contents[1])

    def test_psauron_large_csv_fields(self):
        from altann.cli import validate_psauron

        scores = self.work / "psauron.csv"
        with scores.open("w") as stream:
            stream.write("preamble\n" * 4)
            row = ["fixture"] + [""] * 14
            # Each semicolon-separated frame exceeds csv's usual 128 KiB limit.
            row[9:12] = ["0.5;" * 40_000] * 3
            csv.writer(stream).writerow(row)
        validate_psauron(scores, self.inputs["fasta"])


if __name__ == "__main__":
    unittest.main()
