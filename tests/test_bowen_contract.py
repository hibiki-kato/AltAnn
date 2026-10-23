"""Pinned upstream contract diagnostics and literal processed-score regression.

These fixtures expose upstream discrepancies; they do not designate those
discrepancies as desirable behavior or implement compensating transformations.
"""
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

from test_altann import read_gff, read_report


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data/bowen-contract"
CORE = Path(os.environ.get("ALTANN_CORE", ROOT / "build/altann-core")).resolve()


def native_path(data=DATA):
    scores, predecessors = [], []
    for line in (data / "uniann.trace").read_text().splitlines():
        fields = line.split("\t")
        (scores if fields[1] == "dp" else predecessors).append(list(map(int, fields[2:])))
    state = max(range(7), key=lambda index: scores[-1][index])
    path = [0] * len(predecessors)
    for position in range(len(path) - 1, -1, -1):
        path[position] = state
        state = predecessors[position][state]
    return path


class BowenContractDiagnostics(unittest.TestCase):
    def setUp(self):
        self.metadata = json.loads((DATA / "provenance.json").read_text())
        self.raw = [line.split("\t") for line in (DATA / "coding.txt").read_text().splitlines()]
        self.emissions = [list(map(float, line.split()[1:6]))
                          for line in (DATA / "out.ps.txt").read_text().splitlines()]

    def test_real_exporter_reorders_calibrated_multiclass_probabilities(self):
        metadata = self.metadata
        self.assertTrue(metadata["synthetic_logits"])
        self.assertFalse(metadata["checkpoint_inference"])
        self.assertEqual(metadata["internal_class_order"], ["non_CDS", "frame_0", "frame_1", "frame_2"])
        self.assertEqual(len(self.raw), metadata["length"])
        self.assertEqual((DATA / "exported.tsv").read_text().splitlines()[1:],
                         (DATA / "coding.txt").read_text().splitlines())
        for logits, probabilities in zip(metadata["logit_profiles"], metadata["calibrated_profiles"]):
            values = [math.exp(value / metadata["temperature"]) for value in logits]
            expected = [value / sum(values) for value in values]
            for actual, value in zip(probabilities, expected):
                self.assertAlmostEqual(actual, value, places=14)
        names = ["noncoding", "frame0", "frame1", "frame2"]
        for position, fields in enumerate(self.raw):
            label = metadata["phase_labels"][position]
            target = 0 if label == "." else int(label) + 1
            expected = metadata["calibrated_profiles"][target]
            self.assertEqual(len(fields), 8)
            self.assertEqual(int(fields[1]), position + 1)
            self.assertEqual(fields[2], "+")
            self.assertEqual(fields[7], names[target])
            actual = list(map(float, fields[3:7]))
            self.assertAlmostEqual(sum(actual), 1, places=7)
            for value, expected_value in zip(actual, [expected[i] for i in (1, 2, 3, 0)]):
                self.assertAlmostEqual(value, expected_value, places=7)

    def test_pinned_converter_frame_rotation_diagnostic(self):
        genes = self.metadata["genes"]
        self.assertEqual({gene["start0"] % 3 for gene in genes}, {0, 1, 2})
        self.assertEqual({gene["intron_length"] % 3 for gene in genes}, {0, 1, 2})
        self.assertTrue(self.metadata["minus_reflection_verified"])
        for gene in genes:
            # Actual Bowen labels start 0,1,2 at a phase-0 CDS; second CDS has
            # phase2 and starts offset1. The pinned converter's -1 rotates E.
            self.assertEqual(self.metadata["phase_labels"][gene["start0"]:gene["start0"] + 3], "012")
            second = gene["donor0"] + gene["intron_length"]
            self.assertEqual(self.metadata["phase_labels"][second:second + 3], "120")
            for position, expected_frame in ((gene["start0"] + 8, gene["start0"] % 3),
                                              (second + 8, (gene["start0"] + gene["intron_length"]) % 3)):
                local = int(self.metadata["phase_labels"][position])
                observed = max(range(3), key=lambda i: self.emissions[position][i + 1])
                self.assertEqual((position - local) % 3, expected_frame)
                self.assertEqual(observed, (expected_frame - 1) % 3,
                                 "This diagnostic is tied to the recorded converter source")

    def test_pinned_converter_adds_independent_temperature_transform(self):
        # Noncoding base0 has no stop mask and provides a simple diagnostic.
        probabilities = list(map(float, self.raw[0][3:7]))
        scaled = [1 / (1 + ((1 - value) / (value + 1e-10)) ** (1 / 7))
                  for value in probabilities]
        self.assertGreater(abs(sum(scaled) - 1), .1)
        expected_nc = round(math.log(math.e * scaled[3] + 1e-10), 2)
        self.assertEqual(self.emissions[0][0], expected_nc)
        self.assertGreater(abs(expected_nc - math.log(math.e * probabilities[3] + 1e-10)), .05)

    def test_header_is_an_extra_base_in_pinned_converter(self):
        self.assertEqual(self.metadata["unstripped_converter_rows"], self.metadata["length"] + 1)
        diagnostic = (DATA / "header-diagnostic.txt").read_text().splitlines()
        self.assertEqual(len(diagnostic), 3)
        self.assertNotEqual(diagnostic[0], (DATA / "out.ps.txt").read_text().splitlines()[0])
        self.assertEqual(self.metadata["export_header"][1], "position_1based")


@unittest.skipUnless(CORE.is_file(), "Build altann-core before running integration tests")
class BowenProcessedScoreTests(unittest.TestCase):
    data = DATA

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="altann-bowen-")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.metadata = json.loads((self.data / "provenance.json").read_text())

    def decode(self, name, saved=False, reverse=False, **kwargs):
        output = self.work / (name + ".gff3")
        command = [sys.executable, "-m", "altann", "decode", "--core", str(CORE),
                   "--input", str(self.data), "--k", "8", "--flank", "50", "--output", str(output)]
        command += (["--gff", str(self.data / "uniann.gff3"), "--log", str(self.data / "uniann.trace")]
                    if saved else ["--rerun-viterbi"])
        if reverse:
            command += ["--reverse"]
        subprocess.run(command, cwd=ROOT, capture_output=True, text=True, **kwargs)
        return output

    def test_saved_native_trace_matches_rerun_on_both_output_strands(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                saved = self.decode("saved", saved=True, reverse=reverse, check=True)
                rerun = self.decode("rerun", reverse=reverse, check=True)
                for suffix in ("", ".tsv"):
                    self.assertEqual(Path(str(saved) + suffix).read_bytes(), Path(str(rerun) + suffix).read_bytes())
                rows = read_report(str(saved) + ".tsv")
                self.assertEqual({row["reference_gene_index"] for row in rows}, {"1", "2", "3"})
                self.assertTrue(any(row["classification"] == "alternative_acceptor" for row in rows))
                records = read_gff(saved)
                self.assertEqual({columns[6] for columns, _ in records}, {"-" if reverse else "+"})
                sequence = "".join((self.data / "sequence.fna").read_text().splitlines()[1:])
                complement = str.maketrans("ACGT", "TGCA")
                genome = sequence.translate(complement)[::-1] if reverse else sequence
                transcripts = [attrs["ID"] for columns, attrs in records if columns[2] == "transcript"]
                for transcript in transcripts:
                    cds = sorted((columns for columns, attrs in records
                                  if columns[2] == "CDS" and attrs["Parent"] == transcript),
                                 key=lambda c: int(c[3]), reverse=reverse)
                    count, coding = 0, ""
                    for columns in cds:
                        self.assertEqual(int(columns[7]), (3 - count % 3) % 3)
                        part = genome[int(columns[3]) - 1:int(columns[4])]
                        coding += part.translate(complement)[::-1] if reverse else part
                        count += len(part)
                    self.assertEqual(count % 3, 0)
                    self.assertTrue(coding.startswith("ATG") and coding.endswith("TAA"))
                    self.assertFalse(any(coding[i:i + 3] in ("TAA", "TAG", "TGA")
                                         for i in range(0, len(coding) - 3, 3)))

    def test_literal_scores_match_native_path_without_recalibration(self):
        output = self.decode("scores", check=True)
        path = native_path(self.data)
        emissions = [[struct.unpack("f", struct.pack("f", float(value)))[0]
                      for value in line.split()[1:6]]
                     for line in (self.data / "out.ps.txt").read_text().splitlines()]
        sparse = {name: {int(line.split()[0]): float(line.split()[1])
                         for line in (self.data / f"out.{name}.txt").read_text().splitlines()}
                  for name in ("gt", "ag", "atg", "stop")}
        references = [row for row in read_report(str(output) + ".tsv") if row["identical_to_reference"] == "1"]
        self.assertEqual(len(references), 3)
        for row in references:
            first, last = int(row["window_start_1based"]) - 1, int(row["window_end_1based"]) - 1
            expected = sum(emissions[i][min(path[i], 4)] for i in range(first, last + 1))
            for position in range(first + 1, last + 1):
                source, target = path[position - 1:position + 1]
                if source == target:
                    continue
                if source == 0 or target == 0:
                    kind, coordinate = ("atg" if source == 0 else "stop"), position - 2
                else:
                    kind, coordinate = ("gt" if source < 4 else "ag"), position - 1
                expected += sparse[kind][coordinate]
            self.assertAlmostEqual(float(row["fixed_interval_score"]), expected, places=9)

    def test_extra_base_from_unstripped_header_is_rejected(self):
        # An unstripped exporter header produces L+1 converter rows. AltAnn's
        # existing dense-input validation must reject the out-of-range last row.
        emissions = self.work / "extra-base.txt"
        emissions.write_text((self.data / "out.ps.txt").read_text() +
                             f"{self.metadata['length']}\t0\t0\t0\t0\t0\n")
        inputs = [self.data / "sequence.fna", emissions,
                  *(self.data / f"out.{name}.txt" for name in ("gt", "ag", "atg", "stop"))]
        result = subprocess.run([str(CORE), *map(str, inputs), "--output", str(self.work / "invalid.gff3"),
                                 "--report", str(self.work / "invalid.tsv")], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid or duplicate emission position", result.stderr)


if __name__ == "__main__":
    unittest.main()
