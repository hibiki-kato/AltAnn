"""`altann prepare` against upstream UniAnn f126cc8 (no temperature) outputs."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" / "cnn-input"
OUTPUTS = ("out.ps.txt", "out.gt.txt", "out.ag.txt", "out.atg.txt", "out.stop.txt")


class PrepareCnnTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="altann-prepare-")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)

    def prepare(self, output, coding=DATA / "coding.tsv", scores=DATA / "sites.tsv", check=True):
        command = [sys.executable, "-m", "altann", "prepare", "--fasta", str(DATA / "sequence.fna"),
                   "--coding-frame", str(coding), "--output", str(output)]
        if scores is not None:
            command += ["--scores", str(scores)]
        return subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=check)

    def test_matches_upstream_without_temperature(self):
        output = self.work / "prepared"
        self.prepare(output)
        for name in OUTPUTS:
            with self.subTest(name=name):
                self.assertEqual((output / name).read_bytes(), (DATA / name).read_bytes())
        metadata = json.loads((output / "prepare.json").read_text())
        self.assertIsNone(metadata["temperature"])
        self.assertEqual(metadata["factor"], 1003)
        self.assertEqual(metadata["seqid"], "cnn_fixture")

    def test_skips_rows_of_other_sequences(self):
        # Upstream uniann.sh would score this row on the wrong sequence.
        scores = self.work / "sites.tsv"
        scores.write_text((DATA / "sites.tsv").read_text() + "other_sequence\t105\t+\tdonor\tGT\t0.5\n")
        output = self.work / "other"
        self.prepare(output, scores=scores)
        for name in OUTPUTS[1:]:
            with self.subTest(name=name):
                self.assertEqual((output / name).read_bytes(), (DATA / name).read_bytes())

    def test_coding_only(self):
        output = self.work / "coding-only"
        self.prepare(output, scores=None)
        self.assertEqual((output / "out.ps.txt").read_bytes(), (DATA / "out.ps.txt").read_bytes())
        self.assertFalse((output / "out.gt.txt").exists())

    def test_header_is_optional(self):
        coding = self.work / "headerless.tsv"
        coding.write_text("".join((DATA / "coding.tsv").read_text().splitlines(True)[1:]))
        output = self.work / "headerless"
        self.prepare(output, coding=coding, scores=None)
        self.assertEqual((output / "out.ps.txt").read_bytes(), (DATA / "out.ps.txt").read_bytes())

    def test_rejects_mismatched_tables(self):
        lines = (DATA / "coding.tsv").read_text().splitlines(True)
        cases = {
            "minus-label": [lines[0], lines[1].replace("\t+\t", "\t-\t", 1), *lines[2:]],
            "gap": [lines[0], *lines[2:]],
            "short": lines[:-1],
            "long": [*lines, lines[-1].replace("\t1500\t", "\t1501\t", 1)],
            "seqid": [lines[0], lines[1].replace("cnn_fixture", "other", 1), *lines[2:]],
            "probability": [lines[0], "cnn_fixture\t1\t+\t1.5\t0\t0\t0\tframe0\n", *lines[2:]],
        }
        for name, rows in cases.items():
            with self.subTest(name=name):
                coding = self.work / f"{name}.tsv"
                coding.write_text("".join(rows))
                result = self.prepare(self.work / name, coding=coding, scores=None, check=False)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn("altann: error:", result.stderr)

    def test_rejects_nonempty_output(self):
        output = self.work / "used"
        output.mkdir()
        (output / "keep.txt").write_text("x")
        result = self.prepare(output, check=False)
        self.assertEqual(result.returncode, 2)
        self.assertEqual((output / "keep.txt").read_text(), "x")


if __name__ == "__main__":
    unittest.main()
