"""Latest UniAnn coding-frame scores and native traces, independent of AltAnn."""
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

from test_altann import read_gff, read_report


ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get("ALTANN_CORE", ROOT / "build" / "altann-core")).resolve()
DATA = ROOT / "tests" / "data" / "coding-frame"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def float32(value):
    # Both UniAnn and AltAnn store the converter's emission scores as floats.
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def native_path():
    scores, predecessors = [], []
    for line in (DATA / "uniann.trace").read_text().splitlines():
        fields = line.split("\t")
        (scores if fields[1] == "dp" else predecessors).append(list(map(int, fields[2:])))
    state = max(range(7), key=lambda index: scores[-1][index])
    path = [0] * len(predecessors)
    for position in range(len(path) - 1, -1, -1):
        path[position] = state
        state = predecessors[position][state]
    return path


@unittest.skipUnless(CORE.is_file(), "Build altann-core before running integration tests")
class CodingFrameTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="altann-coding-")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.provenance = json.loads((DATA / "provenance.json").read_text())
        self.sequence = "".join((DATA / "sequence.fna").read_text().splitlines()[1:])
        self.emissions = [list(map(float32, line.split()[1:6]))
                          for line in (DATA / "out.ps.txt").read_text().splitlines()]

    def decode(self, name, reverse=False, saved=False, gff="uniann.gff3"):
        output = self.work / (name + ".gff3")
        command = [sys.executable, "-m", "altann", "decode", "--input", str(DATA),
                   "--core", str(CORE), "--k", "8", "--flank", "50", "--threads", "2",
                   "--output", str(output)]
        command += (["--gff", str(DATA / gff), "--log", str(DATA / "uniann.trace")]
                    if saved else ["--rerun-viterbi"])
        if reverse:
            command += ["--reverse"]
        subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True)
        return output

    def test_fixture_covers_local_frames_and_converter_stop_handling(self):
        raw = [line.split("\t") for line in (DATA / "coding.txt").read_text().splitlines()]
        converted = [line.split("\t") for line in (DATA / "out.ps.txt").read_text().splitlines()]
        self.assertEqual(len(raw), len(self.sequence))
        self.assertEqual(len(converted), len(self.sequence))
        self.assertEqual(len((DATA / "uniann.trace").read_text().splitlines()), 2 * len(self.sequence))
        self.assertEqual({g["start0"] % 3 for g in self.provenance["genes"]}, {0, 1, 2})
        self.assertEqual({g["intron_length"] % 3 for g in self.provenance["genes"]}, {0, 1, 2})
        for position, fields in enumerate(converted):
            self.assertEqual(len(fields), 7)
            self.assertEqual(int(fields[0]), position)
            self.assertEqual(fields[6], self.sequence[position])
            self.assertAlmostEqual(sum(map(float, raw[position][3:])), 1)
        for gene in self.provenance["genes"]:
            for position, frame in ((gene["start0"] + 8, gene["start0"] % 3),
                                    (gene["donor0"] + gene["intron_length"] + 8,
                                     (gene["start0"] + gene["intron_length"]) % 3)):
                probabilities = list(map(float, raw[position][3:6]))
                self.assertEqual(max(range(3), key=lambda i: probabilities[i]), (position - frame - 1) % 3)
                self.assertEqual(max(range(3), key=lambda i: self.emissions[position][i + 1]), frame)
            stop = gene["stop0"] + 2
            frame = gene["stop0"] % 3
            self.assertEqual(self.emissions[stop][frame + 1], -1000000)
            self.assertEqual(self.emissions[stop][4], -2)

    def test_native_trace_replay_matches_rerun_and_all_three_reference_frames(self):
        upstream = [c for c, _ in read_gff(DATA / "uniann.gff3") if c[2] == "CDS"]
        expected = [(a, b, str(phase)) for gene in self.provenance["genes"]
                    for a, b, phase in gene["cds_1based"]]
        self.assertEqual([(int(c[3]), int(c[4])) for c in upstream], [(a, b) for a, b, _ in expected])
        filtered = [(int(c[3]), int(c[4]), c[7]) for c, _ in read_gff(DATA / "uniann.filtered.gff3")
                    if c[2] == "CDS"]
        self.assertEqual(filtered, expected)
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                saved = self.decode("saved", reverse, saved=True)
                rerun = self.decode("rerun", reverse)
                filtered_saved = self.decode("filtered", reverse, saved=True, gff="uniann.filtered.gff3")
                for suffix in ("", ".tsv"):
                    self.assertEqual(Path(str(saved) + suffix).read_bytes(), Path(str(rerun) + suffix).read_bytes())
                    self.assertEqual(Path(str(filtered_saved) + suffix).read_bytes(), Path(str(rerun) + suffix).read_bytes())
                records = read_gff(saved)
                reference_ids = {attrs["ID"] for _, attrs in records if attrs.get("kbest_rank") == "0"}
                actual = sorted((int(c[3]), int(c[4]), c[7]) for c, attrs in records
                                if c[2] == "CDS" and attrs["Parent"] in reference_ids)
                oriented = ([(len(self.sequence) + 1 - b, len(self.sequence) + 1 - a, phase)
                             for a, b, phase in expected] if reverse else expected)
                self.assertEqual(actual, sorted(oriented))
                self.assertEqual({c[6] for c, _ in records}, {"-" if reverse else "+"})

    def test_scores_equal_native_path_sum_without_additional_transformation(self):
        output = self.decode("scores")
        path = native_path()
        sparse = {name: {int(line.split()[0]): float(line.split()[1])
                         for line in (DATA / f"out.{name}.txt").read_text().splitlines()}
                  for name in ("gt", "ag", "atg", "stop")}
        rows = read_report(str(output) + ".tsv")
        for flag, value in (("identical_to_reference", "1"), ("classification", "alternative_acceptor")):
            self.assertEqual({r["reference_gene_index"] for r in rows if r[flag] == value}, {"1", "2", "3"})
        for row in rows:
            gene = self.provenance["genes"][int(row["reference_gene_index"]) - 1]
            self.assertEqual(float(row["actual_start_score"]), sparse["atg"][gene["start0"]])
            if row["identical_to_reference"] == "1":
                first, last = int(row["window_start_1based"]) - 1, int(row["window_end_1based"]) - 1
                expected = sum(self.emissions[i][min(path[i], 4)] for i in range(first, last + 1))
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
            elif row["classification"] == "alternative_acceptor":
                acceptor = gene["donor0"] + gene["intron_length"] - 2
                frame = (gene["start0"] + gene["intron_length"]) % 3
                expected_delta = sparse["ag"][acceptor + 3] - sparse["ag"][acceptor]
                expected_delta += sum(self.emissions[i][4] - self.emissions[i][frame + 1]
                                      for i in range(acceptor + 1, acceptor + 4))
                self.assertAlmostEqual(float(row["score_delta"]), expected_delta, places=9)

    def test_alternative_acceptors_have_valid_cds_and_phase_on_both_strands(self):
        outputs = [self.decode("plus"), self.decode("minus", reverse=True)]
        for reverse, output in enumerate(outputs):
            rows = read_report(str(output) + ".tsv")
            alternatives = [r for r in rows if r["classification"] == "alternative_acceptor"]
            self.assertEqual({r["reference_gene_index"] for r in alternatives}, {"1", "2", "3"})
            genome = self.sequence.translate(COMPLEMENT)[::-1] if reverse else self.sequence
            records = read_gff(output)
            transcripts = [attrs["ID"] for c, attrs in records if c[2] == "transcript"]
            for transcript in transcripts:
                cds = sorted((c for c, attrs in records if c[2] == "CDS" and attrs["Parent"] == transcript),
                             key=lambda c: int(c[3]), reverse=bool(reverse))
                length, coding = 0, ""
                for c in cds:
                    self.assertEqual(int(c[7]), (3 - length % 3) % 3)
                    part = genome[int(c[3]) - 1:int(c[4])]
                    coding += part.translate(COMPLEMENT)[::-1] if reverse else part
                    length += len(part)
                self.assertEqual(length % 3, 0)
                self.assertTrue(coding.startswith("ATG"))
                self.assertTrue(coding.endswith("TAA"))
                self.assertFalse(any(coding[i:i + 3] in ("TAA", "TAG", "TGA")
                                     for i in range(0, len(coding) - 3, 3)))
        plus = [(c[2], len(self.sequence) + 1 - int(c[4]), len(self.sequence) + 1 - int(c[3]), c[7])
                for c, _ in read_gff(outputs[0])]
        minus = [(c[2], int(c[3]), int(c[4]), c[7]) for c, _ in read_gff(outputs[1])]
        self.assertEqual(sorted(plus), sorted(minus))


if __name__ == "__main__":
    unittest.main()
