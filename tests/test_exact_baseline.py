"""--exact-baseline: the global path from exact K=1 metadata-state search, without UniAnn."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get("ALTANN_CORE", ROOT / "build" / "altann-core")).resolve()
DATA = ROOT / "tests" / "data" / "coding-frame"
INPUTS = [DATA / name for name in ("sequence.fna", "out.ps.txt", "out.gt.txt", "out.ag.txt",
                                   "out.atg.txt", "out.stop.txt")]


@unittest.skipUnless(CORE.is_file(), "Build altann-core before running integration tests")
class ExactBaselineTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="altann-exact-baseline-")
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.genes = [gene["cds_1based"] for gene in json.loads((DATA / "provenance.json").read_text())["genes"]]

    def core(self, *options, check=True):
        return subprocess.run([str(CORE), *map(str, INPUTS), *options], capture_output=True, text=True, check=check)

    def test_exact_global_path_matches_the_designed_genes_and_the_coarse_path(self):
        coarse = self.core("--uniann-gff")
        exact = self.core("--uniann-gff", "--exact-baseline")
        self.assertEqual(exact.stdout, coarse.stdout)
        self.assertIn("exact_baseline\tpositions=1740", exact.stderr)
        cds = {}
        for line in exact.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) == 9 and fields[2] == "CDS":
                parent = fields[8].split("Parent=")[1].split(";")[0]
                cds.setdefault(parent, []).append([int(fields[3]), int(fields[4])])  # UniAnn format: no phase
        self.assertEqual(sorted(cds.values()), [[interval[:2] for interval in gene] for gene in self.genes])

    def test_cli_needs_no_uniann_gff_or_trace(self):
        output = self.work / "exact.gff3"
        subprocess.run([sys.executable, "-m", "altann", "decode", "--input", str(DATA), "--core", str(CORE),
                        "--exact-baseline", "--k", "4", "--flank", "50", "--threads", "2",
                        "--output", str(output)], cwd=ROOT, capture_output=True, text=True, check=True)
        metadata = json.loads(Path(str(output) + ".json").read_text())
        self.assertEqual(metadata["reference_mode"], "exact")
        self.assertNotIn("best_gff", metadata["jobs"][0]["inputs"])
        self.assertNotIn("log", metadata["jobs"][0]["inputs"])
        reference = {}
        ranks = {}
        for line in output.read_text().splitlines():
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            attributes = dict(item.split("=", 1) for item in fields[8].split(";"))
            if fields[2] == "transcript":
                ranks[attributes["ID"]] = attributes["kbest_rank"]
            elif fields[2] == "CDS" and ranks[attributes["Parent"]] == "0":
                reference.setdefault(attributes["Parent"], []).append([int(fields[3]), int(fields[4]), int(fields[7])])
        self.assertEqual(sorted(reference.values()), self.genes)

    def test_saved_traces_are_rejected(self):
        result = self.core("--exact-baseline", "--viterbi-log", DATA / "uniann.spliced.trace",
                           "--output", self.work / "x.gff3", "--report", self.work / "x.tsv", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--exact-baseline cannot be combined with --viterbi-log", result.stderr)
        result = subprocess.run([sys.executable, "-m", "altann", "decode", "--input", str(DATA), "--core", str(CORE),
                                 "--exact-baseline", "--log", str(DATA / "uniann.spliced.trace"),
                                 "--output", str(self.work / "y.gff3")], cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--log cannot be combined with --exact-baseline", result.stderr)


if __name__ == "__main__":
    unittest.main()
