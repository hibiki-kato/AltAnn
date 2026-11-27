"""Keep forced-junction decoding compatible with the windowed score source."""

import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from fixtures import ordered_inputs, write_fixture


ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get("ALTANN_CORE", ROOT / "build" / "altann-core")).resolve()


@unittest.skipUnless(CORE.is_file(), "Build altann-core before running integration tests")
class ForcedJunctionTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="altann-forcing-")
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.inputs = write_fixture(self.work / "inputs")
        self.table = self.work / "introns.tsv"
        self.table.write_text("first\t180-240\nsecond\t780-840\n")

    def run_forcing(self, name, *options):
        output = self.work / f"{name}.tsv"
        result = subprocess.run(
            [str(CORE), *ordered_inputs(self.inputs), "--force-junctions", str(self.table),
             "--output", str(output), "--flank", "50", *options],
            capture_output=True, text=True, cwd=self.work,
        )
        return result, output

    def test_forced_introns_recover_both_reference_cds_and_scores(self):
        result, output = self.run_forcing("reference", "--threads", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual([row["cds"] for row in rows],
                         ["60-180,240-453", "660-780,840-1053"])
        for row in rows:
            self.assertEqual(row["status"], "ok")
            self.assertEqual(row["identical_to_reference"], "1")
            self.assertEqual(float(row["delta"]), 0)
            self.assertEqual(row["model_score"], row["reference_score"])

    def test_spilled_parallel_search_matches_in_memory_serial_search(self):
        serial, memory = self.run_forcing("memory", "--threads", "1")
        parallel, spilled = self.run_forcing(
            "spilled", "--threads", "2", "--traceback-memory", "0")
        self.assertEqual(serial.returncode, 0, serial.stderr)
        self.assertEqual(parallel.returncode, 0, parallel.stderr)
        self.assertEqual(memory.read_bytes(), spilled.read_bytes())

    def test_other_output_modes_cannot_silently_override_forcing(self):
        for mode in ("--raw", "--uniann-gff"):
            with self.subTest(mode=mode):
                result, output = self.run_forcing("invalid", mode)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists())
