"""--flank-cds-cap and --flank-junctions* set each locus's window from its own reference."""
import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from fixtures import ordered_inputs, write_fixture


ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get("ALTANN_CORE", ROOT / "build" / "altann-core")).resolve()
PAD = 3000  # noncoding bases added on both sides so flanks are not clipped by the sequence ends


def padded_fixture(directory):
    """The two-gene fixture with PAD noncoding bases before and after it."""
    paths = write_fixture(directory / "base")
    sequence = paths["fasta"].read_text().splitlines()[1]
    padded = {}
    padded["fasta"] = directory / "padded.fa"
    padded["fasta"].write_text(">fixture\n" + "C" * PAD + sequence + "C" * PAD + "\n")
    noncoding = "\t2\t-8\t-8\t-8\t-8\n"
    rows = paths["emissions"].read_text().splitlines()
    padded["emissions"] = directory / "emissions.txt"
    padded["emissions"].write_text(
        "".join(f"{p}{noncoding}" for p in range(PAD)) +
        "".join(f"{int(r.split(chr(9), 1)[0]) + PAD}\t{r.split(chr(9), 1)[1]}\n" for r in rows) +
        "".join(f"{p}{noncoding}" for p in range(PAD + len(sequence), 2 * PAD + len(sequence))))
    for key in ("gt", "ag", "atg", "stop"):
        padded[key] = directory / f"{key}.txt"
        padded[key].write_text("".join(f"{int(r.split()[0]) + PAD}\t{r.split()[1]}\n"
                                       for r in paths[key].read_text().splitlines()))
    return padded


@unittest.skipUnless(CORE.is_file(), "Build altann-core before running integration tests")
class FlankOptionTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="altann-flank-")
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.inputs = padded_fixture(self.work)

    def run_core(self, name, *options, check=True):
        output, report = self.work / f"{name}.gff3", self.work / f"{name}.tsv"
        result = subprocess.run([str(CORE), *ordered_inputs(self.inputs), "--k", "2", "--threads", "1",
                                 "--output", str(output), "--report", str(report), *map(str, options)],
                                capture_output=True, text=True, check=check)
        return result, output, report

    def windows(self, *options):
        """Per locus: reference start/end, window start/end (1-based) and reference CDS intervals."""
        _, output, report = self.run_core("run", *options)
        loci = {}
        with report.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                loci[int(row["reference_gene_index"])] = tuple(int(row[key]) for key in (
                    "reference_start_1based", "reference_end_1based", "window_start_1based", "window_end_1based"))
        cds = {}
        for line in Path(str(output) + ".reference.gff").read_text().splitlines():
            fields = line.split("\t")
            if len(fields) == 9 and fields[2] == "CDS":
                cds.setdefault(fields[8].split("Parent=")[1].split(";")[0], []).append((int(fields[3]), int(fields[4])))
        intervals = [sorted(v) for _, v in sorted(cds.items(), key=lambda item: min(item[1]))]
        return [loci[index] for index in sorted(loci)], intervals

    def test_cds_cap_limits_each_outer_side_to_the_reference_coding_length(self):
        plain, cds = self.windows("--flank", 1000)
        capped, _ = self.windows("--flank", 1000, "--flank-cds-cap")
        coding = [sum(end - start + 1 for start, end in gene) for gene in cds]
        self.assertTrue(all(c < 1000 for c in coding))
        self.assertEqual(plain[0][2], plain[0][0] - 1000)
        self.assertEqual(plain[-1][3], plain[-1][1] + 1000)
        self.assertEqual(capped[0][2], capped[0][0] - coding[0])
        self.assertEqual(capped[-1][3], capped[-1][1] + coding[-1])

    def test_junction_leaving_the_reference_widens_that_side(self):
        start, end, _, _ = self.windows("--flank", 100)[0][-1]
        junction_end = end + 1500  # zero-based exclusive end of [GT, AG + 2)
        table = self.work / "junctions.tsv"
        table.write_text(f"{start + 50}\t{junction_end}\n")
        widened, _ = self.windows("--flank", 100, "--flank-junctions", table, "--flank-junction-pad", 200)
        self.assertEqual(widened[-1][3], junction_end + 200 + 1)
        capped, _ = self.windows("--flank", 100, "--flank-junctions", table, "--flank-junction-pad", 200,
                                 "--flank-junction-cap", 500)
        self.assertEqual(capped[-1][3], end + 500)
        self.assertEqual(widened[0], capped[0])

    def test_read_fraction_drops_weak_junctions(self):
        loci, cds = self.windows("--flank", 100)
        start, end = loci[-1][:2]
        own = (cds[-1][0][1], cds[-1][1][0] - 1)  # the reference's own intron, [GT, AG + 2)
        table = self.work / "junctions.tsv"
        table.write_text(f"{own[0]}\t{own[1]}\t100\n{start + 50}\t{end + 1500}\t10\n")
        strict, _ = self.windows("--flank", 100, "--flank-junctions", table, "--flank-junction-min-fraction", 0.5)
        loose, _ = self.windows("--flank", 100, "--flank-junctions", table, "--flank-junction-min-fraction", 0.05)
        self.assertEqual(strict[-1][3], end + 100)
        self.assertEqual(loose[-1][3], end + 1500 + 1000 + 1)

    def test_invalid_settings_are_rejected(self):
        table = self.work / "plain_junctions.txt"
        table.write_text("10\t100\n")
        mixed = self.work / "mixed_counts.txt"
        mixed.write_text("10\t100\t5\n20\t200\n")
        for name, options, message in (
                ("cap", ["--flank", -1, "--flank-cds-cap"], "--flank-cds-cap needs a nonnegative --flank"),
                ("counts", ["--flank", 100, "--flank-junctions", table, "--flank-junction-min-fraction", 0.5],
                 "needs read counts"),
                ("mixed", ["--flank", 100, "--flank-junctions", mixed], "every row or none")):
            with self.subTest(name=name):
                result, output, _ = self.run_core(name, *options, check=False)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
