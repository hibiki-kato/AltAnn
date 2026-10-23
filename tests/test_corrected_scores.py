"""Corrected upstream compatibility using frozen, synthetic-logit model export."""
import json
import math
from pathlib import Path
import unittest

import test_bowen_contract as contract


class CorrectedBowenProcessedScoreTests(contract.BowenProcessedScoreTests):
    data = Path(__file__).resolve().parent / "data/bowen-corrected"

    def test_header_produces_exactly_one_dense_emission_per_base(self):
        raw = (self.data / "exported.tsv").read_text().splitlines()
        emissions = (self.data / "out.ps.txt").read_text().splitlines()
        self.assertEqual(len(raw), self.metadata["length"] + 1)
        self.assertEqual(len(emissions), self.metadata["length"])
        self.assertEqual([int(line.split()[0]) for line in emissions],
                         list(range(self.metadata["length"])))
        self.assertEqual(raw[0].split("\t")[1], "position_1based")
        self.assertTrue(self.metadata["synthetic_logits"])
        self.assertFalse(self.metadata["checkpoint_inference"])

    def test_emission_peak_matches_local_codon_base_without_rotation(self):
        emissions = [list(map(float, line.split()[1:6]))
                     for line in (self.data / "out.ps.txt").read_text().splitlines()]
        archived = [list(map(float, line.split()[1:6]))
                    for line in (contract.DATA / "out.ps.txt").read_text().splitlines()]
        # Interior codons avoid the intentional stop mask at the final stop base.
        frames = set()
        for gene in self.metadata["genes"]:
            for start, end, phase in gene["cds_1based"]:
                for position in range(start + 7, end - 5):
                    local = int(self.metadata["phase_labels"][position])
                    frame = (position - local) % 3
                    frames.add(frame)
                    self.assertEqual(max(range(3), key=lambda i: emissions[position][i + 1]), frame)
                    self.assertEqual(max(range(3), key=lambda i: archived[position][i + 1]), (frame - 1) % 3)
        self.assertEqual(frames, {0, 1, 2})

    def test_noncoding_scores_use_exported_probabilities_without_temperature7(self):
        raw = [line.split("\t") for line in (self.data / "exported.tsv").read_text().splitlines()[1:]]
        emissions = [list(map(float, line.split()[1:6]))
                     for line in (self.data / "out.ps.txt").read_text().splitlines()]
        for fields, scores in zip(raw, emissions):
            nc = float(fields[6])
            self.assertEqual(scores[0], round(math.log(math.e * nc + 1e-10), 2))
        position = self.metadata["genes"][0]["start0"] + 8
        nc = float(raw[position][6])
        scaled = 1 / (1 + ((1 - nc) / (nc + 1e-10)) ** (1 / 7))
        self.assertGreater(abs(emissions[position][0] - math.log(math.e * scaled + 1e-10)), 1)
        self.assertFalse(self.metadata["additional_temperature_scaling"])

    def test_corrected_fixture_retains_identical_frozen_model_export(self):
        import hashlib
        for name, expected in self.metadata["input_sha256"].items():
            value = (self.data / name).read_bytes()
            self.assertEqual(value, (contract.DATA / name).read_bytes())
            self.assertEqual(hashlib.sha256(value).hexdigest(), expected)
        self.assertEqual(self.metadata["uniann_base_commit"],
                         json.loads((contract.DATA / "provenance.json").read_text())["uniann_commit"])


if __name__ == "__main__":
    unittest.main()
