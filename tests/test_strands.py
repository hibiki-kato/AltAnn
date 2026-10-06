"""Compare one both-strand invocation to independently oriented decoder jobs."""
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from fixtures import write_stranded_fixture
from test_altann import read_gff, read_report
from altann.strands import COMPLEMENT

ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get('ALTANN_CORE', ROOT / 'build' / 'altann-core')).resolve()


@unittest.skipUnless(CORE.is_file(), 'Build altann-core before running integration tests')
class BothStrandsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='altann-strands-')
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.combined, self.plus, self.minus = write_stranded_fixture(self.work / 'inputs')

    def run_cli(self, inputs, name, both=False, reverse=False, extra=(), check=True):
        output = self.work / (name + '.gff3')
        command = [sys.executable, '-m', 'altann', 'decode', '--core', str(CORE),
                   '--k', '8', '--flank', '50', '--threads', '2', '--output', str(output)]
        for key, value in inputs.items():
            command.extend(['--' + key, str(value)])
        command += ['--all-prob'] if both else ['--rerun-viterbi']
        if reverse:
            command += ['--reverse']
        result = subprocess.run([*command, *extra], cwd=ROOT, text=True,
                                capture_output=True, check=check)
        return result, output

    def test_matches_two_legacy_runs_for_each_length_modulo_three(self):
        for length in (1200, 1201, 1202):
            with self.subTest(length=length):
                combined, plus, minus = write_stranded_fixture(self.work / str(length), length)
                # Tables need not be sorted; minus positions run in the
                # opposite order to their oriented decoder coordinates.
                for key in ('emissions', 'gt', 'ag', 'atg', 'stop'):
                    rows = combined[key].read_text().splitlines()
                    random.Random(41).shuffle(rows)
                    combined[key].write_text('\n'.join(rows) + '\n')
                _, forward = self.run_cli(plus, f'plus-{length}')
                _, reverse = self.run_cli(minus, f'minus-{length}', reverse=True)
                _, both = self.run_cli(combined, f'both-{length}', both=True)
                expected = read_gff(forward) + read_gff(reverse)
                actual = read_gff(both)
                self.assertEqual(sorted(tuple(c) for c, _ in expected),
                                 sorted(tuple(c) for c, _ in actual))
                self.assertEqual({c[6] for c, _ in actual}, {'+', '-'})
                self.assertTrue(any(c[2] == 'CDS' and c[6] == '-' and c[7] in ('1', '2')
                                    for c, _ in actual))
                self.assertEqual(read_report(str(both) + '.tsv'),
                                 read_report(str(forward) + '.tsv') + read_report(str(reverse) + '.tsv'))
                ids = [a['ID'] for _, a in actual]
                self.assertEqual(len(ids), len(set(ids)))
                self.assertTrue(all(a.get('Parent', ids[0]) in ids for _, a in actual))

    def test_two_jobs_record_original_inputs_and_one_tsv_header(self):
        _, output = self.run_cli(self.combined, 'metadata', both=True)
        data = json.loads(Path(str(output) + '.json').read_text())
        self.assertTrue(data['both_strands'])
        self.assertEqual(data['reference_mode'], 'rerun')
        self.assertEqual([j['inputs']['strand'] for j in data['jobs']], ['plus', 'minus'])
        for job in data['jobs']:
            for key, path in self.combined.items():
                self.assertEqual(job['inputs'][key], str(path.resolve()))
        report = Path(str(output) + '.tsv').read_text().splitlines()
        self.assertEqual(sum(line.startswith('segment\tstrand\t') for line in report), 1)

    def test_input_directory_discovery(self):
        _, explicit = self.run_cli(self.combined, 'explicit', both=True)
        _, discovered = self.run_cli({'input': self.combined['fasta'].parent}, 'discovered', both=True)
        self.assertEqual(explicit.read_bytes(), discovered.read_bytes())

    def test_short_option_and_whitespace_in_segment_fasta(self):
        _, baseline = self.run_cli(self.combined, 'baseline', both=True)
        fasta = self.combined['fasta']
        sequence = ''.join(fasta.read_text().splitlines()[1:])
        fasta.write_text('>fixture__0-1200-1200\n' +
                         '\n'.join(' '.join(sequence[i:i+60]) for i in range(0, 1200, 60)) + '\n')
        _, output = self.run_cli(self.combined, 'whitespace', both=True, extra=['-a'])
        expected = baseline.read_text().replace('fixture.plus', 'fixture__0-1200-1200.plus').replace(
            'fixture.minus', 'fixture__0-1200-1200.minus')
        self.assertEqual(output.read_text(), expected)

    def test_combined_best_gff_and_segment_mapping(self):
        extra = ['--offset', '1000', '--sequence-length', '3200', '--seqid', 'chr',
                 '--segment-overlap', '0', '--segment-margin', '0']
        _, output = self.run_cli(self.combined, 'segment', both=True, extra=extra)
        rows = read_gff(output)
        ids = {a['ID'] for _, a in rows if a.get('kbest_rank') == '0'}
        best = self.work / 'best.gff'
        best.write_text('##gff-version 3\n' + ''.join('\t'.join(c) + '\n' for c, a in rows
                         if a.get('ID') in ids or a.get('Parent') in ids))
        _, checked = self.run_cli(self.combined, 'checked', both=True,
                                  extra=[*extra, '--gff', str(best)])
        self.assertEqual(output.read_bytes(), checked.read_bytes())

    def test_minus_scores_independently_control_minus_prediction(self):
        for key in ('gt', 'ag', 'atg', 'stop'):
            path = self.combined[key]
            path.write_text(''.join(line for line in path.read_text().splitlines(keepends=True)
                                    if line.startswith('#') or line.split()[1] == '+'))
        _, output = self.run_cli(self.combined, 'no-minus', both=True)
        self.assertEqual({c[6] for c, _ in read_gff(output)}, {'+'})
        _, baseline = self.run_cli(self.plus, 'plus-baseline')
        self.assertEqual(output.read_bytes(), baseline.read_bytes())

    def test_thread_counts_are_byte_identical(self):
        _, serial = self.run_cli(self.combined, 'serial', both=True, extra=['--threads', '1'])
        _, parallel = self.run_cli(self.combined, 'parallel', both=True, extra=['--threads', '4'])
        for suffix in ('', '.tsv'):
            self.assertEqual(Path(str(serial) + suffix).read_bytes(), Path(str(parallel) + suffix).read_bytes())

    def test_malformed_inputs_preserve_all_outputs(self):
        path = self.combined['emissions']
        original = path.read_text()
        rows = original.splitlines()
        cases = [('\n'.join(rows[:-1]) + '\n', 'Missing - emission'),
                 (original + rows[1] + '\n', 'duplicate'),
                 (original.replace('0\t+\t', '0\t?\t', 1), 'Expected position'),
                 (original.replace('0\t+\t', '1200\t+\t', 1), 'Invalid'),
                 (original.replace('0\t+\t2.0', '0\t+\tnan', 1), 'finite'),
                 (original.replace('0\t+\t', '0\t', 1), 'Expected position')]
        for data, message in cases:
            with self.subTest(message=message):
                path.write_text(data)
                output = self.work / 'invalid.gff3'
                for suffix in ('', '.tsv', '.json'):
                    Path(str(output) + suffix).write_text('previous output\n')
                result, _ = self.run_cli(self.combined, 'invalid', both=True, check=False)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                for suffix in ('', '.tsv', '.json'):
                    self.assertEqual(Path(str(output) + suffix).read_text(), 'previous output\n')

    def test_incompatible_trace_and_reverse_options(self):
        for extra in (['--reverse'], ['--log', str(ROOT / 'tests/data/uniann-fixture.trace')]):
            result, output = self.run_cli(self.combined, 'incompatible', both=True, extra=extra, check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('cannot be combined with --all-prob', result.stderr)
            self.assertFalse(output.exists())

    def test_iupac_complement_is_an_involution(self):
        sequence = 'ACGTRYSWKMBDHVNacgtryswkmbdhvn'
        self.assertEqual(sequence.translate(COMPLEMENT).translate(COMPLEMENT), sequence)
        self.assertEqual(sequence.translate(COMPLEMENT), 'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')

    def test_overwrite_guard_includes_stranded_inputs(self):
        original = self.combined['emissions'].read_bytes()
        result, _ = self.run_cli(self.combined, 'overwrite', both=True,
                                  extra=['--output', str(self.combined['emissions'])], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('must not overwrite input', result.stderr)
        self.assertEqual(self.combined['emissions'].read_bytes(), original)
