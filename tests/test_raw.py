"""Decode original-coordinate site scores and six-frame PSAURON probabilities."""

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from fixtures import write_raw_fixture
from test_altann import read_gff, read_report


ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get('ALTANN_CORE', ROOT / 'build' / 'altann-core')).resolve()


@unittest.skipUnless(CORE.is_file(), 'Build altann-core before running integration tests')
class RawInputsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='altann-raw-')
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.combined, self.plus, self.minus = write_raw_fixture(self.work / 'inputs')

    def run_cli(self, inputs, name, both=False, reverse=False, extra=(), check=True):
        output = self.work / (name + '.gff3')
        command = [sys.executable, '-m', 'altann', 'decode', '--core', str(CORE),
                   '--k', '8', '--flank', '50', '--threads', '2', '--output', str(output)]
        for key, value in inputs.items():
            command.extend(['--' + key, str(value)])
        if both:
            command += ['-a']
        if reverse:
            command += ['--reverse']
        result = subprocess.run([*command, *extra], cwd=ROOT, text=True,
                                capture_output=True, check=check)
        return result, output

    def test_matches_independent_oriented_jobs_for_each_length_modulo_three(self):
        for length in (1400, 1401, 1402):
            for columns in (6, 7):
                with self.subTest(length=length, columns=columns):
                    combined, plus, minus = write_raw_fixture(self.work / f'{length}-{columns}',
                                                              length, columns)
                    _, forward = self.run_cli(plus, f'plus-{length}-{columns}')
                    _, reverse = self.run_cli(minus, f'minus-{length}-{columns}', reverse=True)
                    _, both = self.run_cli(combined, f'both-{length}-{columns}', both=True)
                    expected = read_gff(forward) + read_gff(reverse)
                    actual = read_gff(both)
                    self.assertEqual(sorted(tuple(c) for c, _ in expected),
                                     sorted(tuple(c) for c, _ in actual))
                    self.assertEqual({c[6] for c, _ in actual}, {'+', '-'})
                    self.assertTrue(any(c[2] == 'CDS' and c[6] == '-' and c[7] in ('1', '2')
                                        for c, _ in actual))
                    self.assertEqual(read_report(str(both) + '.tsv'),
                                     read_report(str(forward) + '.tsv') +
                                     read_report(str(reverse) + '.tsv'))
                    ids = {a['ID'] for _, a in actual}
                    self.assertEqual(len(ids), len(actual))
                    self.assertTrue(all(a.get('Parent', next(iter(ids))) in ids for _, a in actual))

    def test_input_directory_discovery_and_short_options(self):
        _, explicit = self.run_cli(self.combined, 'explicit', both=True)
        _, discovered = self.run_cli({'input': self.combined['fasta'].parent}, 'discovered', both=True)
        self.assertEqual(explicit.read_bytes(), discovered.read_bytes())
        _, short = self.run_cli({}, 'short', both=True,
                                extra=['-f', str(self.combined['fasta']),
                                       '-p', str(self.combined['psauron']),
                                       '-s', str(self.combined['scores']), '-m', '2.71828182845905'])
        self.assertEqual(explicit.read_bytes(), short.read_bytes())
        self.combined['scores'].rename(self.combined['scores'].with_name('fixture_sites.tsv'))
        _, renamed = self.run_cli({'input': self.combined['fasta'].parent}, 'renamed', both=True)
        self.assertEqual(explicit.read_bytes(), renamed.read_bytes())

    def test_metadata_keeps_original_input_paths_and_one_report_header(self):
        _, output = self.run_cli(self.combined, 'metadata', both=True)
        data = json.loads(Path(str(output) + '.json').read_text())
        self.assertTrue(data['both_strands'])
        self.assertEqual(data['reference_mode'], 'rerun')
        self.assertEqual([job['inputs']['strand'] for job in data['jobs']], ['plus', 'minus'])
        for job in data['jobs']:
            for key, path in self.combined.items():
                self.assertEqual(job['inputs'][key], str(path.resolve()))
        report = Path(str(output) + '.tsv').read_text().splitlines()
        self.assertEqual(sum(line.startswith('segment\tstrand\t') for line in report), 1)

    def test_unrelated_sequence_rows_do_not_change_scores(self):
        _, baseline = self.run_cli(self.combined, 'baseline', both=True)
        with self.combined['scores'].open('a') as stream:
            stream.write('other\t999999\t-\tdonor\tGT\t.99\n')
        _, output = self.run_cli(self.combined, 'unrelated', both=True)
        self.assertEqual(output.read_bytes(), baseline.read_bytes())
        self.assertEqual(Path(str(output) + '.tsv').read_bytes(),
                         Path(str(baseline) + '.tsv').read_bytes())

    def test_donorless_minus_is_skipped(self):
        path = self.combined['scores']
        path.write_text(''.join(line for line in path.read_text().splitlines(keepends=True)
                                if not ('\t-\tdonor\t' in line)))
        _, output = self.run_cli(self.combined, 'no-minus', both=True)
        _, forward = self.run_cli(self.plus, 'plus')
        self.assertEqual(output.read_bytes(), forward.read_bytes())
        self.assertEqual(Path(str(output) + '.tsv').read_bytes(),
                         Path(str(forward) + '.tsv').read_bytes())

    def test_supplied_gff_cannot_add_models_on_a_skipped_strand(self):
        path = self.combined['scores']
        path.write_text(''.join(line for line in path.read_text().splitlines(keepends=True)
                                if not ('\t-\tdonor\t' in line)))
        best = self.work / 'best.gff3'
        best.write_text('fixture\tUniAnn\ttranscript\t5\t250\t.\t-\t.\tID=wrong\n'
                        'fixture\tUniAnn\texon\t5\t250\t.\t-\t.\tID=e;Parent=wrong\n')
        output = self.work / 'invalid-best.gff3'
        output.write_text('previous output\n')
        result, _ = self.run_cli(self.combined, 'invalid-best', both=True,
                                 extra=['--gff', str(best)], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('features for skipped - strand', result.stderr)
        self.assertEqual(output.read_text(), 'previous output\n')

    def test_both_strands_without_positive_rescaled_donors_produce_empty_outputs(self):
        _, output = self.run_cli(self.combined, 'empty', both=True, extra=['--mult', '1'])
        self.assertEqual(read_gff(output), [])
        self.assertEqual(read_report(str(output) + '.tsv'), [])
        metadata = json.loads(Path(str(output) + '.json').read_text())
        self.assertEqual(metadata['jobs'], [])

    def test_donorless_strand_still_requires_valid_reverse_probabilities(self):
        path = self.combined['scores']
        path.write_text(''.join(line for line in path.read_text().splitlines(keepends=True)
                                if not ('\t-\tdonor\t' in line)))
        path = self.combined['psauron']
        lines = path.read_text().splitlines()
        fields = next(csv.reader([lines[4]]))
        fields[14] = 'nan'
        path.write_text('\n'.join(lines[:4]) + '\n' + ','.join(fields) + '\n')
        result, output = self.run_cli(self.combined, 'invalid-skipped', both=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(output.exists())

    def test_reverse_probabilities_control_the_minus_prediction(self):
        _, baseline = self.run_cli(self.combined, 'baseline', both=True)
        path = self.combined['psauron']
        lines = path.read_text().splitlines()
        fields = next(csv.reader([lines[4]]))
        for index in (12, 13, 14):
            fields[index] = ';'.join('0' for _ in fields[index].split(';'))
        with path.open('w', newline='') as stream:
            stream.write('\n'.join(lines[:4]) + '\n')
            csv.writer(stream).writerow(fields)
        _, output = self.run_cli(self.combined, 'no-reverse-coding', both=True)
        original = read_gff(baseline)
        actual = read_gff(output)
        self.assertEqual([c for c, _ in original if c[6] == '+'],
                         [c for c, _ in actual if c[6] == '+'])
        self.assertNotEqual([c for c, _ in original if c[6] == '-'],
                            [c for c, _ in actual if c[6] == '-'])

    def test_malformed_inputs_leave_existing_outputs_untouched(self):
        csv_path, site_path = self.combined['psauron'], self.combined['scores']
        original_csv, original_sites = csv_path.read_text(), site_path.read_text()
        csv_lines = original_csv.splitlines()
        csv_fields = next(csv.reader([csv_lines[4]]))
        cases = []
        for value in ('', 'nan', '1.1', '-.1', 'abc'):
            cells = list(csv_fields)
            cells[14] = value
            cases.append((csv_path, '\n'.join(csv_lines[:4]) + '\n' + ','.join(cells) + '\n'))
        cases.append((csv_path, '\n'.join(csv_lines[:4]) + '\n' +
                      ','.join(csv_fields[:12]) + '\n'))
        cases.append((csv_path, original_csv.replace('fixture,', 'unmatched,', 1)))
        cases.append((csv_path, original_csv + csv_lines[4] + '\n'))
        for position, strand, probability in (('0', '+', '.9'), ('1401', '-', '.9'),
                                               ('4.2', '+', '.9'), ('7', '+', 'nan'),
                                               ('7', '+', '1.1'), ('7', '+', '-.1')):
            cases.append((site_path, original_sites +
                          f'fixture\t{position}\t{strand}\tdonor\tGT\t{probability}\n'))
        for index, (path, data) in enumerate(cases):
            with self.subTest(case=index):
                csv_path.write_text(original_csv)
                site_path.write_text(original_sites)
                path.write_text(data)
                output = self.work / 'invalid.gff3'
                for suffix in ('', '.tsv', '.json'):
                    Path(str(output) + suffix).write_text('previous output\n')
                result, _ = self.run_cli(self.combined, 'invalid', both=True, check=False)
                self.assertNotEqual(result.returncode, 0)
                for suffix in ('', '.tsv', '.json'):
                    self.assertEqual(Path(str(output) + suffix).read_text(), 'previous output\n')

    def test_raw_input_overwrite_guard(self):
        for key in ('psauron', 'scores', 'fasta'):
            with self.subTest(key=key):
                path = self.combined[key]
                original = path.read_bytes()
                result, _ = self.run_cli(self.combined, 'overwrite', both=True,
                                          extra=['--output', str(path)], check=False)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('must not overwrite input', result.stderr)
                self.assertEqual(path.read_bytes(), original)
