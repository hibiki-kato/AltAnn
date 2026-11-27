"""Exact search supports large K and preserves smaller-K raw path prefixes."""
from collections import defaultdict
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from fixtures import ordered_inputs, write_fixture
from test_altann import read_report

ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get('ALTANN_CORE', ROOT/'build/altann-core')).resolve()


@unittest.skipUnless(CORE.is_file(), 'Build altann-core first')
class MetadataStatesTests(unittest.TestCase):
    def test_k100_replays_and_preserves_k10_prefix(self):
        with tempfile.TemporaryDirectory(prefix='altann-exact-', dir='/tmp') as tmp:
            work = Path(tmp)
            inputs = write_fixture(work/'inputs')
            # More than ten frame-preserving acceptor choices exercise actual
            # higher ranks, in addition to the K=100 admission check.
            sequence = list(inputs['fasta'].read_text().splitlines()[1])
            acceptors = [p+offset for offset in (0,600) for p in range(238,299,3)]
            for p in acceptors:
                sequence[p:p+2] = 'AG'
            inputs['fasta'].write_text('>fixture\n'+''.join(sequence)+'\n')
            inputs['ag'].write_text(''.join(f'{p}\t2\n' for p in acceptors))
            results = {}
            for k in (10,100):
                output = work/f'k{k}.gff'
                report = work/f'k{k}.tsv'
                subprocess.run([str(CORE), *ordered_inputs(inputs), '--metadata-states',
                                '--verify-exact', '--raw', '--k',str(k),'--flank','50',
                                '--threads','2','--output',str(output),'--report',str(report)],
                               cwd=work, check=True, capture_output=True, text=True)
                results[k] = read_report(report)
            small = results[10]
            large_prefix = [r for r in results[100] if int(r['rank']) <= 10]
            self.assertTrue(small)
            self.assertEqual(small, large_prefix)
            self.assertTrue(any(int(r['rank']) > 10 for r in results[100]))
            by_locus = defaultdict(list)
            for r in results[100]:
                by_locus[r['reference_gene_index']].append(float(r['fixed_interval_score']))
            for scores in by_locus.values():
                self.assertEqual(scores, sorted(scores, reverse=True))
