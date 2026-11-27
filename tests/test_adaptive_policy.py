import json
import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from altann import cli
from altann.policy import load_adaptive_k_policy
from fixtures import write_fixture

CORE = Path(os.environ.get('ALTANN_CORE', Path(__file__).resolve().parents[1]/'build/altann-core')).resolve()


class AdaptivePolicyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir='/tmp')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.path = self.root / 'policy.json'
        self.linear = {'schema': 'dna-adaptive-k-v1', 'family': 'linear',
                       'base': 4, 'slope': 3, 'cap': 40}

    def load(self, policy):
        self.path.write_text(json.dumps(policy))
        return load_adaptive_k_policy(self.path)

    def test_linear_and_shifted_serialization(self):
        obj, spec = self.load(self.linear)
        self.assertEqual(obj, self.linear)
        self.assertEqual(spec, 'linear,4.0,3.0')
        shifted = dict(schema='dna-adaptive-k-v1', family='shifted_exponential',
                       base=5, scale=20, rate=0.12345678901234566, cap=40,
                       feature='global_UniAnn_CDS_exon_count', rounding='ceil',
                       status='experimental', source='study')
        _, spec = self.load(shifted)
        self.assertEqual(spec, 'shifted-exp,5.0,20.0,0.12345678901234566')

    def test_invalid_values_and_fields(self):
        for change in ({'base': float('nan')}, {'base': float('inf')}, {'base': 0},
                       {'base': True}, {'slope': -1}, {'slope': float('-inf')},
                       {'cap': 1.0}, {'cap': True}, {'cap': 0}, {'cap': 32768},
                       {'feature': 'length'}, {'rounding': 'floor'}, {'slpoe': 1},
                       {'schema': 'other'}, {'family': 'power'}, {'family': []}, {'scale': 1}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.load({**self.linear, **change})
        for field in ('scale', 'rate'):
            for value in (-1, float('nan'), float('inf')):
                policy = dict(schema='dna-adaptive-k-v1', family='shifted_exponential',
                              base=5, scale=1, rate=1, cap=40)
                policy[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    self.load(policy)

    def args(self, extra=()):
        return cli.parser().parse_args(['decode', '--input', str(self.root),
                                        '--output', str(self.root/'result.gff'), *extra])

    def test_default_and_conflict(self):
        args = self.args()
        self.assertIsNone(args.k)
        with patch.object(cli, 'job_from_args', side_effect=RuntimeError('stop')):
            with self.assertRaisesRegex(RuntimeError, 'stop'):
                cli.decode(args)
        self.assertEqual(args.k, 10)
        self.load(self.linear)
        args = self.args(['--adaptive-k-policy', str(self.path), '--k', '39'])
        with self.assertRaisesRegex(ValueError, 'must equal'):
            cli.decode(args)

    def run_mock_decode(self, extra):
        args = self.args(extra)
        fasta = self.root/'seq.fa'; fasta.write_text('>x\nAAA\n')
        job = dict(fasta=fasta, emissions=self.root/'ps', gt=self.root/'gt', ag=self.root/'ag',
                   atg=self.root/'atg', stop=self.root/'stop', psauron=None, scores=None,
                   best_gff=None, log=None, strand='plus')
        commands = []
        def process(command, log, cwd):
            commands.append(command)
            (cwd/'paths.gff').write_text('')
            (cwd/'paths.tsv').write_text('rank\n')
            log.write_text('')
        def convert(raw_path, segment, strand, output, stats, **kwargs):
            output.write_text('##gff-version 3\n')
            return {}
        coords = {'chromosome':'x','offset':0,'length':3,'total':3,'segmented':False}
        with patch.object(cli,'job_from_args',return_value=job), patch.object(cli,'find_core',return_value=Path('/core')), \
             patch.object(cli,'coordinates',return_value=(coords,'x')), \
             patch.object(cli,'run_process',side_effect=process), patch.object(cli,'convert',side_effect=convert):
            cli.decode(args)
        metadata = json.loads(Path(str(args.output)+'.json').read_text())
        return commands[0], metadata

    def test_native_command_and_metadata(self):
        self.load(self.linear)
        cmd, meta = self.run_mock_decode(['--adaptive-k-policy', str(self.path), '--k', '40'])
        self.assertEqual(cmd[cmd.index('--k')+1], 40)
        self.assertEqual(cmd[cmd.index('--adaptive-k')+1], 'linear,4.0,3.0')
        self.assertEqual(meta['adaptive_k_policy'], {'source_path':str(self.path), 'policy':self.linear})

    def test_ordinary_command_unchanged(self):
        cmd, meta = self.run_mock_decode([])
        self.assertEqual(cmd[cmd.index('--k')+1], 10)
        self.assertNotIn('--adaptive-k', cmd)
        self.assertNotIn('adaptive_k_policy', meta)

    @unittest.skipUnless(CORE.is_file(), 'Build altann-core first')
    def test_actual_adaptive_raw_paths_equal_fixed_prefix_per_locus(self):
        inputs = write_fixture(self.root/'inputs')
        # Gene one has no usable donor, so its baseline is single-exon.
        # Gene two has >8 alternative acceptors and a two-exon baseline.
        inputs['gt'].write_text('780\t2\n')
        sequence = list(inputs['fasta'].read_text().splitlines()[1])
        positions = list(range(838,899,3))
        for p in positions:
            sequence[p:p+2] = 'AG'
        inputs['ag'].write_text(''.join(f'{p}\t2\n' for p in positions))
        inputs['fasta'].write_text('>fixture\n'+''.join(sequence)+'\n')
        self.load(dict(schema='dna-adaptive-k-v1',family='shifted_exponential',
                       base=5,scale=20,rate=.1,cap=40))
        outputs = {}
        for name,extra in [('fixed',['--k','40']),('adaptive',['--adaptive-k-policy',str(self.path)])]:
            out = self.root/(name+'.gff3')
            cmd = [sys.executable,'-m','altann','decode','--core',str(CORE),
                   '--metadata-states','--raw','--rerun-viterbi','--flank','50','--threads','2',
                   '--seqid','fixture','--offset','0','--sequence-length','1200','--output',str(out)]
            for key,flag in [('fasta','fasta'),('emissions','emissions'),('gt','gt'),('ag','ag'),('atg','atg'),('stop','stop')]:
                cmd += ['--'+flag,str(inputs[key])]
            subprocess.run(cmd+extra,check=True,capture_output=True,text=True)
            with Path(str(out)+'.tsv').open() as h:
                outputs[name] = list(csv.DictReader(h,delimiter='\t'))
        wanted = {'1':5,'2':8}
        self.assertEqual({r['reference_gene_index'] for r in outputs['adaptive']},set(wanted))
        self.assertTrue(all(int(r['search_k'])==wanted[r['reference_gene_index']] for r in outputs['adaptive']))
        expected = [r for r in outputs['fixed'] if int(r['rank']) <= wanted[r['reference_gene_index']]]
        actual = [{k:v for k,v in r.items() if k!='search_k'} for r in outputs['adaptive']]
        self.assertEqual(actual,expected)
        self.assertTrue(any(int(r['rank'])>8 for r in outputs['fixed']))


if __name__ == '__main__':
    unittest.main()
