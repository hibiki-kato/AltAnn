"""Smoke-test the installed launcher without relying on checkout module paths."""
from pathlib import Path
import subprocess
import sys
import tempfile

from fixtures import ordered_inputs, write_fixture, write_raw_fixture

launcher = Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix="altann-installed-") as work:
    work = Path(work)
    inputs = write_fixture(work / "inputs")
    command = [str(launcher), "decode", "--rerun-viterbi", "--k", "2",
               "--threads", "2", "--flank", "50", "--output", str(work / "out.gff3")]
    for name, value in zip(("fasta", "emissions", "gt", "ag", "atg", "stop"), ordered_inputs(inputs)):
        command.extend(["--" + name, value])
    subprocess.run(command, check=True, cwd=work)
    assert "kbest_origin=uniann_best" in (work / "out.gff3").read_text()
    combined, _, _ = write_raw_fixture(work / 'stranded')
    output = work / 'both.gff3'
    subprocess.run([str(launcher), 'decode', '--input', str(combined['fasta'].parent),
                    '-a', '--k', '2', '--threads', '2', '--flank', '50',
                    '--output', str(output)], check=True, cwd=work)
    strands = {line.split('\t')[6] for line in output.read_text().splitlines()
               if line and not line.startswith('#')}
    assert strands == {'+', '-'}
