#!/usr/bin/env python3
"""Regenerate independent coding-frame fixtures with upstream UniAnn.

Run from AltAnn: python3 tests/generate_coding_fixture.py --uniann ../UniAnn
Requires g++, bash, Perl, and UniAnn/gffread (or --gffread). AltAnn is never
imported or executed here.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def generate(upstream, destination, gffread=None):
    upstream = upstream.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    length = 1740
    sequence = list("C" * length)
    genes = []
    sites = []
    for start, intron_length in zip((60, 661, 1262), (60, 61, 62)):
        donor = start + 121
        acceptor = donor + intron_length - 2
        stop = start + intron_length + 330
        for position, motif in ((start, "ATG"), (donor, "GT"),
                                (acceptor, "AG"), (acceptor + 3, "AG"),
                                (stop, "TAA")):
            sequence[position:position + len(motif)] = motif
        genes.append({"start0": start, "donor0": donor,
                      "intron_length": intron_length, "stop0": stop,
                      "cds_1based": [[start + 1, donor, 0],
                                     [donor + intron_length + 1, stop + 3, 2]]})
        for position, kind, probability in ((start, "start", .99),
                                             (donor, "donor", .99),
                                             (acceptor, "acceptor", .99),
                                             (acceptor + 3, "acceptor", .985),
                                             (stop, "stop", .99)):
            sites.append(f"coding_fixture\t{position + 1}\t+\t{kind}\t.\t{probability}\n")
    fasta = destination / "sequence.fna"
    fasta.write_text(">coding_fixture\n" + "".join(sequence) + "\n")
    (destination / "scores.txt").write_text("".join(sites))
    raw = []
    for position in range(length):
        probabilities = [1e-8] * 4
        selected = 3
        for gene in genes:
            if gene["start0"] <= position < gene["donor0"]:
                frame = gene["start0"] % 3
                selected = (position - frame - 1) % 3
            elif gene["donor0"] + gene["intron_length"] <= position <= gene["stop0"] + 2:
                frame = (gene["start0"] + gene["intron_length"]) % 3
                selected = (position - frame - 1) % 3
        probabilities[selected] = 1 - 3e-8
        raw.append(f"coding_fixture\t{position + 1}\t+\t" +
                   "\t".join(f"{value:.8f}" for value in probabilities) + "\n")
    (destination / "coding.txt").write_text("".join(raw))

    converter = upstream / "scripts" / "convert_coding_frame.pl"
    wrapper = upstream / "scripts" / "uniann.sh"
    source = upstream / "src" / "uniann.cpp"
    gffread = gffread or upstream / "gffread"
    if not gffread.is_file():
        raise FileNotFoundError("Provide --gffread or build UniAnn/gffread before regenerating")
    gffread = gffread.resolve()
    with tempfile.TemporaryDirectory(prefix="uniann-coding-oracle-") as work:
        stage = Path(work)
        binary = stage / "uniann"
        subprocess.run(["g++", "-std=c++17", "-O2", str(source), "-o", str(binary)], check=True)
        for script in (wrapper, converter):
            staged = stage / script.name
            shutil.copyfile(script, staged)
            staged.chmod(0o755)
        (stage / "gffread").symlink_to(gffread)
        env = dict(os.environ, PATH=str(stage) + os.pathsep + os.environ.get("PATH", ""))
        # Remove the wrapper's cache: regeneration must always use the converter.
        (destination / "out.ps.txt").unlink(missing_ok=True)
        subprocess.run(["bash", "-o", "pipefail", str(stage / "uniann.sh"), "-f", str(fasta.resolve()),
                        "-c", str((destination / "coding.txt").resolve()),
                        "-s", str((destination / "scores.txt").resolve())],
                       cwd=destination, env=env, check=True)
        Path(str(fasta) + ".uniann.gff").replace(destination / "uniann.filtered.gff3")
        inputs = [fasta, *(destination / name for name in
                  ("out.ps.txt", "out.gt.txt", "out.ag.txt", "out.atg.txt", "out.stop.txt"))]
        result = subprocess.run([str(binary), *map(str, map(Path.resolve, inputs))],
                                capture_output=True, text=True, check=True)
    (destination / "uniann.gff3").write_text(result.stdout)
    trace = [line for line in result.stderr.splitlines()
             if len(line.split("\t")) == 9 and line.split("\t")[1] in ("dp", "bt")]
    if len(trace) != length * 2:
        raise RuntimeError(f"Expected {length * 2} native trace rows, found {len(trace)}")
    (destination / "uniann.trace").write_text("\n".join(trace) + "\n")
    (destination / "out.err").unlink(missing_ok=True)
    commit = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    version = subprocess.run([str(gffread), "--version"], capture_output=True, text=True, check=True)
    provenance = {"upstream_commit": commit, "length": length, "genes": genes,
                  "upstream_sha256": {str(path.relative_to(upstream)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in (converter, wrapper, source)},
                  "gffread_version": (version.stdout + version.stderr).strip()}
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Generated {length} positions from UniAnn {commit} in {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniann", type=Path, required=True)
    parser.add_argument("--gffread", type=Path)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parent / "data" / "coding-frame")
    args = parser.parse_args()
    generate(args.uniann, args.output.resolve(), args.gffread)
