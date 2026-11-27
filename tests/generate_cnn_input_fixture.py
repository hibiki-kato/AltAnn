#!/usr/bin/env python3
"""Regenerate the CNN-Mamba input fixture with upstream UniAnn f126cc8.

Run from AltAnn: python3 tests/generate_cnn_input_fixture.py --uniann ../UniAnn
f126cc8 is the last upstream converter without a temperature transform
(08eaf15 applies sigmoid(logit(p)/3) to each class). The upstream wrapper,
converter and decoder are exported with git archive, the decoder is compiled,
and uniann.sh runs on ChimAnn's staged form of coding.tsv (header removed,
frame columns rotated left by one for the (n-i-1)%3 formula). AltAnn is never
imported or executed here. Requires git, g++, Bash, Perl and UniAnn/gffread.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tarfile
import tempfile

HEADER = ("chrom\tposition_1based\tstrand\tframe0_probability\tframe1_probability\t"
          "frame2_probability\tnoncoding_probability\tpredicted_state\n")
SITE_HEADER = "chrom\tpos\tstrand\ttype\tmotif\tprob\n"
UPSTREAM_FILES = ("scripts/uniann.sh", "scripts/convert_coding_frame.pl", "src/uniann.cpp")
OUTPUTS = ("out.ps.txt", "out.gt.txt", "out.ag.txt", "out.atg.txt", "out.stop.txt")


def inputs(destination, seqid="cnn_fixture", length=1500):
    rng = random.Random(20261127)
    bases = [rng.choice("ACGT") for _ in range(length)]
    for position in range(400, 430):
        bases[position] = bases[position].lower()  # soft-masked bases are upper-cased
    bases[700] = "N"
    sequence = "".join(bases)
    (destination / "sequence.fna").write_text(f">{seqid} oriented fixture\n" + sequence[:800] + "\n"
                                              + sequence[800:] + "\n")
    edge = {10: (0, 0, 0, 1), 11: (1, 0, 0, 0), 12: (0, 1, 0, 0), 13: (0, 0, 1, 0),
            14: (1e-12, 1e-12, 1e-12, 1 - 3e-12), 15: (0.25, 0.25, 0.25, 0.25)}
    rows = [HEADER]
    for n in range(length):
        if n in edge:
            values = edge[n]
        else:
            weights = [rng.random() ** 3 for _ in range(4)]
            total = sum(weights)
            values = [weight / total for weight in weights]
        state = ("frame0", "frame1", "frame2", "noncoding")[max(range(4), key=values.__getitem__)]
        rows.append(f"{seqid}\t{n + 1}\t+\t" + "\t".join(f"{value:.8g}" for value in values) + f"\t{state}\n")
    (destination / "coding.tsv").write_text("".join(rows))
    sites = [SITE_HEADER]
    for kind, motif in (("donor", "GT"), ("acceptor", "AG"), ("start", "ATG"), ("stop", "TAA")):
        for position in sorted(rng.sample(range(1, length + 1), 25)):
            sites.append(f"{seqid}\t{position}\t+\t{kind}\t{motif}\t{rng.random():.8g}\n")
    sites += [f"{seqid}\t100\t+\tdonor\tGT\t0.995\n",  # guarantees a positive donor maximum
              f"{seqid}\t101\t+\tacceptor\tAG\t0\n",
              f"{seqid}\t102\t+\tstart\tATG\t1\n",
              f"{seqid}\t103\t+\tdonor\tGT\t0.2\t0.75\n",  # seventh column replaces the sixth
              f"{seqid}\t104\t-\tdonor\tGT\t0.999\n"]  # other strand rows are ignored
    # Upstream also scores + rows naming another sequence; AltAnn skips them (tested separately).
    (destination / "sites.tsv").write_text("".join(sites))


def staged_coding(source, destination):
    with source.open() as handle, destination.open("w") as output:
        for number, line in enumerate(handle, 1):
            if number == 1 and line == HEADER:
                continue
            row = line.rstrip("\n").split("\t")
            row[3:6] = [row[4], row[5], row[3]]
            output.write("\t".join(row) + "\n")


def generate(upstream, commit, destination, gffread):
    destination.mkdir(parents=True, exist_ok=True)
    inputs(destination)
    with tempfile.TemporaryDirectory(prefix="uniann-cnn-oracle-") as work:
        work = Path(work)
        archive = work / "upstream.tar"
        subprocess.run(["git", "-C", str(upstream), "archive", "--output", str(archive), commit, *UPSTREAM_FILES],
                       check=True)
        with tarfile.open(archive) as handle:
            handle.extractall(work / "src")
        stage = work / "stage"
        stage.mkdir()
        subprocess.run(["g++", "-std=c++17", "-O2", str(work / "src/src/uniann.cpp"), "-o", str(stage / "uniann")],
                       check=True)
        for name in ("uniann.sh", "convert_coding_frame.pl"):
            shutil.copyfile(work / "src/scripts" / name, stage / name)
            (stage / name).chmod(0o755)
        (stage / "gffread").symlink_to(gffread)
        run = work / "run"
        run.mkdir()
        shutil.copyfile(destination / "sequence.fna", run / "sequence.fna")
        shutil.copyfile(destination / "sites.tsv", run / "sites.tsv")
        staged_coding(destination / "coding.tsv", run / "coding.tsv")
        env = dict(os.environ, PATH=str(stage) + os.pathsep + os.environ.get("PATH", ""))
        subprocess.run(["bash", "-o", "pipefail", str(stage / "uniann.sh"), "-f", "sequence.fna",
                        "-c", "coding.tsv", "-s", "sites.tsv", "-n"], cwd=run, env=env, check=True)
        for name in OUTPUTS:
            shutil.copyfile(run / name, destination / name)
        hashes = {name: hashlib.sha256((work / "src" / name).read_bytes()).hexdigest() for name in UPSTREAM_FILES}
    full = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", commit], text=True).strip()
    provenance = {"upstream_commit": full, "upstream_sha256": hashes,
                  "staging": "header removed; frame columns rotated (1, 2, 0) as ChimAnn does for (n-i-1)%3",
                  "wrapper_arguments": ["-f", "sequence.fna", "-c", "coding.tsv", "-s", "sites.tsv", "-n"],
                  "expected_outputs": {name: hashlib.sha256((destination / name).read_bytes()).hexdigest()
                                       for name in OUTPUTS}}
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Generated {destination} from UniAnn {full}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniann", type=Path, required=True)
    parser.add_argument("--commit", default="f126cc8")
    parser.add_argument("--gffread", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "data" / "cnn-input")
    args = parser.parse_args()
    gffread = (args.gffread or args.uniann / "gffread").resolve()
    if not gffread.is_file():
        raise SystemExit("Provide --gffread or a runnable UniAnn/gffread")
    generate(args.uniann.resolve(), args.commit, args.output.resolve(), gffread)
