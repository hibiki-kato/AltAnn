#!/usr/bin/env python3
"""Replay frozen synthetic-logit Bowen export through working-tree UniAnn.

Uses only the Python standard library, g++, Bash, Perl and gffread. No neural
model inference, real S. pombe predictions or AltAnn decoder is involved.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


ORIGINAL = Path(__file__).resolve().parent / "data/bowen-contract"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate(uniann, destination, gffread=None):
    uniann, destination = uniann.resolve(), destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    original = json.loads((ORIGINAL / "provenance.json").read_text())
    input_names = ("sequence.fna", "scores.txt", "exported.tsv")
    for name in input_names:
        shutil.copyfile(ORIGINAL / name, destination / name)
    converter = uniann / "scripts/convert_coding_frame.pl"
    wrapper = uniann / "scripts/uniann.sh"
    source = uniann / "src/uniann.cpp"
    gffread = (gffread or uniann / "gffread").resolve()
    length = original["length"]
    with tempfile.TemporaryDirectory(prefix="bowen-corrected-") as temporary:
        stage = Path(temporary)
        binary = stage / "uniann"
        subprocess.run(["g++", "-std=c++17", "-O2", str(source), "-o", str(binary)], check=True)
        for script in (wrapper, converter):
            copied = stage / script.name
            shutil.copyfile(script, copied)
            copied.chmod(0o755)
        (stage / "gffread").symlink_to(gffread)
        env = dict(os.environ, PATH=str(stage) + os.pathsep + os.environ.get("PATH", ""))
        # Exercise regeneration of stale cached scores as well as header parsing.
        (destination / "out.ps.txt").write_text("stale emission cache\n")
        subprocess.run(["bash", str(stage / "uniann.sh"), "-f", str(destination / "sequence.fna"),
                        "-c", str(destination / "exported.tsv"), "-s", str(destination / "scores.txt")],
                       cwd=destination, env=env, check=True)
        Path(str(destination / "sequence.fna") + ".uniann.gff").replace(destination / "uniann.filtered.gff3")
        inputs = [destination / name for name in
                  ("sequence.fna", "out.ps.txt", "out.gt.txt", "out.ag.txt", "out.atg.txt", "out.stop.txt")]
        native = subprocess.run([str(binary), *map(str, inputs)], text=True, capture_output=True, check=True)
        (destination / "uniann.gff3").write_text(native.stdout)
        trace = [line for line in native.stderr.splitlines()
                 if len(line.split("\t")) == 9 and line.split("\t")[1] in ("dp", "bt")]
        if len(trace) != 2 * length:
            raise RuntimeError("Incomplete native trace")
        (destination / "uniann.trace").write_text("\n".join(trace) + "\n")
    (destination / "out.err").unlink(missing_ok=True)
    if len((destination / "out.ps.txt").read_text().splitlines()) != length:
        raise RuntimeError("Headered exporter input must produce exactly one emission per base")
    provenance = {key: original[key] for key in
                  ("synthetic_logits", "checkpoint_inference", "length", "genes", "temperature",
                   "internal_class_order", "phase_labels", "minus_reflection_verified", "bowen_commit")}
    provenance.update({
        "description": "Frozen synthetic-logit real-exporter TSV replayed through corrected working-tree UniAnn; not real S. pombe model predictions",
        "uniann_base_commit": original["uniann_commit"],
        "uniann_working_tree_modified": True,
        "input_sha256": {name: sha256(ORIGINAL / name) for name in input_names},
        "original_provenance_sha256": sha256(ORIGINAL / "provenance.json"),
        "original_source_sha256": original["source_sha256"],
        "source_sha256": {"uniann/" + str(path.relative_to(uniann)): sha256(path)
                          for path in (converter, wrapper, source)},
        "gffread_sha256": sha256(gffread),
        "emission_rows": length,
        "header_included_in_input": True,
        "additional_temperature_scaling": False,
        "absolute_frame_formula": "(zero_based_position - local_phase) % 3",
    })
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Generated corrected synthetic-logit contract for {length} bases in {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniann", type=Path, required=True)
    parser.add_argument("--gffread", type=Path)
    parser.add_argument("--output", type=Path, default=ORIGINAL.parent / "bowen-corrected")
    args = parser.parse_args()
    generate(args.uniann, args.output, args.gffread)
