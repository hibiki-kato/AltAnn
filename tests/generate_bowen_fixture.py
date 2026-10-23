#!/usr/bin/env python3
"""Generate a pinned Bowen -> UniAnn contract diagnostic, without AltAnn.

Requires numpy, pandas, torch, g++, Perl, and UniAnn/gffread. Neural network
weights and GPU inference are deliberately not used: logits are synthetic.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def generate(bowen, uniann, destination):
    bowen, uniann = bowen.resolve(), uniann.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    original = Path(__file__).resolve().parent / "data" / "coding-frame"
    geometry = json.loads((original / "provenance.json").read_text())
    for name in ("sequence.fna", "scores.txt"):
        shutil.copyfile(original / name, destination / name)
    length, genes = geometry["length"] - 40, geometry["genes"]
    fasta = destination / "sequence.fna"
    fasta_lines = fasta.read_text().splitlines()
    fasta.write_text(fasta_lines[0] + "\n" + "".join(fasta_lines[1:])[:length] + "\n")
    records = [(start - 1, end, phase) for gene in genes
               for start, end, phase in gene["cds_1based"]]
    prepare = bowen / "src/cnn_mamba/prepare_droso.py"
    # Execute the actual function unchanged, without importing its pysam-dependent
    # data preparation entry point. The extracted AST needs only NumPy.
    function = next(node for node in ast.parse(prepare.read_text()).body
                    if isinstance(node, ast.FunctionDef) and node.name == "build_phase_track")
    scope = {"np": np}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(prepare), "exec"), scope)
    labels = scope["build_phase_track"](records, length, "+")
    minus_records = [(length - end, length - start, phase) for start, end, phase in records]
    minus = scope["build_phase_track"](minus_records, length, "-")
    if not np.array_equal(labels, minus[::-1]):
        raise RuntimeError("Bowen plus/minus phase-track reflection differs")
    targets = np.where(labels >= 0, labels + 1, 0)
    profiles = np.full((4, 4), 1e-8, dtype=np.float64)
    np.fill_diagonal(profiles, 1 - 3e-8)
    logits = np.log(profiles)[targets]
    calibration = bowen / "results/human_joint_frame_beta010_v1/temperature_scaling/temperature.json"
    temperature = json.loads(calibration.read_text())["temperature"]
    scaling = bowen / "src/cnn_mamba/temperature_scaling.py"
    spec = importlib.util.spec_from_file_location("bowen_temperature", scaling)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    probabilities = module.apply_temperature_numpy(logits, temperature)
    exporter = bowen / "src/cnn_mamba/export_joint_frame_tsv.py"
    converter = uniann / "scripts/convert_coding_frame.pl"
    wrapper = uniann / "scripts/uniann.sh"
    source = uniann / "src/uniann.cpp"
    with tempfile.TemporaryDirectory(prefix="bowen-contract-") as temporary:
        stage = Path(temporary)
        np.save(stage / "probabilities.npy", probabilities)
        np.save(stage / "argmax.npy", probabilities.argmax(axis=1).astype(np.uint8))
        exported = destination / "exported.tsv"
        exported.unlink(missing_ok=True)
        Path(str(exported) + ".complete.json").unlink(missing_ok=True)
        subprocess.run([sys.executable, str(exporter), "--probabilities", str(stage / "probabilities.npy"),
                        "--argmax", str(stage / "argmax.npy"), "--output", str(exported),
                        "--chrom", "coding_fixture", "--chunk-size", "317"], check=True)
        # Only remove the header; retain the exporter's extra predicted_state column.
        lines = exported.read_text().splitlines(keepends=True)
        (destination / "coding.txt").write_text("".join(lines[1:]))
        sentinel = json.loads(Path(str(exported) + ".complete.json").read_text())
        Path(str(exported) + ".complete.json").unlink()
        header_result = subprocess.run(["perl", str(converter), str(destination / "sequence.fna")],
                                       input="".join(lines), text=True, capture_output=True, check=True)
        (destination / "header-diagnostic.txt").write_text(
            "\n".join(header_result.stdout.splitlines()[:3]) + "\n")
        binary = stage / "uniann"
        subprocess.run(["g++", "-std=c++17", "-O2", str(source), "-o", str(binary)], check=True)
        for script in (wrapper, converter):
            copied = stage / script.name
            shutil.copyfile(script, copied)
            copied.chmod(0o755)
        (stage / "gffread").symlink_to(uniann / "gffread")
        env = dict(os.environ, PATH=str(stage) + os.pathsep + os.environ.get("PATH", ""))
        (destination / "out.ps.txt").unlink(missing_ok=True)
        subprocess.run(["bash", "-o", "pipefail", str(stage / "uniann.sh"),
                        "-f", str(destination / "sequence.fna"), "-c", str(destination / "coding.txt"),
                        "-s", str(destination / "scores.txt")], cwd=destination, env=env, check=True)
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
    provenance = {
        "synthetic_logits": True, "checkpoint_inference": False,
        "length": length, "genes": genes, "temperature": temperature,
        "internal_class_order": list(module.CLASS_ORDER),
        "logit_profiles": np.log(profiles).tolist(),
        "calibrated_profiles": module.apply_temperature_numpy(np.log(profiles), temperature).tolist(),
        "phase_labels": "".join(str(int(value)) if value >= 0 else "." for value in labels),
        "minus_reflection_verified": True,
        "export_header": lines[0].rstrip("\n").split("\t"),
        "export_probability_sum_error": sentinel["maximum_probability_sum_error"],
        "unstripped_converter_rows": len(header_result.stdout.splitlines()),
        "bowen_commit": subprocess.check_output(["git", "-C", str(bowen), "rev-parse", "HEAD"], text=True).strip(),
        "uniann_commit": subprocess.check_output(["git", "-C", str(uniann), "rev-parse", "HEAD"], text=True).strip(),
        "source_sha256": {"bowen/" + str(path.relative_to(bowen)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (prepare, scaling, exporter, calibration)},
    }
    provenance["source_sha256"].update({"uniann/" + str(path.relative_to(uniann)): hashlib.sha256(path.read_bytes()).hexdigest()
                                        for path in (converter, wrapper, source)})
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Generated synthetic-logit Bowen contract for {length} positions in {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bowen", type=Path, required=True)
    parser.add_argument("--uniann", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "data/bowen-contract")
    args = parser.parse_args()
    generate(args.bowen, args.uniann, args.output.resolve())
