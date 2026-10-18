#!/usr/bin/env python3
"""Evaluate coding predictions with the project's frozen S. pombe conventions.

Gene matches use transcript OR strict CDS-as-exon trmap matches. Exact CDS
counts use chromosome, strand, and all ordered CDS intervals, independent of
model IDs. Candidate pools are explicitly labeled; no candidates are selected.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


@dataclass
class Model:
    chromosome: str
    strand: str
    identifier: str
    attributes: dict = field(default_factory=dict)
    exons: list = field(default_factory=list)
    cds: list = field(default_factory=list)

    @property
    def structure(self):
        return self.chromosome, self.strand, tuple(sorted(set(self.cds)))

    @property
    def gene(self):
        attrs = self.attributes
        return (
            self.chromosome,
            self.strand,
            attrs.get("kbest_locus") or attrs.get("geneID")
            or attrs.get("Parent") or self.identifier,
        )


def read_models(path, chromosomes=()):
    models = {}
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected nine GFF3 fields")
            chrom, _, kind, start, end, _, strand, _, text = fields
            if chromosomes and chrom not in chromosomes:
                continue
            if kind not in {"mRNA", "transcript", "exon", "CDS"}:
                continue
            attrs = dict(item.split("=", 1) for item in text.split(";") if "=" in item)
            ids = [attrs.get("ID")] if kind in {"mRNA", "transcript"} else attrs.get("Parent", "").split(",")
            for identifier in filter(None, ids):
                # Independent chromosome/strand runs can reuse transcript IDs.
                key = chrom, strand, identifier
                model = models.setdefault(key, Model(chrom, strand, identifier))
                if kind in {"mRNA", "transcript"}:
                    model.attributes.update(attrs)
                else:
                    interval = int(start), int(end)
                    if interval[0] > interval[1] or interval[0] < 1:
                        raise ValueError(f"{path}:{number}: invalid GFF3 interval")
                    (model.cds if kind == "CDS" else model.exons).append(interval)
    return {key: value for key, value in models.items() if value.cds or value.exons}


def write_models(models, path, cds_only=False):
    """Serialize safe unique IDs for trmap, retaining all CDS-chain boundaries."""
    lookup = {}
    with Path(path).open("w") as handle:
        handle.write("##gff-version 3\n")
        for index, model in enumerate(models.values()):
            if cds_only and not model.cds:
                continue
            identifier = f"model_{index}"
            lookup[identifier] = model
            coordinates = sorted(set(model.cds if cds_only else model.exons or model.cds))
            handle.write(
                f"{model.chromosome}\tevaluation\ttranscript\t{coordinates[0][0]}\t"
                f"{coordinates[-1][1]}\t.\t{model.strand}\t.\tID={identifier}\n"
            )
            for start, end in coordinates:
                handle.write(
                    f"{model.chromosome}\tevaluation\texon\t{start}\t{end}\t.\t"
                    f"{model.strand}\t.\tParent={identifier}\n"
                )
    return lookup


def parse_trmap(text):
    pairs = set()
    query = None
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0].startswith(">"):
            query = fields[0][1:]
        elif fields[0] == "=":
            if query is None or len(fields) < 6:
                raise ValueError("Malformed trmap match record")
            pairs.add((query, fields[5]))
    return pairs


def run_trmap(trmap, reference, query, strict=False):
    command = [str(trmap), "-c", "="]
    if strict:
        command.append("--strict-match")
    result = subprocess.run(command + [str(reference), str(query)], capture_output=True, text=True, check=True)
    return parse_trmap(result.stdout)


def percentage(numerator, denominator):
    return 100.0 * numerator / denominator if denominator else None


def evaluate(reference, query, pairs, candidate=False):
    reference_genes = {model.gene for model in reference.values()}
    query_genes = {model.gene for model in query.values()}
    matched_reference_genes = {reference[ref_id].gene for _, ref_id in pairs}
    matched_query_genes = {query[query_id].gene for query_id, _ in pairs}
    reference_cds = {model.structure for model in reference.values() if model.cds}
    query_cds = {model.structure for model in query.values() if model.cds}
    exact = reference_cds & query_cds
    origins = Counter(model.attributes.get("kbest_origin", "unspecified") for model in query.values())
    return {
        "interpretation": "candidate_pool" if candidate else "baseline_annotation",
        "precision_scope": "candidate_pool" if candidate else "baseline_annotation",
        "reference_genes": len(reference_genes),
        "query_gene_loci": len(query_genes),
        "matched_reference_genes": len(matched_reference_genes),
        "matched_query_gene_loci": len(matched_query_genes),
        "gene_sensitivity_percent": percentage(len(matched_reference_genes), len(reference_genes)),
        "gene_precision_percent": percentage(len(matched_query_genes), len(query_genes)),
        "reference_coding_transcripts": sum(bool(model.cds) for model in reference.values()),
        "query_coding_transcripts": sum(bool(model.cds) for model in query.values()),
        "query_transcripts": len(query),
        "reference_unique_cds": len(reference_cds),
        "query_unique_cds": len(query_cds),
        "exact_match_unique_cds": len(exact),
        "exact_match_reference_transcripts": sum(model.structure in exact for model in reference.values()),
        "exact_match_query_transcripts": sum(model.structure in exact for model in query.values()),
        "cds_sensitivity_percent": percentage(len(exact), len(reference_cds)),
        "cds_precision_percent": percentage(len(exact), len(query_cds)),
        "alternative_isoform_transcripts": sum(
            model.attributes.get("kbest_origin") == "local_kbest"
            for model in query.values()
        ),
        "kbest_origins": dict(origins),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--predictions", action="append", required=True, metavar="LABEL=PATH")
    parser.add_argument("--output", required=True, type=Path, help="Output JSON report")
    parser.add_argument("--chromosome", action="append", default=[])
    parser.add_argument("--candidate-label", action="append", default=[], help="Mark a prediction set as an unselected candidate pool")
    parser.add_argument("--trmap", default="trmap")
    args = parser.parse_args()
    predictions = {}
    for value in args.predictions:
        label, separator, path = value.partition("=")
        if not separator or not label or label in predictions:
            parser.error("--predictions requires a unique LABEL=PATH")
        predictions[label] = Path(path)
    if set(args.candidate_label) - predictions.keys():
        parser.error("--candidate-label must name a --predictions label")
    models = read_models(args.reference, args.chromosome)
    if not any(model.cds for model in models.values()):
        parser.error("Reference has no CDS models on the requested chromosomes")
    report = {
        "reference": str(args.reference.resolve()),
        "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        "chromosomes": sorted(args.chromosome or {model.chromosome for model in models.values()}),
        "definitions": {
            "gene_match": "transcript trmap -c '=' OR CDS-as-exons trmap -c '=' --strict-match",
            "reference_gene": "geneID or Parent, scoped by chromosome and strand",
            "query_gene_locus": "kbest_locus, then geneID/Parent, then transcript ID; scoped by chromosome and strand",
            "exact_cds": "unique (chromosome, strand, sorted complete CDS intervals), no coordinate tolerance; no partial/pseudo filtering beyond the supplied frozen reference",
            "sensitivity": "matched reference genes / reference genes; exact unique CDS matches / reference unique CDS",
            "precision": "matched query gene loci / query loci; exact unique CDS matches / query unique CDS; candidate pool precision is not final annotation precision",
        },
        "predictions": {},
    }
    structures = {}
    with tempfile.TemporaryDirectory(prefix="spom_evaluation_") as directory:
        work = Path(directory)
        reference_lookup = write_models(models, work / "reference.gff3")
        write_models(models, work / "reference.cds.gff3", cds_only=True)
        for label, path in predictions.items():
            query = read_models(path, args.chromosome)
            lookup = write_models(query, work / "query.gff3")
            write_models(query, work / "query.cds.gff3", cds_only=True)
            pairs = run_trmap(args.trmap, work / "reference.gff3", work / "query.gff3")
            pairs |= run_trmap(args.trmap, work / "reference.cds.gff3", work / "query.cds.gff3", strict=True)
            result = evaluate(reference_lookup, lookup, pairs, label in args.candidate_label)
            result.update(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            report["predictions"][label] = result
            structures[label] = {model.structure for model in query.values() if model.cds}
    reference_structures = {model.structure for model in models.values() if model.cds}
    report["cds_set_comparisons"] = {
        label: {
            other: {"shared_unique_cds": len(keys & values),
                    "additional_unique_cds": len(keys - values),
                    "additional_exact_reference_cds": len((keys - values) & reference_structures),
                    "missing_unique_cds": len(values - keys)}
            for other, values in structures.items() if other != label
        }
        for label, keys in structures.items()
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["predictions"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
