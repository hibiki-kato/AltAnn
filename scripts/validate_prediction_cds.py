#!/usr/bin/env python3
"""Read-only GFF3 structure, spliced CDS biology, and serialized score audit.

Uses one-based inclusive GFF coordinates on the original genome. CDS pieces
are joined in transcript direction WITHOUT removing phase bases at each exon:
those bases complete the codon begun in the preceding exon. This is a validity
audit, not an annotation precision/recall evaluation or a candidate filter.
"""

import argparse
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote


COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
STOPS = {"TAA", "TAG", "TGA"}


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.resolve()), "sha256": digest.hexdigest()}


def load_genome(path):
    sequences, name, parts = {}, None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
                if name in sequences:
                    raise ValueError("Duplicate FASTA identifier: " + name)
            else:
                if name is None and line.strip():
                    raise ValueError("FASTA sequence before header")
                parts.append("".join(line.split()))
    if name is not None:
        sequences[name] = "".join(parts).upper()
    if not sequences:
        raise ValueError("Empty genome FASTA")
    return sequences


def origin(record):
    attrs = record["attrs"]
    if attrs.get("kbest_origin") == "local_kbest":
        return "candidates"
    if attrs.get("kbest_origin") == "uniann_best":
        return "baseline"
    return "candidates" if attrs.get("kbest_rank", "0") != "0" else "baseline"


def validate(genome, gff):
    genome, gff = Path(genome), Path(gff)
    sequences = load_genome(genome)
    records, ids, structural = [], {}, []

    def row_error(record, error):
        record["errors"].add(error)
        structural.append({"line": record["line"], "feature": record["kind"],
                           "id": record["attrs"].get("ID"), "error": error})

    with gff.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) != 9:
                structural.append({"line": line_number, "error": "expected_9_columns"})
                continue
            attrs, attr_errors = {}, []
            for field in row[8].split(";"):
                if field in ("", "."):
                    continue
                if "=" not in field:
                    attr_errors.append("malformed_attribute")
                    continue
                key, value = field.split("=", 1)
                if key in attrs:
                    attr_errors.append("duplicate_attribute")
                attrs[key] = unquote(value) if key != "Parent" else value
            record = {"line": line_number, "chrom": unquote(row[0]), "kind": row[2],
                      "strand": row[6], "phase": row[7], "attrs": attrs,
                      "parents": [unquote(parent) for parent in attrs.get("Parent", "").split(",") if parent],
                      "errors": set(), "start": None, "end": None}
            records.append(record)
            for error in attr_errors:
                row_error(record, error)
            if row[5] != ".":
                try:
                    if not Decimal(row[5]).is_finite():
                        row_error(record, "invalid_GFF_score")
                except InvalidOperation:
                    row_error(record, "invalid_GFF_score")
            try:
                record["start"], record["end"] = int(row[3]), int(row[4])
            except ValueError:
                row_error(record, "noninteger_coordinates")
            if record["chrom"] not in sequences:
                row_error(record, "unknown_genome_sequence")
            elif record["start"] is not None and not (1 <= record["start"] <= record["end"] <= len(sequences[record["chrom"]])):
                row_error(record, "coordinates_outside_original_genome")
            if record["kind"] in ("gene", "mRNA", "transcript", "exon", "CDS") and record["strand"] not in ("+", "-"):
                row_error(record, "invalid_strand")
            if record["kind"] == "CDS" and record["phase"] not in ("0", "1", "2"):
                row_error(record, "invalid_CDS_phase")
            elif record["kind"] != "CDS" and record["phase"] != ".":
                row_error(record, "phase_on_non_CDS_feature")
            if record["kind"] in ("CDS", "exon") and not record["parents"]:
                row_error(record, "missing_Parent")
            identifier = attrs.get("ID")
            if not identifier and record["kind"] in ("gene", "mRNA", "transcript"):
                row_error(record, "missing_ID")
            if identifier:
                if identifier in ids:
                    previous = ids[identifier]
                    # GFF3 permits a discontinuous CDS to repeat its feature ID.
                    if not (record["kind"] == previous["kind"] == "CDS"
                            and (record["chrom"], record["strand"], record["parents"])
                            == (previous["chrom"], previous["strand"], previous["parents"])):
                        row_error(record, "duplicate_ID")
                else:
                    ids[identifier] = record

    children = defaultdict(list)
    for record in records:
        for parent in record["parents"]:
            children[parent].append(record)
            if parent not in ids:
                row_error(record, "missing_parent_feature")
                continue
            ancestor = ids[parent]
            if (record["chrom"], record["strand"]) != (ancestor["chrom"], ancestor["strand"]):
                row_error(record, "parent_chromosome_or_strand_mismatch")
            if record["start"] is not None and ancestor["start"] is not None:
                if not (ancestor["start"] <= record["start"] <= record["end"] <= ancestor["end"]):
                    row_error(record, "child_outside_parent")
            if record["kind"] in ("CDS", "exon") and ancestor["kind"] not in ("mRNA", "transcript"):
                row_error(record, "invalid_transcript_parent_type")
            if record["kind"] in ("mRNA", "transcript") and ancestor["kind"] not in ("gene", "pseudogene"):
                row_error(record, "invalid_gene_parent_type")
        identifier = record["attrs"].get("ID")
        pending, visited = list(record["parents"]), set()
        while pending:
            parent = pending.pop()
            if parent == identifier:
                row_error(record, "cyclic_Parent_relationship")
                break
            if parent not in visited and parent in ids:
                visited.add(parent)
                pending.extend(ids[parent]["parents"])

    groups, failures = defaultdict(Counter), []
    counts = Counter({"transcripts": 0, "valid_transcripts": 0, "failed_transcripts": 0})
    transcripts = [record for record in records if record["kind"] in ("mRNA", "transcript")]
    references = {}
    for record in transcripts:
        attrs = record["attrs"]
        if origin(record) == "baseline" and "kbest_locus" in attrs and "kbest_score" in attrs:
            references[attrs["kbest_locus"]] = attrs["kbest_score"]

    for record in transcripts:
        identifier = record["attrs"].get("ID", "line:" + str(record["line"]))
        attrs, strand = record["attrs"], record["strand"]
        group = groups[origin(record) + ":" + strand]
        group["transcripts"] += 1
        counts["transcripts"] += 1
        errors = set(record["errors"])
        child_records = children[identifier]
        errors.update(error for child in child_records for error in child["errors"])
        cds = [child for child in child_records if child["kind"] == "CDS"]
        exons = [child for child in child_records if child["kind"] == "exon"]
        internal_stops = []
        group["CDS_segments"] += len(cds)
        if not cds:
            errors.add("missing_CDS")
        can_extract = cds and strand in ("+", "-") and all(
            child["chrom"] == record["chrom"] and child["strand"] == strand
            and child["start"] is not None and child["chrom"] in sequences
            and 1 <= child["start"] <= child["end"] <= len(sequences[child["chrom"]]) for child in cds)
        if can_extract:
            cds.sort(key=lambda child: child["start"], reverse=strand == "-")
            pieces, cumulative, segments = [], 0, []
            intervals = sorted((child["start"], child["end"]) for child in cds)
            if any(left[1] >= right[0] for left, right in zip(intervals, intervals[1:])):
                errors.add("overlapping_CDS")
            for child in cds:
                if child["phase"] != str((-cumulative) % 3):
                    errors.add("CDS_phase_discontinuity")
                if exons and not any(exon["start"] is not None and exon["start"] <= child["start"]
                                     and child["end"] <= exon["end"] for exon in exons):
                    errors.add("CDS_outside_exon")
                piece = sequences[child["chrom"]][child["start"] - 1:child["end"]]
                pieces.append(piece if strand == "+" else piece.translate(COMPLEMENT)[::-1])
                segments.append((cumulative, cumulative + len(piece), child["start"], child["end"]))
                cumulative += len(piece)
            sequence = "".join(pieces)
            if len(sequence) % 3:
                errors.add("CDS_length_not_divisible_by_3")
            if not sequence.startswith("ATG"):
                errors.add("non_ATG_start")
            if sequence[-3:] not in STOPS:
                errors.add("missing_terminal_stop")
            for index in range(0, len(sequence) - 3, 3):
                if sequence[index:index + 3] not in STOPS:
                    continue
                positions = []
                for offset in range(index, index + 3):
                    for lo, hi, start, end in segments:
                        if lo <= offset < hi:
                            positions.append(start + offset - lo if strand == "+" else end - offset + lo)
                            break
                step = 1 if strand == "+" else -1
                kind = ("split_across_junction" if any(right - left != step for left, right in zip(positions, positions[1:]))
                        else "within_exon")
                group["internal_stop_codons:" + kind] += 1
                internal_stops.append({"codon": sequence[index:index + 3], "spliced_CDS_0based_offset": index,
                                       "genomic_positions_1based": positions, "kind": kind})
            if internal_stops:
                errors.add("internal_stop")
            if set(sequence) - set("ACGT"):
                errors.add("ambiguous_CDS_bases")
            group["CDS_bases_checked"] += len(sequence)

        score_keys = ("kbest_score", "kbest_reference_score", "kbest_delta")
        if any(key in attrs for key in score_keys):
            group["score_triples_checked"] += 1
            if not all(key in attrs for key in score_keys):
                errors.add("incomplete_score_attributes")
            else:
                try:
                    score, reference, delta = (Decimal(attrs[key]) for key in score_keys)
                    if not all(value.is_finite() for value in (score, reference, delta)):
                        errors.add("nonfinite_score")
                    else:
                        with localcontext() as context:
                            context.prec = max(50, sum(len(value.as_tuple().digits) for value in (score, reference, delta))
                                               + max(value.adjusted() for value in (score, reference, delta))
                                               - min(value.as_tuple().exponent for value in (score, reference, delta)) + 5)
                            if score - reference != delta:
                                errors.add("serialized_score_delta_mismatch")
                        if origin(record) == "baseline" and (score != reference or delta != 0):
                            errors.add("baseline_reference_score_mismatch")
                        if attrs.get("kbest_locus") in references and reference != Decimal(references[attrs["kbest_locus"]]):
                            errors.add("locus_reference_score_mismatch")
                except InvalidOperation:
                    errors.add("malformed_score")
        if errors:
            group["failed_transcripts"] += 1
            counts["failed_transcripts"] += 1
            for error in sorted(errors):
                group["error:" + error] += 1
            failures.append({"transcript": identifier, "chromosome": record["chrom"], "strand": strand,
                             "origin": origin(record), "errors": sorted(errors)})
            if internal_stops:
                failures[-1]["internal_stops"] = internal_stops
        else:
            group["valid_transcripts"] += 1
            counts["valid_transcripts"] += 1

    counts["GFF_feature_rows"] = len(records)
    counts["structural_diagnostics"] = len(structural)
    return {"purpose": "GFF3/CDS/serialized-score validity audit; not precision/recall evaluation",
            "inputs": {"genome": fingerprint(genome), "gff": fingerprint(gff)},
            "sequence_lengths": {name: len(sequence) for name, sequence in sequences.items()},
            "counts": dict(counts), "by_origin_and_strand": {key: dict(value) for key, value in sorted(groups.items())},
            "structural_diagnostics": structural, "failed_models": failures,
            "all_checks_pass": not structural and not failures}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in (args.genome.resolve(), args.gff.resolve()):
        parser.error("--output must differ from the input genome and GFF")
    result = validate(args.genome, args.gff)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"counts": result["counts"], "by_origin_and_strand": result["by_origin_and_strand"],
                      "all_checks_pass": result["all_checks_pass"]}, indent=2))
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
