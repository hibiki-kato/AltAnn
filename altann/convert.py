#!/usr/bin/env python3
"""Export all distinct gene-local k-best transcripts in genomic coordinates.

Delta is the path score minus the UniAnn best reference path score over the
same local interval, read from the required sibling TSV report. The reference
and candidate scores use 17 significant digits. References take precedence
for matching intron chains; rank 0 denotes the actual UniAnn reference.
This is a path-level score, not an independently scored transcript.
"""
from collections import Counter
import csv
from decimal import Decimal, localcontext
import json
from itertools import chain
from pathlib import Path
import re
from urllib.parse import quote, unquote

SEGMENT = re.compile(r"^(.*)__(\d+)-(\d+)-(\d+)$")


def attrs(row):
    return dict(item.split('=', 1) for item in row[8].split(';') if '=' in item)


def owned(start, end, offset, stop, total, overlap, margin):
    midpoint = (start + end) // 2
    lo = 1 if offset == 0 else offset + (overlap // 2) + 1
    hi = total if stop == total else stop - (overlap // 2)
    inner_lo = 1 if offset == 0 else offset + margin + 1
    inner_hi = total if stop == total else stop - margin
    return lo <= midpoint <= hi and start >= inner_lo and end <= inner_hi


def reference_scores(input_path):
    report = Path(input_path).with_suffix('.tsv')
    scores = {}
    with report.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'reference_gene_index', 'reference_fixed_interval_score'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Missing reference columns in {report}')
        for row in reader:
            index = int(row['reference_gene_index'])
            if index < 1:
                raise ValueError(f'Invalid reference gene index: {index}')
            locus = f'locus{index}'
            score = Decimal(row['reference_fixed_interval_score'])
            if not score.is_finite():
                raise ValueError(f'Non-finite reference score: {locus}')
            if locus in scores and scores[locus] != score:
                raise ValueError(f'Inconsistent reference score: {locus}')
            scores[locus] = score
    return scores


def convert(input_path, segment, strand, output_path, stats_path, *, chromosome,
            offset, length, total, segmented=False, k=10, short_cds_filter=False,
            overlap=4_000_000, margin=20_000, best_gff=None):
    """Map sequence-oriented paths to genomic GFF3 and retain unique intron chains.

    Coordinates in raw decoder files are one-based inclusive. Reverse-oriented
    input is reflected within its segment before adding the zero-based offset.
    """
    stop = offset + length
    references = reference_scores(input_path)
    transcripts, children, scores = {}, [], {}
    counts = Counter(input_transcripts=0, emitted_transcripts=0, duplicate_transcripts=0,
                     short_cds_transcripts=0, boundary_transcripts=0, emitted_exons=0,
                     emitted_cds=0, emitted_rank1=0, emitted_references=0,
                     emitted_local_candidates=0)
    reference_path = Path(str(input_path) + '.reference.gff')
    with reference_path.open() as ref_handle, open(input_path) as candidate_handle:
        handle = chain(ref_handle, candidate_handle)
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9:
                raise ValueError('Expected nine GFF columns')
            if row[2] not in ('transcript', 'mRNA', 'exon', 'CDS'):
                continue
            if not 1 <= int(row[3]) <= int(row[4]) <= length:
                raise ValueError(f'Coordinates outside segment: {row}')
            info = attrs(row)
            if row[2] in ('transcript', 'mRNA'):
                tid = info['ID']
                if tid in transcripts:
                    raise ValueError(f'Duplicate raw ID: {tid}')
                parts = tid.split('.')
                locus = parts[0]
                is_reference = parts[1] == 'reference'
                rank = 0 if is_reference else int(parts[1].removeprefix('k'))
                if not is_reference and not 1 <= rank <= k:
                    raise ValueError(f'Invalid local rank: {rank}')
                gene = parts[2]
                score = Decimal(row[5])
                if not score.is_finite():
                    raise ValueError('Non-finite path score')
                key = (locus, rank)
                if key in scores and scores[key] != score:
                    raise ValueError(f'Inconsistent path score: {key}')
                scores[key] = score
                transcripts[tid] = dict(row=row, exons=[], cds=[], locus=locus,
                                        rank=rank, gene=gene, score=score,
                                        origin='uniann_best' if is_reference else 'local_kbest')
            else:
                children.append((info['Parent'], row))
    for parent, row in children:
        transcripts[parent]['exons' if row[2] == 'exon' else 'cds'].append(row)
    counts['input_transcripts'] = len(transcripts)
    counts['input_loci'] = len({t['locus'] for t in transcripts.values()})
    # No-candidate loci have a sidecar reference but no TSV candidate rows.
    reference_loci = {t['locus'] for t in transcripts.values() if t['rank'] == 0}
    if not set(references).issubset(reference_loci):
        raise ValueError('TSV reference loci missing from reference GFF')
    missing = {t['locus'] for t in transcripts.values() if t['rank'] > 0} - references.keys()
    if missing:
        raise ValueError(f'Missing UniAnn reference scores: {sorted(missing)}')
    for t in transcripts.values():
        if t['rank'] == 0:
            if t['locus'] in references and t['score'] != references[t['locus']]:
                raise ValueError(f"Reference GFF/TSV score mismatch: {t['locus']}")
            references[t['locus']] = t['score']
    counts['input_reference_transcripts'] = sum(t['rank'] == 0 for t in transcripts.values())
    if best_gff:
        validate_best(best_gff, transcripts, chromosome, strand, offset, length)
    seen, output, ranks = set(), ['##gff-version 3'], Counter()
    for tid, t in sorted(transcripts.items(), key=lambda item: item[1]['rank']):
        exons = sorted((int(r[3]), int(r[4])) for r in t['exons'])
        if not exons or any(a[1] >= b[0] for a, b in zip(exons, exons[1:])):
            raise ValueError(f'Invalid exon chain: {tid}')
        introns = tuple((a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:]))
        key = (t['locus'], introns)
        if key in seen:
            counts['duplicate_transcripts'] += 1
            continue
        seen.add(key)
        if short_cds_filter and len(t['exons']) <= 2 and sum(int(r[4]) - int(r[3]) + 1 for r in t['cds']) <= 200:
            counts['short_cds_transcripts'] += 1
            continue
        def transform(row):
            row = row.copy()
            start, end = int(row[3]), int(row[4])
            if strand == 'minus':
                start, end = length - end + 1, length - start + 1
            row[0], row[3], row[4] = quote(chromosome, safe='._-:'), str(start + offset), str(end + offset)
            row[5], row[6] = '.', '-' if strand == 'minus' else '+'
            return row
        row = transform(t['row'])
        if segmented and not owned(int(row[3]), int(row[4]), offset, stop, total, overlap, margin):
            counts['boundary_transcripts'] += 1
            continue
        reference = references[t['locus']]
        with localcontext() as ctx:
            ctx.prec = max(50, len(t['score'].as_tuple().digits) + len(reference.as_tuple().digits) +
                           abs(t['score'].as_tuple().exponent - reference.as_tuple().exponent) + 5)
            delta = t['score'] - reference
        new_id = quote(f'{segment}.{strand}.kbest.{tid}', safe='._-')
        row[8] = ';'.join(f'{k}={v}' for k, v in [
            ('ID', new_id), ('kbest_score', format(t['score'], 'f')),
            ('kbest_reference', 'uniann_best'), ('kbest_origin', t['origin']),
            ('kbest_reference_score', format(reference, 'f')), ('kbest_delta', format(delta, 'f')),
            ('kbest_rank', t['rank']), ('kbest_locus', quote(f"{segment}.{strand}.kbest.{t['locus']}", safe='._-')),
            ('kbest_gene_in_path', quote(t['gene'], safe='._-'))])
        output.append('\t'.join(row))
        for kind, rows in [('exon', t['exons']), ('CDS', t['cds'])]:
            for index, child in enumerate(rows, 1):
                child = transform(child)
                child[8] = f'ID={new_id}.{kind}{index};Parent={new_id}'
                output.append('\t'.join(child))
        counts['emitted_transcripts'] += 1
        counts['emitted_rank1'] += t['rank'] == 1
        counts['emitted_references'] += t['rank'] == 0
        counts['emitted_local_candidates'] += t['rank'] > 0
        counts['emitted_exons'] += len(t['exons'])
        counts['emitted_cds'] += len(t['cds'])
        ranks[str(t['rank'])] += 1
    Path(output_path).write_text('\n'.join(output) + '\n')
    result = dict(counts, segment=segment, strand=strand, emitted_by_rank=dict(ranks))
    Path(stats_path).write_text(json.dumps(result, indent=2) + '\n')
    return result


def validate_best(path, transcripts, chromosome, strand, offset, length):
    """Require supplied (possibly filtered) GFF models to match the baseline.

    UniAnn output may use local forward coordinates even for a reverse input.
    Accept that convention or already mapped genomic coordinates, but require
    every supplied model to use one consistent coordinate convention.
    """
    models, children, coding = {}, [], []
    populated = False
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            if len(row) != 9:
                raise ValueError(f'Expected nine GFF columns in {path}')
            populated = True
            info = attrs(row)
            if row[2] in ('mRNA', 'transcript'):
                models[info['ID']] = []
            elif row[2] == 'exon':
                for parent in info['Parent'].split(','):
                    children.append((parent, row))
            elif row[2] == 'CDS':
                parents = info.get('Parent', '').split(',')
                if not all(parents):
                    raise ValueError(f'CDS has no Parent in {path}')
                for parent in parents:
                    coding.append((parent, row))
    # The original UniAnn binary emits only CDS records, grouped by gene
    # Parent. Its wrapper adds transcripts and exons. Validate either form
    # against the corresponding baseline features without inventing exons.
    cds_only = not models and bool(coding)
    if cds_only:
        models = {parent: [] for parent, _ in coding}
        if children:
            raise ValueError(f'Exon has no transcript in {path}')
        children = coding
    for parent, row in children:
        if parent not in models:
            raise ValueError(f'Exon has no transcript in {path}: {parent}')
        models[parent].append((unquote(row[0]), int(row[3]), int(row[4]), row[6]))
    if not models:
        if populated:
            raise ValueError(f'Best GFF has no transcript or CDS models: {path}')
        return  # An entirely filtered UniAnn GFF is a valid empty subset.
    if any(not exons for exons in models.values()):
        raise ValueError(f'Best GFF validation requires exon features: {path}')
    local, genomic = set(), set()
    local_seqids = {t['row'][0] for t in transcripts.values() if t['rank'] == 0}
    for t in transcripts.values():
        if t['rank'] != 0:
            continue
        features = t['cds'] if cds_only else t['exons']
        exons = sorted((int(r[3]), int(r[4])) for r in features)
        local.add(tuple(exons))
        mapped = [(length-b+1, length-a+1) if strand == 'minus' else (a,b)
                  for a,b in exons]
        genomic.add(tuple(sorted((chromosome, a+offset, b+offset,
                                   '-' if strand == 'minus' else '+') for a,b in mapped)))
    matches_local = all(tuple(sorted((a,b) for _,a,b,_ in e)) in local
                        and all(c in local_seqids and s in ('+', '.') for c,_,_,s in e)
                        for e in models.values())
    matches_genomic = all(tuple(sorted(e)) in genomic for e in models.values())
    if not (matches_local or matches_genomic):
        raise ValueError(f'Supplied best GFF differs from the reconstructed baseline: {path}')
