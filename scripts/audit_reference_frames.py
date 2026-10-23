#!/usr/bin/env python3
"""Audit CDS biology and oriented frame coordinates from an existing reference.

This validates annotation conventions; it does not produce model predictions
or sensitivity/precision measurements. All coordinates in the input GFF are
one-based inclusive, while stored internal intervals are zero-based half-open.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote


COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
STOPS = {'TAA', 'TAG', 'TGA'}


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def audit(genome, reference, chromosomes=None):
    sequences, identifier, chunks = {}, None, []
    with genome.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if identifier:
                    sequences[identifier] = ''.join(chunks).upper()
                identifier, chunks = line[1:].split()[0], []
            else:
                chunks.append(''.join(line.split()))
    if identifier:
        sequences[identifier] = ''.join(chunks).upper()
    selected = set(chromosomes) if chromosomes else set(sequences)
    if selected - sequences.keys():
        raise ValueError('Unknown chromosome: ' + ', '.join(sorted(selected - sequences.keys())))
    metadata, coding = {}, defaultdict(list)
    with reference.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            attrs = {k: unquote(v) for k, v in
                     (part.split('=', 1) for part in row[8].split(';') if '=' in part)}
            if 'ID' in attrs:
                metadata[attrs['ID']] = attrs
            if row[2] == 'CDS':
                for parent in attrs['Parent'].split(','):
                    coding[parent].append((row[0], row[6], int(row[3]) - 1,
                                           int(row[4]), int(row[7]), attrs))
    counts, exclusions, failures, examples, valid_cases = Counter(), Counter(), [], [], []
    frame_histogram = Counter()
    for transcript, records in coding.items():
        if records[0][0] not in selected:
            exclusions['out_of_scope'] += 1
            continue
        attributes = [metadata.get(transcript, {}), *(record[5] for record in records)]
        for parent in metadata.get(transcript, {}).get('Parent', '').split(','):
            attributes.append(metadata.get(parent, {}))
        reason = next((key for key in ('partial', 'pseudo', 'exception')
                       if any((bool(attrs.get(key)) if key == 'exception'
                               else attrs.get(key) == 'true') for attrs in attributes)), None)
        if reason:
            exclusions[reason] += 1
            continue
        chrom, strand = records[0][:2]
        sequence = sequences[chrom]
        length = len(sequence)
        records.sort(key=lambda record: record[2], reverse=strand == '-')
        errors, pieces, cumulative, intron_bases, probes = [], [], 0, 0, 0
        oriented = []
        for _, _, start, end, phase, _ in records:
            lo, hi = (start, end) if strand == '+' else (length - end, length - start)
            oriented.append((lo, hi, phase))
        initial_frame = oriented[0][0] % 3
        for index, (record, (lo, hi, phase)) in enumerate(zip(records, oriented)):
            if phase != (-cumulative) % 3:
                errors.append('phase_discontinuity')
            if index:
                gap = lo - oriented[index - 1][1]
                if gap < 0:
                    errors.append('overlapping_CDS')
                intron_bases += gap
            expected_frame = (initial_frame + intron_bases) % 3
            # Labels come from the recorded GFF phase, independently of the
            # ATG coordinate and physical intron lengths defining the E frame.
            for offset in sorted({0, min(1, hi - lo - 1), min(2, hi - lo - 1),
                                  (hi - lo) // 2, hi - lo - 1}):
                local_phase = ((3 - phase) + offset) % 3
                if (lo + offset - local_phase) % 3 != expected_frame:
                    errors.append('frame_disagreement')
                probes += 1
            piece = sequence[record[2]:record[3]]
            pieces.append(piece if strand == '+' else piece.translate(COMPLEMENT)[::-1])
            cumulative += hi - lo
        cds = ''.join(pieces)
        if not cds.startswith('ATG'):
            errors.append('non_ATG_start')
        if cds[-3:] not in STOPS:
            errors.append('missing_terminal_stop')
        if len(cds) % 3:
            errors.append('CDS_length_not_divisible_by_3')
        if any(cds[index:index + 3] in STOPS for index in range(0, len(cds) - 3, 3)):
            errors.append('internal_stop')
        counts['eligible_' + strand] += 1
        if errors:
            failures.append({'transcript': transcript, 'chromosome': chrom,
                             'strand': strand, 'errors': sorted(set(errors))})
            continue
        counts['valid_' + strand] += 1
        counts['validated_CDS_segments'] += len(records)
        counts['frame_probes'] += probes
        counts['old_mapping_disagrees_at_probes'] += probes
        frame_histogram[f'{strand}:E{initial_frame}'] += 1
        case = {'transcript': transcript, 'chromosome': chrom, 'strand': strand,
                'chromosome_length': length, 'start_frame': initial_frame,
                'CDS_length': len(cds), 'start_codon': cds[:3], 'stop_codon': cds[-3:],
                'CDS_1based': [[record[2] + 1, record[3], record[4]] for record in records],
                'oriented_CDS_0based': [list(interval) for interval in oriented]}
        if len(examples) < 12 or transcript in ('rna-NM_001018182.4', 'rna-NM_001355965.2'):
            examples.append(case)
        key = strand, initial_frame
        if (key not in {(item['strand'], item['start_frame']) for item in valid_cases}
                or transcript == 'rna-NM_001355965.2'):
            low, high = min(record[2] for record in records), max(record[3] for record in records)
            segment = sequence[low:high]
            if strand == '-':
                segment = segment.translate(COMPLEMENT)[::-1]
            prefix = 30 + oriented[0][0] % 3
            case = dict(case, prefix_bases=prefix, sequence='C' * prefix + segment + 'C' * 30,
                        local_CDS_0based=[[lo - oriented[0][0] + prefix,
                                          hi - oriented[0][0] + prefix, phase]
                                         for lo, hi, phase in oriented])
            valid_cases.append(case)
    result = {'purpose': 'reference coordinate audit, not model evaluation',
              'inputs': {name: {'path': str(path.resolve()), 'sha256': sha256(path)}
                         for name, path in (('genome', genome), ('reference', reference))},
              'sequence_lengths': {name: len(seq) for name, seq in sequences.items()},
              'chromosome_scope': sorted(selected),
              'counts': dict(counts), 'excluded': dict(exclusions),
              'start_frames': dict(frame_histogram), 'failed_models': failures,
              'verified_examples': examples,
              'correct_mapping': '(oriented_position - local_codon_offset) % 3',
              'old_mapping': '(oriented_position - local_codon_offset - 1) % 3'}
    return result, valid_cases


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fixture', type=Path)
    parser.add_argument('--chromosome', action='append', help='Repeat to restrict the audit scope.')
    args = parser.parse_args()
    result, cases = audit(args.genome, args.reference, args.chromosome)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    if args.fixture:
        args.fixture.parent.mkdir(parents=True, exist_ok=True)
        args.fixture.write_text(json.dumps({'inputs': result['inputs'], 'cases': cases}, indent=2) + '\n')
    print(json.dumps({'counts': result['counts'], 'excluded': result['excluded'],
                      'failed_models': len(result['failed_models'])}))
