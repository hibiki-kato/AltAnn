"""Validate and serialize opt-in adaptive K policies."""
import json
import math
from pathlib import Path


def load_adaptive_k_policy(path):
    """Return the original policy object and its native decoder specification."""
    try:
        policy = json.loads(Path(path).read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f'Invalid adaptive K policy JSON: {exc}') from exc
    if not isinstance(policy, dict):
        raise ValueError('Adaptive K policy must be a JSON object')
    family = policy.get('family')
    family_fields = {'linear': ('slope',), 'shifted_exponential': ('scale', 'rate')}
    if not isinstance(family, str) or family not in family_fields:
        raise ValueError('Unsupported adaptive K policy family')
    allowed = {'schema', 'family', 'base', 'cap', 'feature', 'rounding', 'status', 'source'} | set(family_fields[family])
    unknown = set(policy) - allowed
    if unknown:
        raise ValueError(f'Unknown adaptive K policy fields: {", ".join(sorted(unknown))}')
    if policy.get('schema') != 'dna-adaptive-k-v1':
        raise ValueError('Unsupported adaptive K policy schema')
    if policy.get('feature', 'global_UniAnn_CDS_exon_count') != 'global_UniAnn_CDS_exon_count':
        raise ValueError('Unsupported adaptive K policy feature')
    if policy.get('rounding', 'ceil') != 'ceil':
        raise ValueError('Unsupported adaptive K policy rounding')
    cap = policy.get('cap')
    if type(cap) is not int or not 1 <= cap <= 32767:
        raise ValueError('Adaptive K policy cap must be an integer between 1 and 32767')
    def number(name, minimum):
        value = policy.get(name)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f'Adaptive K policy {name} must be a finite number >= {minimum}')
        try:
            value = float(value)
        except OverflowError as exc:
            raise ValueError(f'Adaptive K policy {name} must be finite') from exc
        if not math.isfinite(value) or value < minimum:
            raise ValueError(f'Adaptive K policy {name} must be a finite number >= {minimum}')
        return repr(value)
    parts = ['linear' if family == 'linear' else 'shifted-exp', number('base', 1)]
    parts.extend(number(name, 0) for name in family_fields[family])
    return policy, ','.join(parts)
