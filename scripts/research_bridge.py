"""Bounded, local evidence runner for the EncapsulatorP arithmetic interface.

No model calls, remote jobs, credentials, publication, or theorem promotion.
Run: python scripts/research_bridge.py --limit 4096 --tower-repo ../towers
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from exact_arithmetic import (NTT_MODULUS, NTT_ROOT, chi4, exact_goldbach_counts,
    integer, mobius_prefix, prime_indicator, signed_prime_pair_sum, unit_pair_count)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_state(root, paths):
    root = Path(root).resolve()
    def git(*args):
        result = subprocess.run(['git', '-C', str(root), *args], capture_output=True,
                                text=True, timeout=5, check=True)
        return result.stdout.strip()
    try:
        commit = git('rev-parse', 'HEAD')
        dirty = bool(git('status', '--porcelain', '--untracked-files=normal'))
    except (OSError, subprocess.SubprocessError):
        commit, dirty = None, None
    return {'commit': commit, 'dirty': dirty,
            'sha256': {p: digest(root / p) for p in paths}}


def beta_interval(terms, weight=2):
    """Rational alternating-series enclosure for beta(weight), weight >= 1."""
    terms = integer(terms, 'terms', 1, 1024)
    weight = integer(weight, 'weight', 1, 8)
    partial = sum((Fraction((-1)**k, (2*k+1)**weight) for k in range(terms)), Fraction())
    adjacent = partial + Fraction((-1)**terms, (2*terms+1)**weight)
    return min(partial, adjacent), max(partial, adjacent)


def character_identity(terms, weight):
    """Exact finite Li_s(i) split, with 2*terms summands and rational parts."""
    terms = integer(terms, 'terms', 1, 1024)
    weight = integer(weight, 'weight', 1, 8)
    real = sum((Fraction((0, -1, 0, 1)[(n-1) % 4], n**weight)
                for n in range(1, 2*terms+1)), Fraction())
    imag = sum((Fraction(chi4(n), n**weight)
                for n in range(1, 2*terms+1)), Fraction())
    eta = sum((Fraction((-1)**(n-1), n**weight) for n in range(1, terms+1)), Fraction())
    beta = sum((Fraction((-1)**k, (2*k+1)**weight) for k in range(terms)), Fraction())
    return {'weight': weight, 'terms': 2*terms,
            'real_matches_minus_eta_over_two_power_s': real == -eta / 2**weight,
            'imag_matches_beta': imag == beta}


def mertens_summary(limit):
    mu, prefix = mobius_prefix(limit)
    # Compare squared ratios with integer cross products, not floating sqrt.
    peak = 1
    violations = []
    for n in range(1, limit+1):
        if int(prefix[n])**2 * peak > int(prefix[peak])**2 * n:
            peak = n
        if int(prefix[n])**2 > n:
            violations.append(n)
    # Independent divisor-sum identity checks for a fixed bounded prefix.
    checked = min(limit, 512)
    inverse_ok = all(sum(int(mu[d]) for d in range(1,n+1) if n % d == 0) == int(n == 1)
                     for n in range(1, checked+1))
    return {'limit': limit, 'M_at_limit': int(prefix[-1]), 'peak_n': peak,
            'peak_M_squared_over_n': str(Fraction(int(prefix[peak])**2, peak)),
            'M_squared_gt_n_count': len(violations), 'first_violation': next(iter(violations), None),
            'mobius_divisor_identity_checked_through': checked, 'identity_passed': inverse_ok,
            'status': 'exact_finite_computation', 'universal_bound_claimed': False}


def formal_readiness(manifest):
    """Readiness observation only; never accepts a self-reported proof."""
    if not manifest:
        return {'status': 'unavailable', 'formal_proof_verified': False}
    if not isinstance(manifest, dict) or not isinstance(manifest.get('backends', {}), dict):
        raise ValueError('manifest and its backends must be JSON objects')
    backends = manifest.get('backends', {})
    fixture = (manifest.get('trust_level') == 'test_fixture'
               or 'fixture' in str(manifest.get('corpus_version', '')).lower()
               or 'fixture' in str(backends.get('kernel_runner', '')).lower())
    return {'status': 'fixture_only' if fixture else 'requires_independent_kernel_replay',
            'observed_corpus_version': manifest.get('corpus_version'),
            'reported_production_ready': manifest.get('production_ready'),
            'formal_proof_verified': False,
            'source_claim_status_changed': False}


def run(limit=4096, tower_repo=None, validatron_manifest=None):
    limit = integer(limit, 'limit', 4, 100_000)
    if limit % 2:
        raise ValueError('limit must be even')
    start = time.perf_counter()
    flags = prime_indicator(limit)
    counts = exact_goldbach_counts(limit, flags)
    checked_targets = list(range(4, min(limit, 4096)+1, 2))
    # Above the exhaustive prefix, select a deterministic coverage sample.
    if limit > 4096:
        checked_targets = sorted(set(checked_targets +
            [2*int(x) for x in np.linspace(2049, limit//2, 32)]))
    count_checks = []
    signed_checks = []
    for n in checked_targets:
        direct = sum(bool(flags[p] and flags[n-p]) for p in range(2, n//2+1))
        count_checks.append(direct == int(counts[n]))
        expected = (-1)**(n//2+1)*(int(counts[n])-int(n == 4))
        signed_checks.append(signed_prime_pair_sum(n, flags) == expected)
    root = Path(__file__).resolve().parents[1]
    report = {
        'schema_version': 'encapsulatorp-evidence-1',
        'purpose': 'Reproducible arithmetic kernels and explicit cross-scope proof obligations',
        'parameters': {'limit': limit, 'oracle_exhaustive_through': min(limit, 4096),
                       'oracle_target_count': len(checked_targets)},
        'sources': {'goldbach': source_state(root, ['scripts/exact_arithmetic.py',
                                                    'scripts/research_bridge.py'])},
        'environment': {'python': sys.version.split()[0], 'numpy': np.__version__},
        'goldbach': {'method': 'integer_ntt', 'modulus': NTT_MODULUS, 'primitive_root': NTT_ROOT,
                    'coefficient_bound': limit+1, 'coefficient_bound_below_modulus': limit+1 < NTT_MODULUS,
                    'counts_sha256': hashlib.sha256(counts.astype('<i8').tobytes()).hexdigest(),
                    'oracle_checks_passed': all(count_checks), 'count_at_limit': int(counts[limit]),
                    'zero_even_targets_in_range': int(np.count_nonzero(counts[4::2] == 0)),
                    'universal_conjecture_claimed': False},
        'modular': {'modulus': 30, 'ordered_unit_pair_counts':
                    {str(n): unit_pair_count(n, 30) for n in range(0, 30, 2)},
                    'actual_prime_existence_implied': False},
        'character_bridge': {'identity': 'W4(N)=(-1)^(N/2+1)*(rG(N)-1[N=4])',
                             'direct_checks_passed': all(signed_checks),
                             'independent_predictive_signal': False},
        'catalan': {'finite_character_decompositions': [character_identity(128,s) for s in (2,3,4)],
                    'beta2_enclosure': [str(x) for x in beta_interval(128)],
                    'enclosure_method': 'exact_rational_alternating_series',
                    'brun_tower_identity_claimed': False},
        'mertens': mertens_summary(limit),
        'validatron': formal_readiness(validatron_manifest),
    }
    if tower_repo is not None:
        tower_repo = Path(tower_repo).resolve()
        sys.path.insert(0, str(tower_repo / 'scripts'))
        solver = importlib.import_module('boundary_residues.constructive')
        if Path(solver.__file__).resolve() != tower_repo/'scripts/boundary_residues/constructive.py':
            raise ValueError('tower module resolved outside the selected checkout')
        report['sources']['towers'] = source_state(tower_repo, [
            'scripts/boundary_residues/constructive.py', 'notes/real_tower_theorem.md'])
        report['tower'] = {'base': 1.11, **solver.solve_real_signflip(1.11),
                           'parameter_fitted_to_goldbach_target': False}
    else:
        report['tower'] = {'status': 'not_run', 'reason': 'explicit tower checkout required'}
    checks = [all(count_checks), all(signed_checks), report['mertens']['identity_passed']]
    checks += [v for item in report['catalan']['finite_character_decompositions']
               for k,v in item.items() if k.endswith('matches_beta') or k.startswith('real_matches')]
    report['computation_checks_passed'] = all(checks)
    report['formal_proof_verified'] = False
    report['elapsed_seconds'] = time.perf_counter()-start
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=4096)
    parser.add_argument('--tower-repo', type=Path)
    parser.add_argument('--validatron-manifest', type=Path)
    parser.add_argument('--out', type=Path, default=Path('outputs/research_bridge.json'))
    args = parser.parse_args()
    manifest = json.loads(args.validatron_manifest.read_text()) if args.validatron_manifest else None
    if manifest is not None and not isinstance(manifest, dict):
        parser.error('manifest must be a JSON object')
    report = run(args.limit, args.tower_repo, manifest)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(args.out), 'checks_passed': report['computation_checks_passed'],
                      'formal_proof_verified': False, 'elapsed_seconds': report['elapsed_seconds']}))
    return 0 if report['computation_checks_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
