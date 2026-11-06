# A reproducible interface for EncapsulatorP research

This is an executable first integration, not a combined theorem or a production
Validatron proof service. It keeps experiments useful while making their inputs,
scope, and missing proof obligations inspectable.

## Run locally

Install the pinned repository requirements, then run:

```bash
python -m unittest discover -s tests -v
python scripts/research_bridge.py --limit 4096 --out outputs/research_bridge.json
# Optional: use an explicitly selected tower checkout containing constructive.py.
python scripts/research_bridge.py --limit 4096 --tower-repo ../towers --out outputs/research_bridge.json
```

`--validatron-manifest path.json` records readiness from an exported manifest.
It does not contact Validatron, accept a proof attestation, or change the status
of any source mathematical claim. A fixture manifest remains fixture-only.
No model, credentials, cloud jobs, or publication are needed.

The runner accepts even limits from 4 to 100,000. It records Python/NumPy
versions, source commit and dirty status, source-file hashes, exact-count hash,
parameters and checks. Runtime and working-tree state can differ between runs;
the same arithmetic inputs must yield the same count hash. This is provenance,
not a signed independent attestation.

## Exact counts and compatibility

`scripts/exact_arithmetic.py` replaces rounded floating FFT counts with an
integer NTT modulo 998244353, primitive root 3, transform length at most 2^23.
For a binary vector of length L, every ordinary convolution coefficient is at
most L. Here L < 998244353, so each modular coefficient uniquely equals the
integer coefficient: no rounding or multi-modulus reconstruction is needed.
Intermediate products fit signed int64. The backend rejects limits above
4,194,303; the evidence runner deliberately uses a smaller budget.

For even N the unordered count is `(ordered[N] + is_prime[N/2]) // 2`.
The historical `exact_goldbach_counts_fft` API remains as a compatibility name
and delegates to the NTT. Heuristic density convolutions can still use floating
FFT. Supplying a custom prime indicator is a caller assertion, not a primality
certificate. Tests compare trial-division counts, direct integer convolution,
and modular enumeration; the capsule directly enumerates targets through 4096
and 32 deterministic samples beyond it.

## An explicit Goldbach–Catalan character connection

Let chi4(n) be zero for even n, +1 for n=1 mod 4, and -1 for n=3 mod 4.
Let rG(N) count unordered prime pairs p+q=N, and let W4(N) weight each pair
by chi4(p)chi4(q). For every even N >= 4,

    W4(N) = (-1)^(N/2+1) * (rG(N) - 1[N=4]).

Proof: outside (2,2), both primes are odd. When N=0 mod 4 their classes are
opposite and the product is -1; when N=2 mod 4 the classes match and the product
is +1. The pair (2,2) contributes zero. This elementary identity requires no
Goldbach existence assumption. No novelty or formal proof-assistant claim is
made. The weighted sum is determined by rG and N mod 4, so this specific weight
contains no additional predictive information. That is a constraint on a useful
research direction, not a rejection of richer weights.

The same character gives beta(s)=sum chi4(n)/n^s and, for Re(s)>1,

    Li_s(i) = -2^(-s) eta(s) + i beta(s).

Splitting even and odd indices proves it. Beta(2) is Catalan's constant.
The runner checks exact finite rational decompositions for s=2,3,4 and an
alternating-series rational enclosure for beta(2). The signed Ruelle–Mayer
determinant proposal in the recovered Catalan draft needs its own operator and
trace identities; evaluating one observable is not a determinant proof.

## Mertens, local constraints, and towers

The exact baseline computes mu and M(n)=sum mu(k), checks the divisor-sum
identity, and uses integer `M(n)^2 > n` comparisons. It is not a reimplementation
or audit of the unavailable Mertens counterexample source.

For m=product p^e, ordered unit pairs a+b=N mod m number

    product p^(e-1) * (p-1 if p divides N, otherwise p-2).

This follows by excluding a=0 and a=N modulo p, lifting, and CRT. It counts
local possibilities, not actual primes. This kernel continues the September
2026 EncapsulatorP development pack.

The optional tower adapter runs the real sign-flipping solver at a stated base
(1.11), with no fit to a Goldbach target. Its contraction theorem is in the
tower repository's `notes/real_tower_theorem.md`. Its float iteration estimate
excludes rounding error; it is not an interval certificate or an identity with
Brun's constant.

## Next integration obligations

1. Specify a proposed common shift: domain, state, update, observable, and exact
   equality or error bound. Distinguish Euler numbers/Dirichlet beta from other
   meanings of "Euler". No common operator has yet been established.
2. Compare richer character or tower observables on held-out ranges against
   count, target size, and modular baselines. Record parameters chosen before
   evaluation and the information required to compute them.
3. Connect real Validatron source to a pinned Lean/mathlib environment. Replay
   the exact statement independently, record allowed axioms and source lineage,
   and preserve source claim status separately from verification status.
4. Automate bounded capsules and candidate PRs only after those contracts exist.
   A model proposal, a passing finite experiment, and a checked proof are three
   different deliverables. Deployments and pipeline execution remain separate.
