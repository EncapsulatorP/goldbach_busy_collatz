# Boundary Relation

This note defines a **boundary relation** between generalized Collatz maps and
generalized Goldbach systems in the narrow sense that is actually defensible in
this repo: they can be compared through a shared modular finite-state
framework. That is a statement about method and operator structure, not about
one problem implying the other.

## Generalized Collatz Side

Take a normalized generalized Collatz map

`T_{a,b,d}(n) = (a n + b) / d^{v_d(a n + b)}`

on an admissible arithmetic domain, typically with `gcd(a, d) = gcd(b, d) = 1`
and a parity or residue restriction that keeps the map closed.

For a fixed modulus `m`, this induces a residue transition graph or transfer
operator

`K_T^(m) : f(r) -> sum_s k_T(r, s) f(s)`

on residue classes `r, s mod m`.

What matters here is not trajectory truth yet, but the modular transport
structure:

- which residue classes communicate
- with what multiplicity or weight
- whether the induced residue dynamics has drift, recurrence, or concentration

## Generalized Goldbach Side

For Goldbach, the primary arithmetic object is not a map but a representation
count. In this repo the working objects are:

- `r_G(N)`: exact unordered Goldbach count
- `h(N)`: heuristic count
- `z_h(N)`: normalized residual
- `rho_m(N) = N mod m`: residue coordinate

Fix a modulus `m` and a residual discretization `z -> z_bucket`. Then the data
induces a residue-state field

`A_G^(m)(t, r) in Z`

where `t` is a window index in `N` and `r` is a residue class. In the new
automata script, this field is realized as a cellular automaton on the `rho`
lattice whose cell state is the dominant `z_bucket` in each residue/window
cell.

This gives a Goldbach-side local operator

`K_G^(m) : neighborhood(r-1, r, r+1) -> next_state(r)`

estimated empirically from the residual field.

## Boundary Relation Definition

For a fixed modulus `m`, say that a generalized Collatz family `T_{a,b,d}` and
a generalized Goldbach family `G` are in **boundary relation at modulus `m`**
if both admit comparison through operators on the same residue state space

`R_m subseteq Z / mZ`

after projection to modular observables.

Concretely, the boundary relation consists of:

1. A shared residue index set `R_m`.
2. A Collatz operator `K_T^(m)` on `R_m`.
3. A Goldbach residue automaton or transfer operator `K_G^(m)` on `R_m`.
4. A comparison layer built from operator-level quantities such as:
   - stationary distributions
   - spectral radius or gap
   - transition entropy
   - recurrence or return structure
   - residue-class drift

Written schematically:

`T_{a,b,d}  ~_boundary,m  G`

means:

`pi_m(T_{a,b,d})` and `pi_m(G)` live on the same modular state skeleton, so
their operators can be compared.

It does **not** mean:

- Collatz truth implies Goldbach truth
- Goldbach residue effects control Collatz stopping behavior
- any theorem is transferred across the boundary

## Why The Relation Is Mostly Framework-Level

In this repo, the honest bridge is modularity plus operator methods.

- **Modularity:** both problems produce nontrivial residue-class structure at
  small moduli.
- **Framework:** both can be represented as finite-state systems after modular
  projection.
- **Methods:** both invite transfer operators, entropy summaries, and local
  neighborhood rules.

So if a similarity appears, the default interpretation should be:

`shared modular skeleton -> similar coarse operator statistics`

not:

`deep arithmetic equivalence`

## Practical Repo Consequence

This relation is strong enough to justify:

- using residue-state automata on the Goldbach side
- comparing Goldbach `rho` dynamics with Collatz residue dynamics at the level
  of transfer operators
- speaking about a shared "boundary" or "interface" language

It is not strong enough to justify:

- claiming that Goldbach residuals explain Collatz trajectories
- treating the Collatz framing as evidence for Goldbach
- treating the Goldbach automaton as a proof object

That is the intended meaning of the boundary relation in this workspace.
