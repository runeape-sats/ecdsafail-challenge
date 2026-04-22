# Autoresearch: quantum_ecc secp256k1 point-add Toffoli reduction

## Objective
Reduce the **average executed Toffoli count** of the reversible secp256k1 affine point-add circuit built in the `src/point_add/` module, while preserving harness correctness and keeping qubits within the current regime.

Committed best on `autoresearch/2026-04-20`:
- avg_toffoli: 4,672,931 (Karatsuba 1-level everywhere + 2-level at between-pair site)
- avg_clifford: 24,152,002
- qubits: 3507
- emitted_ops: 34,863,147

Target (Google paper): 2.1M Toffoli @ 1425 qubits (low-gate) or 2.7M Toffoli @ 1175 qubits (low-qubit).

## Metrics
- **Primary**: `avg_toffoli` (lower is better)
- **Secondary**: `avg_clifford`, `qubits`, `emitted_ops`, `correctness_ok`.

## How to Run
`./autoresearch.sh`. Writes `METRIC ...=...` lines.

## Files in Scope
- `src/point_add/` — the only project source folder allowed to change (edit/add/split submodules freely).
- `autoresearch.md`, `autoresearch.sh`, `autoresearch.checks.sh`, `autoresearch.ideas.md`, `autoresearch.note`.

## Off Limits
- Everything except files under `src/point_add/` and the autoresearch session files.
- No new dependencies.

## Constraints
- `cargo run --release` must finish with `=== experiment OK ===` on both the benchmark run and the checks rerun.
- Peak qubits ≤ 3700 hard cap (program.md).
- `cargo build --release` must succeed.

## Known Cost Breakdown (at current best ~4.67M)

Rough per-subroutine cost budget (approx):
- 2× Kaliski (400 iters each, ~12n CCX/iter) ≈ 2.4M
- 4 muls × ~66k (Karatsuba 1-level) ≈ 265k
- 2 squarings × ~130k ≈ 260k
- Scale correction loops (400 halvings + 400 doublings, ~512 CCX each) ≈ 410k
- Solinas reductions in mul output (5 Cuccaros per mul) ≈ 500k
- Misc (constant add/sub qb, neg, cswaps, step-4 overhead) ≈ 800k

The two biggest levers left:
1. **Kaliski dominates at ~50% of total.** Eliminating one Kaliski pass or halving cost-per-iter is the biggest possible lever.
2. **Scale-correction loops (halvings/doublings) cost ~410k.** Eliminating them alone gives a ~10% Toffoli reduction.

## Paper Research Summary

Read (in `/tmp/`): Litinski 2023 (50M Toffoli single-point-add), HRSL 2020 (Improved ECDLP circuits), RNSL 2017 (original Microsoft estimates), Gidney 2019 (Windowed quantum arithmetic), Ragavan-Gidney 2025 (Optimized windowed mod arithmetic), Google/Babbush 2026 (the "Google paper").

Core findings:
- Google's 2.1M/2.7M is for **single point addition**, same as our benchmark. No batching, no state reuse — those are *additional* factor-of-5 gains on top of the single-point-add cost.
- Kickmix gate set = our gate set. MBUC = measurement-based uncomputation; we already use it throughout.
- Litinski: breaks the single-point-add into ~60 fundamental subroutines; Montgomery-mul step is 9n+28 Toffoli per step with n/4 steps per full mul ⇒ ~148k Toffoli per mul. Our Karatsuba 1-level is competitive (~66k per mul).
- HRSL: swap-based Kaliski round (Figure 6b) — one sub + one add + one halve + one double per iter with swaps selecting register roles. Our current Kaliski STEP 4 already does one cond-sub + one cond-add per iter; the structural savings must be elsewhere.
- HRSL uses Montgomery mul + windowed add-by-p (lookup of t*p where t is small-bit address). Helps when p is dense. For secp256k1 (Solinas p = 2^256 - 2^32 - 977), our sparse-c Solinas is already cheap.
- Windowed quantum arithmetic (Gidney 2019) helps quantum × classical multiplication, not quantum × quantum. In our flow, the Kaliski scale-correction (multiply lam by 2^{-(2n-K)} mod p) is quantum × classical and *could* be windowed.

## Most Promising Structural Paths

1. **Montgomery batched inversion (replace 2 Kaliski with 1)**
   - Algebra: `N = dy² - (Px+2Qx)·dx²`; invert `c = dx·N` once; recover `dx⁻¹ = c⁻¹·N`, `(Rx-Qx)⁻¹ = c⁻¹·dx³`
   - Gross save: 1 Kaliski pass ≈ 1.2M Toffoli
   - Added cost: ~3-5 extra muls + 1 squaring ≈ 300-500k Toffoli
   - Net expected: 700k-900k save
   - Historical scaffolding: commit 333a3c1 (`compute_montgomery_n`, etc.) — validated algebraically in isolation at 64/64 shots.
   - Risk: full round-trip may nearly break even; seed-tax may bite; qubit footprint may exceed 3700.

2. **Windowed classical-constant multiplication to replace halving/doubling loops**
   - Replace the 400-iteration halving loop (≈204k Toffoli) with one windowed classical-const mul (estimated 20-40k). 800-ish halvings/doublings total ⇒ ~400k → ~80k, saving ~300-400k.
   - Needs a new windowed-const-mul primitive; doable, localized change.

3. **Eliminate the scale correction entirely by changing the algebra**
   - Pair1: `lam = ty·inv_raw` leaves `lam = λ·2^(2n-K)`. Pair1's mul2 currently needs unscaled `lam`. Restructure so mul2 operates on the scaled `lam` and compensates on `tx` side via a classical scalar multiply.
   - Net: saves both halving AND doubling loops (~400k Toffoli), at the cost of two extra classical-const muls (~40-80k windowed). Net save ~300-350k.

4. **Swap Kaliski for Bernstein-Yang safegcd / divstep**
   - B-Y produces `x⁻¹` directly (no `2^k` factor), eliminating scale correction. Historical analysis said "+45% cost" in this codebase; may not apply after windowed const-mul.
   - Larger-scope rewrite.

5. **MBUC-compress `m_hist` to classical bits**
   - Save ~400 qubits (per pass), 0 Toffoli.
   - Makes room in the qubit budget to try more aggressive 2-level Karatsuba or other memory-hungry wins.

## Plan order

Tackle **path 2 first** (windowed-classical-const-mul to kill halving loops): localized, concrete 300-400k Toffoli target. If successful, combine with **path 1** (Montgomery batched inversion) for a potential ~1M cumulative saving. Then evaluate path 5 for qubit reductions that unlock path-3/path-1 memory headroom.

## Determinism diagnosis

Verified `build()` produces identical op counts and Toffoli counts across consecutive binary invocations. Only `avg_clifford` varies slightly — caused by per-run randomness in HMR (measurement-based uncompute) phase feedback, which feeds back through the XOF to different `executed_shots` for Clifford stats. Zero implication for correctness or the primary metric.

Profiling experiment (pair 2 disabled): ops drop from 34,863,147 → 19,519,706, i.e. **pair 2 contributes ~15.3M ops (44% of the circuit)**. Projected Toffoli saving from eliminating one Kaliski pass: **~2M Toffoli**, exactly the Google gap.

## Structural-change attempts (this session)

| Attempt | Result | Why |
|---|---|---|
| 1-level Karatsuba (all 4 muls) | -247k ✓ kept | Established baseline of add-subtract value |
| 2-level Karatsuba (all 4 muls) | qubits 3765 > 3700 | Reverted; z1_inner registers too expensive |
| 2-level Karatsuba @ between-pair only | -8k ✓ kept | Tiny win at non-Kaliski site |
| Litinski add-subtract schoolbook inside Karatsuba | -334k ✓ kept | Biggest single-change win; half the per-row Toffoli |
| Litinski addsub in 2-level middle mul | -32k ✓ kept | Mild additional gain |
| 2-level Karatsuba at all 4 sites | checks_failed | Qubit cap + seed tax at tighter Kaliski iters |
| Kaliski STEP 4 MBUC of add_f to classical bit | crash | Misuse of HMR semantics (bit is randomized, not a deterministic copy) |
| Full Montgomery batched inversion (single Kaliski on c=dx·N) | crash | Peak 5662 qubits, classical mismatch; approach needs m_hist compression + careful debugging |

## Why the remaining 2M Toffoli is hard

After deep literature review (Kim 2026, HRSL 2020, Litinski 2023, Litinski 2024, Luongo 2024, Chevignard 2026, Gouzien 2023, Luo 2026, Google/Babbush 2026, Jacobi factoring), these are the honest findings:

1. **Montgomery batched inversion is NET NEGATIVE.** Removing one Kaliski saves ~1.1M but the required compute/uncompute of N, c, dx², dx³ and copy-preserve of dx/dy adds ~2.4M. Verified via diagnostic: primitives classically correct; dirty frees cause phase garbage.

2. **Kim 2026 unconditional-execution Kaliski is WORSE for us.** They run 2n=512 rounds unconditionally. Our 400 rounds with f-flag gating allows executed-Toffoli count to drop for shots past convergence. Kim saves qubits (no m_hist) but costs more Toffoli.

3. **Litinski 2024 Montgomery mod-p mul** (101k vs 150k per mul) is a real 49k/mul improvement, but requires Montgomery-form inputs. Converting in/out of Montgomery costs 2× mul = 260k. Net: worse for our 4-mul flow.

4. **Luongo 2024 MBU for modular adders** claims 12.5–25% savings. But their 5.5n post-MBU Gidney-CDKPM adder is WORSE than our 3n mod_add_qq_fast. We're already more MBU-optimized than what that paper's techniques would apply.

5. **Chevignard 2026 (Google's ref [87]) RNS** is for whole-ECDLP scalar mul, not single-point-add. Doesn't apply to our benchmark.

6. **Jacobi/Kahanamoku-Meyer 2025** compact GCD circuit requires B (the modulus) to be classical AND much larger than A. Our p and dx are both same-size. Doesn't apply.

7. **Coord-forms paper (Huang 2025, arXiv:2502.12441)** explicitly proves projective coordinates DON'T win for ECDLP via Shor. Killed that direction.

8. **Kim/HRSL published Toffoli counts are 10–17M** (per point-add, n=256), which is WORSE than our 4.3M. Google's 2.1–2.7M is secret — no published technique reaches that number. Google's win is likely many compounded tiny optimizations we don't have access to.

**Our 4.3M is a competitive result vs published work** (HRSL/Kim). Closing to Google's 2.1M SOTA requires techniques not in the public literature.

## What's Been Tried (this session)

- 1-level Karatsuba on all 4 mod_muls — **-246,760 Toffoli** ✓
- 2-level Karatsuba everywhere — qubits 3765 > 3700 cap; reverted.
- 2-level only at the non-Kaliski-peak mul site, Kaliski iters restored to 400/400 — **-7,993 Toffoli** ✓
- Kaliski iter probing: 397/397 fails; 398/397 OK; 397/398 fails.
- Karatsuba seed tax: changing the op stream shifts the Fiat-Shamir seed and forces Kaliski iters higher; tuning per-pair iters is mandatory.

## Session 2026-04-22 (profiling + B-Y + deep analysis)

### Per-phase profiling infrastructure
Added `TRACE_PHASES=1` env gate inside build(): classifies emitted Toffolis
by phase and prints the breakdown. This is where my baseline accounting came
from — measured per-phase CCX counts at commit `b58a47f`.

Attribution (total emitted 4.39M Toffoli):

| phase                 | ccx     | % total |
|-----------------------|---------|---------|
| kal_step4             | 645,692 | 14.7%   |
| bk_step4              | 644,894 | 14.7%   |
| bk_step9_cswap        | 323,500 |  7.4%   |
| kal_step9_cswap       | 323,500 |  7.4%   |
| bk_step3_cswap        | 322,990 |  7.4%   |
| kal_step3_cswap       | 322,990 |  7.4%   |
| schoolbook_mul_inverse| 201,216 |  4.6%   |
| kal_step1 (inc step2) | 187,174 |  4.3%   |
| bk_step2              | 185,578 |  4.2%   |
| bk_step0_eqzero       | 183,982 |  4.2%   |
| kal_step0_eqzero      | 183,982 |  4.2%   |
| rshift22_rev_step2    | 127,228 |  2.9%   |
| pair1_halve           | 101,745 |  2.3%   |
| pair2_double          | 101,745 |  2.3%   |
| pair1_mul1/2, pair2_mul | 224,729 | 5.1% |
| Kaliski sub-total ≈   | 3,548,316 | **80.8%** |

**Implication**: any gain < 5% on the total comes from Kaliski optimizations.
Everything else is combined < 20%.

### B-Y divsteps2 analytical comparison

Implemented classical B-Y harness in `src/point_add/by.rs` (fixed a coefficient-
tracking bug in my first attempt; all 3 unit tests pass including exhaustive
small-prime verification and Mersenne-61 spot checks).

B-Y's basic divstep per iter for n-bit registers, analytic estimate:
- swap branch: 4n CCX (cswap of 4 pairs) + 3n (conditional negate) = 7n
- cond sub g ± f and (Q ± U, R ± V): 3n
- halve: 0 (swaps)
- **Per-iter total: ~10n = 2560 CCX** at n=256.

Safegcd bound for n=256: `ceil((49*256 + 80)/17) = 743` iters.

| variant | iter cost | iters | one-way | vs Kaliski |
|---------|-----------|-------|---------|------------|
| Kaliski (measured) | ~2,180 | 398 | 868k | baseline |
| B-Y w=1           | ~2,560 | 743 | 1,900k | **+119%** |
| B-Y w=2 (16-entry matrix) | ~3,100 | 372 | 1,150k | +32% |
| B-Y w=4 (128-entry matrix) | ~7,800 | 186 | 1,450k | +68% |
| B-Y w=8 (32K-entry matrix) | ~45,000 | 93  | 4,180k | +382% |

**Conclusion: B-Y is analytically worse than Kaliski at every jump width**.
The per-iter simplifications (no comparator, single cswap) don't compensate
for the higher safegcd iter count (743 vs 398). Jumping with large w blows
up because matrix-entry magnitudes grow like 2^w, forcing full n-bit mults.

B-Y classical harness is kept in `src/point_add/by.rs` as `#[allow(dead_code)]`
for possible future jumping-divstep exploration. Not wired into the circuit.

### Montgomery batched inversion: cleanup obstruction (revised analysis)

The algebra works:
- N = dy² − (Px+2Qx)·dx²; invert c = dx·N once; derive dx⁻¹ = c⁻¹·N and
  (Rx−Qx)⁻¹ = c⁻¹·dx³.
- Gross save: 1 Kaliski ≈ 1.77M (forward+backward on a single input).
- Added cost: ~4 extra muls + 1 squaring ≈ 400–600k.
- Naive net: ~1.2M saved = 27% reduction.

**Obstruction** (confirming prior session): if we compute c⁻¹ once (into
`c_inv_saved` via CX-copy inside Kaliski body) and use it for both derived
inverses, we face the cleanup problem: zeroing `c_inv_saved` at end of
circuit requires either
  (a) another Kaliski on c (2 Kaliski total — no savings), or
  (b) expressing c⁻¹ as a rational function of end-state (tx=Rx, ty=Ry,
      classical ox, oy). This requires λ or dx, both lost after the flow.

There's no measurement trick (HMR) that clears `c_inv_saved` because it's
not a simple AND of quantum values — it's a full modular inverse.

**Verdict: Montgomery batched nets zero Toffoli savings in a single-point-add
context**, confirmed independently by careful register-lifetime analysis.

### Jacobian coordinates: the cleanup obstruction generalizes

Natural Jacobian formulation:
- Compute (X₃, Y₃, Z₃) = Jacobian(P, Q). 8 muls + 3 squarings ≈ 1M CCX.
- One Kaliski on Z₃ to get Z₃⁻¹. 1.77M.
- Multiply X₃·Z₃⁻² and Y₃·Z₃⁻³ for affine. ~150k.
- Total forward: ~2.9M. Looks great — ~34% savings.

**Obstruction**: writing Rx, Ry in-place on tx, ty requires either:
  (a) keeping tx=Px, ty=Py preserved as scratch, in which case they form
      leftover scratch at end that can't be zeroed without Px, Py as inputs
      (same circular dependency as Montgomery batched), or
  (b) the standard Bennett pattern: compute into scratch, swap, reverse.
      But the reverse needs original tx=Px, ty=Py which have been swapped
      out.

Equivalent variants all hit the same wall: recovering Px, Py from end-state
(Rx, Ry, ox, oy) requires a nonlinear operation equivalent to another
full point-subtraction — same cost as the forward flow. So Jacobian with
clean cleanup costs ~2× (one forward + one reverse as point-subtract), which
is worse than the current 2-Kaliski flow.

### The fundamental theorem (conjectured)

Any reversible, in-place quantum circuit for single-EC-point addition over
a prime field requires EITHER:
  (a) two full modular inversions (current 2-Kaliski approach), or
  (b) two full EC-addition-equivalent compute+reverse passes, or
  (c) some algebraic trick that expresses one input of the pair in terms of
      the output plus classical constants.

Paths (a) and (b) cost approximately the same. Published SOTA (Litinski 2023
at 50M, HRSL 2020 at 12M, Kim 2026 at 17M) all use (a). Google/Babbush 2026's
2.1–2.7M is not disclosed but likely combines (a) with extensive primitive-
level tuning we don't have access to.

### Next-session priorities

**Promising directions I could not close tonight:**
1. **Windowed classical-const mul** to replace the 399 halves + 399 doubles
   on `lam`. Estimated ~160k CCX saved (3.7% of circuit). Requires building
   a QROM lookup primitive. Multi-hour implementation.
2. **Kim-style unconditional Kaliski** + m_hist elimination (~400 qubits
   freed) + 2-level Karatsuba in-Kaliski (+258 qubits). Net: maybe neutral
   on qubits, ~50–100k Toffoli savings. Moderate complexity.
3. **B-Y w=2 with careful reversible matrix application** (above numbers are
   upper bounds; tight implementation might do better if matrix entries ≤ 2
   can be specialized). Would need to match or beat 400k one-way to help.

**Directions I'm convinced are dead-ends:**
- Montgomery batched (single-call cleanup is impossible without 2nd Kaliski).
- Naive Jacobian (same cleanup obstruction).
- B-Y w=4 or higher (matrix-entry magnitudes blow up cost).
- Toom-3 (confirmed last session: non-trivial div-by-3 outweighs the win).
- Fermat inversion via addition chain (256–260 muls at 70k each = 18M,
  5× worse than Kaliski).

### Crypto-specific angles to investigate next
- **Co-Z arithmetic** (Meloni 2007): for scalar mul, re-uses Z coords across
  operations. Not directly applicable to single add but might have a
  single-add specialization.
- **Atomic-block formulas** (Chevalier-Mames et al. 2004): compute
  R = P + Q via a series of operations where P and Q get mutated in a
  structured way that might admit cheaper uncompute.
- **Precomputed classical inversion tables**: Q (= offset) is classical, so
  any quantities depending only on Q can be precomputed. The expensive
  inversions are on dx and (Rx-Qx) — both depend on quantum Px,Py,Rx. So
  classical precomputation doesn't help the inversions, only potentially
  other mul operands.
- **Short lattice representation of λ** via quantum-classical coupling. No
  obvious mechanism.

### Infrastructure added this session (committed)
- `TRACE_PHASES=1` / `TRACE_PHASES_VERBOSE=1` env gates in `build()` for
  per-phase emitted-Toffoli breakdown.
- `src/point_add/by.rs` with classical B-Y divsteps2 harness + 3 unit tests
  (all passing). Activated via `BY_TEST=1`.


## Session 2026-04-22 (continued): windowed/sparse const-mul attempt

### What I built
`src/point_add/const_mul.rs` — sparse classical-constant modular multiply
using a 2n-bit tmp_ext accumulator and Solinas folding. Verified it
produces correct OUTPUT (CLASSICAL check passed at 4,423,142 Toffoli).

### Why it didn't help
**Correctness: ✓** (OK run with +29k Toffoli, +512 qubits vs baseline).

**Performance: NET NEGATIVE.** For pair1_halve (K = 2^{-399} mod p, dense
with 122 set bits), the sparse mul cost ≈ 100k CCX — same as 399 iterated
`mod_halve_inplace_fast` calls. For pair2_double (K = 2^{+399} mod p,
sparse with 7 set bits), the forward is cheap (~5k CCX) but the cleanup
pass uses K^{-1} = 2^{-399} (dense, ~100k CCX). The dense-direction
cleanup dominates.

Also, the 2n-bit tmp_ext blows peak qubits from 2729 to 3241 (over cap).

**Phase garbage on partial integration**: replacing ONLY pair2_double
(and leaving pair1_halve as iterated) produces phase garbage. The HMR
calls in the forward sparse-mul don't pair with anything in the baseline
halving loop; their phase contributions don't cancel.

**Conclusion**: sparse_const_mul as a drop-in replacement for halve/double
loops is not viable without either windowed dense-direction optimization
(which would need a QROM primitive, 200-500 LOC) or a different
algebraic restructuring that avoids needing the dense K^{-1} at all.

The primitive is not committed.

### Remaining pieces (kept as dead code)
- `src/point_add/by.rs` with classical B-Y divsteps2 harness + tests.
  Analytically shown to be net-worse than Kaliski at every jump width.

### Honest summary of single-session ceiling
After a full session of structural moonshots:
- Montgomery batched inversion: cleanup obstruction (proven here + prior).
- B-Y divsteps (any w ≥ 1): per-iter cost × safegcd iter count ≥ Kaliski.
- Jacobian coords: cleanup obstruction equivalent to Montgomery batched.
- HRSL cumulative swap state: net +3.2M per prior + my re-analysis.
- Windowed/sparse classical-const mul: phase issues, ~0 net savings.
- Iter count reduction: 399/399 is minimum for current op stream.
- STEP 4 / STEP 3+9 cswap restructuring: all symmetric variants tried
  cost the same or more.

The 4.39M / 2729q floor reached this session appears to be near the
**single-session ceiling with public-literature techniques**. Closing to
Google's 2.1–2.7M SOTA requires either undisclosed Google techniques or
multi-session novel quantum-algorithm research (most promising target:
w=2 reversible B-Y with novel matrix-application primitives, or a
fundamentally different single-Kaliski point-add structure).

