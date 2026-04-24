# SOTA rebuild plan

North Star: **2.1M–2.7M Toffoli, 1175–1425 qubits, single secp256k1 point-add,
exact affine `(Rx, Ry) <- (Px + Qx, Py + Qy)`, `(ox, oy)` classical bit input.**

Current live: **4,180,502 Toffoli, 2716 qubits**.

Gap: -35% Toffoli, -57% qubits. Not reachable by micro-optimizations.

## What we already validated classically (kept artifacts)
- `single_inv_numeric.rs`: single-inversion point-add formula works on
  200/200 secp256k1 trials. Three concrete variants.
- `kim_proto.rs`: Kim-style wide-r unconditional 2n-round Kaliski, followed
  by a single classical `× 2^{-2n}` unscale, produces the true modular
  inverse on 200/200. End-to-end `kim_style_end_to_end_point_add_passes_200_trials`
  confirms the full point-add works under this scheme.
- `luo_proto.rs`: Luo register sharing saves 632–888 qubits on the
  inversion state alone, but inversion-savings-alone ceiling is ~1828q —
  still misses the 1175q target by ~716q.
- `lit_tricks.md`: SOTA chassis is still exact six-step affine add (Litinski
  2023, Google 2026, validated by arXiv 2506.03318). Gap to SOTA is in the
  subroutine implementations, not in the top-level formula.

## Architecture we need to build
A new inversion module that at its core is:
1. **Unconditional**, `2n` fixed rounds — no termination flag, no `f` qubit,
   no `m_hist`.
2. **Wide `r,s`** registers — `r` is `(2n+1)`-bit, `s` is `(2n+1)`-bit.
   Postpone mod-p reduction to a single pass at the end.
3. **Register-shared** — `r_i` and `t_i` cohabit the same `n+2` qubit
   register under the bit-length invariant `bitlen(r_{i-1}) + bitlen(t_i) <= n+1`.
   (Luo's sharpest form.)
4. **Shared workspace with mul** — after the inversion completes, the
   inversion's auxiliary registers (now known to hold `{0, 1, p}`) get
   cleared with n X gates and immediately reused as multiplication tmp
   register. (HRSL Fig. 8b.)
5. **Single inversion per point-add** (not two). Write the scaffold as
   exact affine one-inv chain with reversible Bennett-style uncompute.

## Decomposition

```
pub fn build() -> Vec<Op> {
    // Inputs as today: tx (2n qubits), ty (2n qubits), ox (n bits), oy (n bits).
    //
    // 1. tx -= ox    (inplace mod_sub_qb)           -> tx = dx
    // 2. ty -= oy    (inplace mod_sub_qb)           -> ty = dy
    //
    // 3. Invert dx into a fresh output register `inv_dx = n qubits` using
    //    Kim/Luo unified inversion (see below).
    //    Body:
    //      a) copy tx into v_w (n+2-bit shared reg) via CX
    //      b) initialize u := p (classical bits flipped into u register)
    //      c) 2n fixed unconditional Kaliski rounds (no f, no m_hist)
    //      d) at termination r holds raw = ±dx^{-1} * 2^{2n} in wide form
    //      e) one final classical `× 2^{-2n} mod p` on r (windowed const-mul)
    //      f) copy r into `inv_dx` via CX
    //      g) emit_inverse of the whole inversion block
    //         -> u, v_w, r all return to 0, tx untouched, inv_dx holds 1/dx
    //
    // 4. Compute lam into a fresh register: mod_mul_write_into_zero_acc(lam, ty, inv_dx)
    //
    // 5. Exact affine body:
    //      tx := dx - lam^2
    //      tx += 3*Qx
    //      tx := -tx                            (== λ² - dx - 3Qx)
    //      tx += 2*Qx                           (== Rx - Qx)
    //      ty += lam * tx                       (ty = dy + λ·(Rx-Qx) = dy - (Ry+Qy))
    //      ...and the exact ty correction per Litinski step 2/5/6.
    //
    // 6. Uncompute lam = dy * inv_dx via emit_inverse on the mul. 0-qubit cost.
    //
    // 7. Uncompute inv_dx by re-running the inversion + inverse of step 3f.
    //    Bennett style: two inversion calls total, same as HRSL.
    //
    // 8. Final tx += Qx to land Rx.
}
```

## Operation budget (target, decomposed)

Target: 2.7M Toffoli total for 1175q variant.

Allowed per subroutine, roughly:
- inversion (x 2, for Bennett uncompute): ~2 * 900k = 1.8M Toffoli
- λ mul (dy * inv_dx) + uncompute: ~2 * 150k = 300k
- λ² squaring + uncompute: ~2 * 130k = 260k
- λ·(Rx-Qx) mul + uncompute: ~2 * 150k = 300k
- other (constant adds, negations, flag mgmt, fix-up): ~40k
- total: ~2.7M ✓

Our current inversion is ~2M per pass (and we do 2 passes). Target is
900k per pass via Kim unconditional + Luo register sharing.

## Per-round Kaliski budget (unconditional, wide-r, shared)
- cond swap on `(u, v_w)` width n+2: 3(n+2) CCX
- cond swap on `(r, s)` width n+2: 3(n+2) CCX
- `v := v - u` width n+2 (Cuccaro): 2n CCX
- `s := s + r` width n+2 (Cuccaro): 2n CCX
- `v := v / 2` shift: 0 CCX
- `r := 2r` shift: 0 CCX
- branch flag bit computation + uncompute (m_i on fly): ~4 CCX

Per round ~12n CCX. 2n rounds = 2n * 12n = 24n² CCX.
At n=256: 24 * 65536 = 1,572,864 Toffoli per inversion. 2× Bennett = 3.15M.
Over budget by ~450k.

Need to shave per-round to ~7n CCX. That matches HRSL's reported
"swap-based" Kaliski where one sub + one add per round, plus the
cswaps. 7n × 2n = 14n² = ~918k per inversion, 2× = 1.84M. Budget fits.

## Phase A/B/C status (SESSION 2026-04-24)

**Phase A COMPLETE**: reversible Kim iteration with matching backward.
Round-trip clean at n=256 over full 2n=512 iters.

**Phase B COMPLETE**: init + 2n forward produces correct ±x^-1*2^{2n}.

**Phase C COMPLETE**: full reversible `kim_inv(x, out)` with all state
cleanup, verified 3/3 trials via Simulator.
- Toffoli cost: 2,530,240
- Peak qubits (standalone): 4102
- `x` unchanged, `out` holds ± x^-1 * 2^{2n} mod p, no ancilla leaks.

**Phase D NEXT**: wire `kim_inv` into a single-inversion `build_sota()`.
Estimated budget:
- 1 kim_inv call: 2.53M
- λ via dy * inv_dx (+ uncompute): ~300k
- λ² squaring (+ uncompute): ~260k
- λ * (Rx-Qx) mul (+ uncompute): ~300k
- classical 2^{-2n} unscale (if not absorbed): ~200k
- exact affine ops: ~40k
- **projected total: ~3.6M Toffoli** (vs current 4.18M, −14%)

This beats the baseline but isn't SOTA alone. Additional levers needed
to reach 2.7M:
- Luo register sharing on Kim state: qubit reduction, ~small Toffoli cost.
- HRSL workspace reuse: -258 peak qubits + unlocks karatsuba at mul sites.
- Absorb 2^{-2n} unscale into the first mul's Solinas to save ~200k.

Combined projection: 2.9-3.1M, within striking distance of SOTA (2.7M).

## Ordering plan

Phase A: **Build and classically verify the new inversion primitive**
at the reversible-builder level (not yet wired into `build()`).
  - target: `src/point_add/kim_inv_circuit.rs`
  - a dry-run harness: build the circuit, run it through `Simulator` on
    random inputs, check that the n-bit output register equals `x^{-1} mod p`.
  - this bypasses the live build entirely while we work.

Phase B: **Add register sharing (Luo).** Model `r,t` in a shared register.
  - target: in-place within the above module.
  - classical harness verifies output unchanged.

Phase C: **Wire into `build()` behind `KIM_INV=1` env gate.**
  - keep the live 2-Kaliski scaffold as default.
  - measure. If Toffoli drops, keep; remove gate, make default.

Phase D: **Switch to single inversion per point-add** in `build()`.
  - affine 6-step chain from the literature.
  - classical replay already done in single_inv_numeric Strategy C.

Phase E: **Workspace reuse** (HRSL Fig 8b). Share freed inversion-internal
registers as mul tmp during the λ multiplication.

Phase F: **Windowed classical-quantum adds** for the constant ops (Qx, Qy).
Low-impact, late cleanup.

## First concrete deliverable

`src/point_add/kim_inv_circuit.rs` implementing `kim_invert(b, x, inv_out)`:
- produces `inv_out <- x^{-1} mod p`, reversibly
- `x` unchanged
- all internal state returns to |0⟩
- classically tested on random secp256k1 inputs

Not wired into `build()` yet. No primary metric change at this step.

## Hard structural finding from Phase A analysis

Narrow-r-per-round Kim (= the existing `kaliski_iteration_bulk_prefix3`
pattern) cannot extend to the full 2n=512 rounds without re-introducing
the v_w==0 detection (step 0 + f_flag + m_hist), because the per-round
`mod_double_inplace_fast` on a narrow `r` bakes in termination-dependent
scaling. Empirically the safe upper bound of bulk_prefix3 is ~315, and
past that the scaffold breaks because v_w hits 0 inside the bulk region.

Therefore Kim's true win REQUIRES wide r, s (the 2n+1-bit accumulator)
and postponed modular reduction. With wide r:
- once v_w reaches 0, the wide-r shift-left just appends zeros,
- no termination tracking needed,
- no f_flag, no m_hist, no step 0,
- single mod-p reduction pass at the end.

Consequence: Phase A is NOT a small modification of the existing
`kaliski_iteration`. It is a new inversion primitive with different
register widths (r: 2n+1 bits, s: 2n+1 bits). The per-round arithmetic
changes from narrow-mod-p-per-round to wide-shift-and-add with no
modular reduction inside the loop.

Literature-estimated per-round cost (HRSL Fig 7b, Kim Alg 2):
- 1 cswap on wide (u,v) pair: 3(n+2) CCX
- 1 cswap on wide (r,s) pair: 3(n+2) CCX
- `v := v - u` wide: 2(n+2) CCX
- `s := s + r` wide: 2(n+2) CCX
- 1-bit flag generation + uncompute: ~4 CCX
Total ~10(n+2) ≈ 10n per round, 2n rounds ≈ 20n² ≈ 1.31M per inversion.

Two inversions (Bennett) = 2.62M + ~0.3M for muls + ~0.1M misc ≈ 3.0M.

With single-inversion scaffold (Strategy C) + Luo register sharing:
close to Google's 2.7M at ~1800q. Still not 1175q but within range.
