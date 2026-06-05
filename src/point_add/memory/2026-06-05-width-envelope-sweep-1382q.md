# 2026-06-05 Width-envelope sweep + first peak-qubit drop to 1382 (4 promotions)

Session arc on the COMPARE_BITS=52 lineage (alexander-sei's 488afae opened it).
The frontier was very hot early (robertkodra/factory-sagar/alexander-sei grinding
single notches every ~10-30 min), then cooled mid-session. Four promotions landed:

| submission | change | score |
|---|---|---|
| 6936bcd | re-stack `KAL_FOLD 24->23` on the compare52 base (it had reverted) | 2,111,711,630 |
| 11d1423 | `ACTIVE_ITERATIONS 259->258` + `WIDTH_SLOPE 1004->1005` (shared island) | **2,092,653,422** |
| 79d7ee1 | `WIDTH_SLOPE 1005->1006` | 2,091,901,614 |

(11d1423 also via factory-sagar's intermediate 78696aa = `KAL_DOUBLE 24->23`.)

## The two big lessons

1. **Dropping a whole GCD iteration also drops peak qubits.**
   `DIALOG_GCD_ACTIVE_ITERATIONS 259 -> 258` cut ~3,446 executed Toffoli AND moved
   the peak width **1390 -> 1382** (the dropped iteration's scratch row is no
   longer live). First peak reduction in this lineage. This is worth far more than
   a carry-window bit and re-opens the width x peak product. 258->257 on the new
   base is sparse (13+5) — past convergence for many inputs — but a margin/iter
   combo that sheds another scratch row is the highest-value next target.

2. **Break sets can CANCEL across levers -> stack two savings under one island.**
   On the 1382 base, `slope1005` alone = 4+3, `active258` alone = 4+3, but
   **`slope1005 + active258` = 3+3** (gentler than either). Always measure the
   combo break count, not just the sum. (Did NOT recur on the next base:
   slope1006+active257 = 26+18, purely additive. Lever-set dependent.)

## Density map (classical+phase breaks at inherited nonce, by base)

- compare52 base (488afae): kal_f23 = 0+2 (island nonce 55, ~15s).
- 1382 base (11d1423): slope1006 = 6+6 (island nonce 13555, ~30min/13.6k nonces),
  active257 = 13+5, margin9 = 9+4, applyc19 = 10+6. No cancelling combos.

Phase breaks hurt density more than the classical count alone suggests
(6+6 took ~13.6k nonces vs 3+3 at ~3.6k and 0+2 at ~55).

## The wall: WIDTH_MARGIN 10 -> 9

The biggest remaining single lever (~-4,184 T, likely another peak drop) but a
sparse island: 9+ breaks, >1/14000 by local brute search (killed at 14k nonces).
Brute tail-nonce search at ~7-8 nonce/s on 11 threads runs out of runway here.

**To break this wall, build the both-factor classical convergence pre-filter**
(alexander-sei's 488afae note + robertkodra clearly use one). Recipe:
- For each test input derive the TWO GCD inputs per point-add: `dx = Px - Qx mod p`
  and `c = Qx - Rx mod p` where `Rx = E.x`, `E = target + offset` (expected point,
  so Rx is free from the precomputed expected point). Checking only `dx`
  undercounts the hard rate ~2x and yields false-clean islands.
- Run the truncated-width binary GCD classically (model `active_width(step)` from
  the WIDTH_MARGIN/SLOPE/ACTIVE envelope + body-carry band) and reject any input
  whose register width overflows the envelope (= a value error = hard input).
- A nonce is clean iff NONE of its 9024 inputs is hard on EITHER factor.
- Must be bit-exact with the quantum truncation; validate by reproducing known
  islands (e.g. nonce 13555 clean, neighbors dirty) before trusting it.
`src/point_add/kaliski_classical_replay.rs` has the FULL-WIDTH recurrence as a
starting point, but it does NOT yet model the width envelope — that's the work.

## Tooling

`src/bin/island_search_jac.rs` — full-validation tail-nonce searcher, bit-exact
with eval_circuit (startup self-check vs affine k*G). ~7-8 nonce/s on 11 threads,
sim-bound. Usage:
`ISLAND_LOWQ=1 <KNOB=val ...> ISLAND_THREADS=11 ./target/release/island_search_jac <start> <count> [step]`
Reports `CLEAN nonce=N avg_toffoli=T ...` on the first clean island. Note the
avg_toffoli is exact but peak qubits must be read from the official `ecdsafail run`
(the searcher runs LOWQ) — that's how the 1382 drop was first observed.
