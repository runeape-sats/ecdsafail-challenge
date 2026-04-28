# Structural SOTA plan — 2026-04-28

User directive: stop local tuning. This file is the working structural model for
matching Google's secp256k1 point-add frontier, not a micro-optimization list.

## 0. Current measured point

Current committed baseline after the last in-flight retune:

- **4,132,750 Toffoli**
- **2716 qubits**
- exact / phase-clean under the 24-seed gate and checks

Google/ZKP targets:

- low-qubit: **2.7M Toffoli @ 1175 qubits**
- low-gate: **2.1M Toffoli @ 1425 qubits**

So the real gaps are:

- **−1.43M Toffoli** to low-qubit
- **−2.03M Toffoli** to low-gate
- **−1541 qubits** to low-qubit

## 1. The Toffoli gap is one inversion-sized object

Measured decomposition of the current design:

| component | cost |
|---|---:|
| Kaliski invocation #1 (`with_kal_inv_raw`, fwd+body+bwd scale excluded) | ~1.60M |
| Kaliski invocation #2 | ~1.59M |
| non-Kaliski scaffold (muls, scale correction, Solinas, constants) | ~0.94M |
| total | ~4.13M |

Therefore a SOTA-grade design must do one of exactly two things:

1. **Delete one full inversion invocation.** Current primitives then land at
   roughly `4.13M - 1.60M = 2.53M`, already matching Google's low-qubit
   Toffoli target.
2. **Keep two inversions but make each ~45% cheaper.** Need per-invocation
   cost `~1.60M -> ~0.90M`, i.e. save ~700k per Kaliski invocation.

Anything that does not attack one of these two objects cannot match SOTA.
A 50k improvement is useful only if it is a stepping stone toward one of
these two structural changes.

## 2. Why the one-inversion route keeps failing

The map is in-place:

```text
(Px, Py) -> (Rx, Ry) = (Px, Py) + (Qx, Qy)
```

A one-inversion affine formula exists classically (Strategy C, `w = dx^3`),
but reversible cleanup is the hard part. The issue is not the algebraic
formula for `Rx,Ry`; it is zeroing the old input information.

The obstruction in every attempted one-inversion schedule is:

- To uncompute the inversion input (`dx`, `dx^3`, or a product containing it),
  the circuit needs a live source of `dx` after `tx` has become `Rx`.
- Reconstructing `dx = Px-Qx` from `(Rx,Ry,Qx,Qy)` is exactly point subtraction
  by the classical point `Q`.
- Point subtraction has the same denominator `(Rx-Qx)` that the current second
  Kaliski invocation inverts.

So the naive one-inversion design just moves the second inversion into the
cleanup path.

### Strategy C re-estimate at the current baseline

Classically correct formula:

```text
dx = Px-Qx, dy = Py-Qy, w=dx^3
v  = dy^2 - dx^2(Px+Qx)
Rx = v * dx * w^-1
Ry = (dy(dx^2 Qx - v) - w Qy) * w^-1
```

Cost at current 407/403 Kaliski settings:

| block | estimated Toffoli |
|---|---:|
| compute `dx^2, w=dx^3` | 250-300k |
| one Kaliski invocation on `w` | 1.60M |
| compute `dy^2, v, Rx, Ry` with Bennett-clean temps | 1.7-2.0M |
| uncompute `w, dx^2` | 250-300k |
| misc/scale | 100-200k |
| **total** | **3.9-4.4M** |

This is not a SOTA route unless the cleanup can be made triangular/in-place
instead of Bennett-clean. That triangular schedule is the unresolved research
problem.

## 3. Therefore the most credible Toffoli route is a cheaper Kaliski

Current Kaliski cost sources, across both forward and backward and both
invocations:

| substructure | total cost | SOTA relevance |
|---|---:|---|
| step 3 + step 9 cswaps | ~1.0M+ | biggest public lever |
| step 4 controlled `v-=u; s+=r` | ~1.0M | second biggest lever |
| comparator / eqzero / flags | ~0.4M | moderate |
| scale correction loops | ~0.2M | not enough alone |

To reach 2.7M with two inversions, we need to remove/replace **most of the
cswap + step4 cost**. This points to a jumped/windowed/divstep-style Kaliski,
not to local adder swaps.

## 4. Candidate structural programs

### Program A — jumped/windowed Kaliski (highest Toffoli relevance)

Batch `t` binary-GCD microsteps into one matrix update:

```text
[u']   1/2^t [a b] [u]
[v'] =       [c d] [v]

[r']          [A B] [r]
[s'] =        [C D] [s]
```

For `t=2..4`, coefficients are small. The hope is that one matrix-selected
update costs less than `t` copies of:

- two full cswap layers, and
- one step4 sub/add layer.

Target economics:

| window | current cost for t steps | needed jumped cost | result |
|---:|---:|---:|---|
| t=2 | ~2× current step | <1.2× current step | ~40% Kaliski win |
| t=3 | ~3× | <1.8× | ~40% win |
| t=4 | ~4× | <2.4× | ~40% win |

This is exactly the size of win required for SOTA with two inversions.

**Fast invalidation criterion:** if reversible matrix application needs more
than ~2 q-q adds/subs plus one controlled shift per microstep, it cannot beat
current Kaliski. If a low-coefficient t=2 or t=3 schedule can be synthesized
with one scratch n-register and <=~1.5n Toffoli/step equivalent, it is live.

Fresh survey from `kaliski_jump.rs` / scratch replay (10k inputs for the
low-bit table, 2k inputs for quick t sweep):

| t | distinct matrices | max coeff | mean log2 coeff | mean mixed rows |
|---:|---:|---:|---:|---:|
| 2 | 13 | 4 | 1.79 | 0.85 |
| 3 | 41 | 8 | 2.58 | 1.12 |
| 4 | 125 | 16 | 3.29 | 1.34 |
| 6 | 1133 | 64 | 4.71 | 1.63 |

Low-bit lookup is **not** enough: at `w=8,t=4`, each low-bit class still has
`mean 4.49` possible matrices and up to `16`; at `w=8,t=6`, mean `9.46`, max
`62`. A strengthened executable invalidation
`initial_gt_window_classifier_not_approx_good_enough` adds one full comparator
bit (`u>v`) to the low-bit key and still sees a disjoint-sample majority error
of about **60%** for `w=8,t=4`. So a one-comparator window is not even close to
1% approximate correctness.

A real jumped Kaliski must either compute the whole comparator sequence
coherently, use a Bernstein-Yang/divstep variable that avoids full comparisons,
or synthesize a matrix application whose cost beats the per-step loop despite
those predicates.

Positive qubit-side result: `window_hint_bits_can_compress_history_but_not_select_matrix_alone`
records the actual matrix choice as a small per-window hint instead of per-step
history. On 5k sampled trajectories with key `(low8(u), low8(v), u>v)`:

| t | max matrices/key | hint bits/window | total hint bits |
|---:|---:|---:|---:|
| 4  | 8  | 3 | 306 |
| 8  | 23 | 5 | 255 |
| 16 | 34 | 6 | 156 |

So window hints can plausibly save 100-250 history qubits versus `m_hist`, but
only if a selected matrix can be applied cheaper than replaying microsteps. This
is qubit-structural, not yet a Toffoli route.

Next concrete work: extend `kaliski_jump.rs` from matrix observation into a
matrix-application lower-bound/synthesis search for `t=2,3`, with comparator
cost included. If t=2 cannot beat two ordinary iterations after paying for the
second-step branch computation, kill this path quickly and move to BY/divstep.

### Program B — triangular one-inversion schedule (highest payoff, highest risk)

Goal: use Strategy C or B2 but avoid Bennett-clean fresh outputs. A successful
schedule must satisfy:

1. Kaliski input (`dx` or `dx^3`) remains uncomputable after output mutation.
2. `tx,ty` are transformed in-place, not via fresh `(rx,ry)` registers.
3. Any copied slope/inverse is phase-uncomputed from live state without
   inverting `(Rx-Qx)`.

Fast invalidation criterion: if the schedule ever contains both
`old dx` and `new Rx` as independent live n-bit values after Kaliski backward,
it has already lost; zeroing one from the other is point-subtraction.

### Program C — Kaliski representation rewrite for qubits

Even a cheaper Kaliski will not hit 1175q with the current state layout.
Required qubit reductions:

| source | potential saving |
|---|---:|
| `m_hist` compression/elimination | 407q |
| fold Kaliski input copy into `tx`/scratch | 256q |
| fold `r` output into next multiplier/output scratch | 256q |
| low-workspace step4 / venting only at peak sites | 200-260q |
| register sharing / length-tracked tail | 100-300q |

The easy `m_i` start-state formula is **not sufficient** because iter-end
state does not expose that fingerprint. A real qubit breakthrough needs either:

- a self-cleaning Kaliski body whose inverse branch is recoverable from end
  state, or
- pebbling/checkpoint recomputation, or
- Luo-style length/location registers.

This is qubit-structural, but by itself does not solve Toffoli.

### Program D — coset/padded arithmetic only after a long region exists

Exact padded/coset add chains cross over at ~12 repeated additions and save
~44% by 256 additions, but cost +500-800 qubits in the current prototype.
Current affine scaffold has no long non-Kaliski add region. Coset becomes
relevant only if Program A creates windowed/batched arithmetic regions or if
we accept a larger representation rewrite.

## 5. Decision rule going forward

Do not pursue a code change unless it plausibly satisfies at least one:

- **Toffoli structural:** can save >=0.5M by deleting an inversion or cutting
  Kaliski per-step cost.
- **Qubit structural:** can save >=400q without >0.5M Toffoli regression.
- **Fast falsification:** conclusively kills a tempting structural path so we
  do not waste another session on it.

Immediate next target: **Program A**, because it is the only public-ish path
whose economics can plausibly produce the missing 1.4M Toffoli while staying
inside the exact harness.
