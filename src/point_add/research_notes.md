# Research notes — inversion moonshots inside `src/point_add/`

Session: 2026-04-22 (continued, moonshot-only work).

This file keeps all moonshot literature / classical-analysis work under
`src/point_add/`, per the current scope rules.

## Deliverable 1 (classical B-Y on secp256k1) — confirmed

Implemented classical `divstep2` reference and modular-inverse recovery in
`src/point_add/by.rs`, then ran a 10,000-input secp256k1 survey.

Results:

| metric | value |
|---|---|
| theoretical bound `⌈(49·256 + 57)/17⌉` | 742 |
| observed minimum iters | 502 |
| observed maximum iters | 567 |
| observed mean iters | 531.01 |
| max `|δ|` observed | 20 |
| modinv matches (vs Fermat) | 10,000 / 10,000 |

Interpretation:
- The BY safegcd upper bound is pessimistic by ~24% on secp256k1 inputs.
- However, this is **not enough** to save plain B-Y: the per-iter reversible
  cost is still too high relative to Kaliski.

## Deliverable 2 (algorithm-space survey) — corrected final version

### 1. Kaliski almost-inverse (baseline)
- Classical ref: Burton S. Kaliski Jr., “The Montgomery inverse and its
  applications,” IEEE Trans. Computers 44(8), 1995.
- Quantum / reversible refs:
  - Roetteler–Naehrig–Svore–Lauter 2017, arXiv:1706.06752.
  - Häner–Roetteler–Soeken 2020, arXiv:2001.09580 / ePrint 2020/077.
- Iterations in our tuned circuit: 399.
- Measured per-iter reversible cost: ~2180 CCX.
- Per-pass cost: ~1.81M CCX.

### 2. Bernstein–Yang divstep2 (w = 1)
- Ref: Bernstein–Yang 2019, ePrint 2019/266.
- Reversible implementation: unpublished / would be novel.
- Empirical iterations on secp256k1: max 567, mean 531.
- Per-iter reversible estimate: 10–12n CCX.
- Conclusion: still worse than Kaliski.

### 3. Bernstein–Yang jumpdivsteps2 (w > 1)
- Ref: Bernstein–Yang 2019, Figure 10.2 / §10.
- Reversible implementation: unpublished / would be novel.

#### 3a. Corrected matrix-growth result
A previous version of the jump survey undercounted the scaled transition
matrix. After fixing it, the 100,000-sample survey now shows the **full
scaled** transition matrices do hit the theoretical `2^w` growth.

Corrected survey over 100,000 random low-word states:

| w | max observed `|entry|` | max log2 | mean log2 | theoretical max log2 |
|---|---:|---:|---:|---:|
| 4  | 16    | 4.00  | 2.03 | 4  |
| 8  | 256   | 8.00  | 4.28 | 8  |
| 12 | 4096  | 12.00 | 6.34 | 12 |
| 16 | 65536 | 16.00 | 8.19 | 16 |

Interpretation:
- The **maximum** entry size really does hit the full `2^w` growth.
- So a faithful reversible matrix-apply must still handle `w`-bit classical
  coefficients.
- That restores the pessimistic reversible cost model: batching by `w` does
  not automatically beat Kaliski.

#### 3b. Exact matrix-family compression result
Even if entries hit `2^w`, a quantum QROM implementation might still benefit
if the number of **distinct** transition matrices is tiny compared to the raw
state space. I measured this exactly for all low-word states with
`delta ∈ [-20, 20]`, odd `f_low`, and arbitrary `g_low`.

Results:

| w | total states | distinct matrices | compression factor |
|---|---:|---:|---:|
| 4 | 5,248 | 656 | 8× |
| 6 | 83,968 | 2,624 | 32× |
| 8 | 1,343,488 | 10,496 | 128× |

Pattern:
- compression factor = `2^(w−1)` exactly on the observed range.
- equivalently, distinct matrix count appears to scale like `2^(w+2)`.

This does **not** rescue full jumped B-Y by itself, but it is a strong sign
that *compressed local transition classes* are real and exploitable.

#### 3c. Updated verdict on jumped B-Y
Full jumped B-Y still looks too expensive as a drop-in replacement, because:
- matrix entries hit the full `2^w` growth,
- full coefficient tracking would still need to carry those `w`-bit entries,
- cleanup is all-new machinery.

But the compression result changes the local-batching story.

### 4. Montgomery inverse (Savaş–Koç)
- Classical ref: Savaş–Koç 2000, “The Montgomery modular inverse revisited.”
- Quantum / reversible refs: effectively same family as RNSL/HRSL Kaliski.
- Conclusion: not a distinct win over Kaliski in our setting.

### 5. Lehmer-style GCDs
- Classical refs: Lehmer 1938; Jebelean 1993.
- Reversible implementation: unpublished / novel.
- Main issue: runtime matrix selection depends on quantum data, so a faithful
  reversible implementation needs a QROM keyed by top bits. No concrete,
  literature-backed reversible cost win established yet.
- Still potentially interesting as novel research, but now less grounded than
  a compressed Kaliski-local batching route, because we have exact empirical
  class-compression data for the latter.

### 6. Fermat / addition-chain inversion
- Standard classical method; discussed in cryptographic resource estimates.
- Prime-field reversible cost is far too large (hundreds of multiplications).
- Not competitive.

### 7. Itoh–Tsujii
- Only for GF(2^n), not GF(p).
- Not applicable to secp256k1.

## Stronger result: coefficient-side compression matches (u, v) compression

A remaining risk in the hybrid Kaliski-jump idea was that even if the `(u, v)`
window transition family compressed well, the coefficient-side `(r, s)`
transforms might explode and ruin the QROM story.

I derived the per-case coefficient matrices directly from the implemented
`kaliski_iteration` logic:

- UEven: `(r, s) -> (r, 2s)`
- VEven: `(r, s) -> (2r, s)`
- UGtV : `(r, s) -> (r+s, 2s)`
- VGtU : `(r, s) -> (2r, r+s)`

Then I ran the same exact 10,000-input window survey for those coefficient-side
matrices.

**Result:** the `(r, s)` side compresses **identically** to the `(u, v)` side.

| w | t | distinct uv mats | distinct rs mats | max `|uv|` | max `|rs|` | mean mats/class |
|---|---:|---:|---:|---:|---:|---:|
| 6 | 4 | 125  | 125  | 16 | 16 | 4.506 |
| 8 | 4 | 125  | 125  | 16 | 16 | 4.493 |
| 8 | 6 | 1133 | 1133 | 64 | 64 | 9.461 |

This is the strongest empirical evidence so far that **hybrid Kaliski-jump**
is a coherent moonshot and not just a half-broken idea.

## Current best moonshot conclusion

**Conclusion: `hybrid Kaliski-jump is the bet.`**

This is now stronger than the previous statement.

### Why full B-Y replacement is not the best bet
Full BY jumpdivsteps2 still has two major problems:
1. matrix entries hit the full `2^w` growth;
2. coefficient tracking and cleanup are all-new machinery.

So a *full* B-Y replacement remains very high-risk.

### Why the histogram result matters
The exact histogram shows there are vastly fewer distinct local transition
matrices than raw low-word states. That suggests a more focused route:

> keep Kaliski's global state machine and cleanup structure,
> but replace short local runs of the `(u, v_w)` update path with
> **compressed pre-batched transition classes**.

This attacks the actual hot path while preserving the machinery that we already
know is reversible and correct.

## New classical proposal: hybrid Kaliski-jump

### Model
Standard Kaliski / binary almost-inverse update on `(u, v)` has four branch
cases:

```text
if u even:                   (u, v) ← (u/2, v)
elif v even:                 (u, v) ← (u, v/2)
elif u > v:                  (u, v) ← ((u-v)/2, v)
else:                        (u, v) ← (u, (v-u)/2)
```

Each step is a linear map with a shared `1/2` factor. Over `t` steps we get
an integer 2×2 matrix `P_t` with

```text
(u_t, v_t)^T = (1 / 2^t) · P_t · (u_0, v_0)^T.
```

The classical question is: along actual secp256k1 trajectories, keyed by low
`w` bits of `(u, v)`, how many distinct `P_t` arise? If small, we can imagine
QROM-selecting those classes instead of executing the per-step parity / compare /
cswap / sub / halve sequence.

### Empirical hybrid Kaliski-window survey
I added `src/point_add/kaliski_jump.rs` and sampled actual Kaliski trajectories
for 10,000 random secp256k1 inputs. Windows overlap (advance one step, observe
`t`-step lookahead), because that's the runtime use-case.

The survey now tracks **three** transition alphabets per `(w, t)`:
1. state-side `(u, v_w)` matrices,
2. coefficient-side `(r, s)` matrices,
3. the **joint pair** `(uv_mat, rs_mat)` — i.e. the actual object a reversible
   batched primitive would need to look up.

Results:

| w | t | distinct `uv` mats | distinct `rs` mats | distinct joint pairs | max `|entry|` | mean mats/class | max mats/class |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | 4 | 125 | 125 | 125 | 16 | 4.506 | 16 |
| 8 | 4 | 125 | 125 | 125 | 16 | 4.493 | 16 |
| 8 | 6 | 1,133 | 1,133 | 1,133 | 64 | 9.461 | 62 |

Interpretation:
- For **t = 4**, the **entire joint transition family** is only **125** matrices.
- For **t = 6**, the joint family is still only **1,133** matrices.
- The coefficient-side `(r, s)` compression is not merely similar — it is
  effectively identical to `(u, v_w)` on the sampled trajectories.
- Most importantly, the **joint pair count equals the per-side count** in the
  sampled data. That means the state-side and coefficient-side transforms appear
  to be locked together, so a single compressed lookup may suffice.

This is now the strongest empirical structural lead in the project.

## Proposed next sessions

### P1. Enumerate exact branch-class representatives for `t = 4`
For the 125 observed 4-step matrices, enumerate:
- a canonical representative branch sequence,
- exact parity/ordering preconditions,
- exact `(u, v)` low-bit regions mapping to each matrix.

This is the step needed before any reversible QROM design.

### P2. Build a reversible cost model for a 125-matrix QROM
Now that the matrix alphabet is only ~125 elements for `t=4`, the next work is
not abstract algorithmics but concrete reversible cost accounting:
- raw lookup,
- compressed-class lookup,
- select-swap QROM,
- matrix-apply on 256-bit regs,
- cleanup interaction with existing `m_hist`.

### P3. Decide whether `t=4` or `t=6` is the sweet spot
`t=4` gives only 125 matrices with max coefficient 16.
`t=6` gives 1,133 matrices with max coefficient 64.
Need to compare fewer batches vs. larger matrix-apply cost.

## Bottom line

The strongest current research judgement is:

> The best moonshot is **not** full B-Y replacement.
> The best moonshot is **hybrid Kaliski-jump batching** over short windows,
> because the exact local transition family is very small on both the state
> side `(u, v_w)` and the coefficient side `(r, s)`.

That's still novel research, but unlike the other moonshots, it now has
clear empirical support directly tied to the 81%-of-budget hot path.
