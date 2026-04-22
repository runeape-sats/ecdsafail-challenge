# Research notes — inversion moonshots inside `src/point_add/`

Session: 2026-04-22 (continued, moonshot-only work).

This file replaces the deleted top-level `research_notes.md` and keeps all
moonshot literature / classical-analysis work under `src/point_add/`, per the
current scope rules.

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

**Correction to the previous session's optimism:**
I fixed the jump-matrix survey code. The earlier result claiming entries were
much smaller than `2^w` was wrong because the matrix accumulation was wrong.

Corrected survey over 100,000 random low-word states:

| w | max observed `|entry|` | max log2 | mean log2 | theoretical max log2 |
|---|---:|---:|---:|---:|
| 4  | 16    | 4.00  | 2.03 | 4  |
| 8  | 256   | 8.00  | 4.28 | 8  |
| 12 | 4096  | 12.00 | 6.34 | 12 |
| 16 | 65536 | 16.00 | 8.19 | 16 |

Interpretation:
- The **maximum** entry size really does hit the full `2^w` growth.
- That means a faithful reversible implementation must still budget `w`-bit
  classical coefficients for matrix-apply.
- So the reversible matrix-apply cost scales like `w · n`, which cancels the
  `1/w` reduction in iteration count.
- **Conclusion: jumped B-Y does not beat Kaliski.**

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
- Still potentially interesting as novel research, but much less grounded than
  the previously hoped-for jumped B-Y path.

### 6. Fermat / addition-chain inversion
- Standard classical method; discussed in cryptographic resource estimates.
- Prime-field reversible cost is far too large (hundreds of multiplications).
- Not competitive.

### 7. Itoh–Tsujii
- Only for GF(2^n), not GF(p).
- Not applicable to secp256k1.

## Deliverable 3 — final conclusion

**Conclusion: `no known algorithmic path remains open — the gap to Google's
SOTA likely requires undisclosed techniques.`**

This is the corrected answer after fixing the jumpdivstep survey.

Why:
1. Plain B-Y is still too expensive even with empirical 567-iteration max.
2. Jumped B-Y loses because transition matrices really do hit `2^w` growth,
   so matrix-apply scales as `w · n` and cancels the batching gain.
3. Montgomery batched inversion remains blocked by the cleanup obstruction:
   zeroing the saved inverse state needs another inversion.
4. Jacobian / projective approaches hit the same cleanup obstruction in the
   single-point-add reversible setting.

So the remaining routes are not “known better algorithms,” but rather:
- a fundamentally new reversible inversion algorithm,
- a hybrid Kaliski batching method not in the literature,
- or hidden/undisclosed engineering from the Google result.

## Proposals for future sessions (still moonshots)

### Proposal P1: hybrid Kaliski jump-batching
Batch small windows of the Kaliski state machine itself, instead of replacing
it wholesale with B-Y. This avoids the full matrix-growth problem and preserves
our existing cleanup structure.

### Proposal P2: novel reversible Lehmer
Still possible, but it is much less justified by current evidence than it
looked when the jumpdivstep survey was wrong.

### Proposal P3: QROM-backed classical-constant multiplication
Still worth doing for a few-percent win if it can be composed into a larger
hybrid, but not a direct SOTA path by itself.

## Bottom line

After correcting the classical moonshot analysis, the strongest honest summary is:

> We do not currently know a literature-backed algorithmic replacement for
> Kaliski that plausibly beats our existing 4.39M-Toffoli circuit at n=256.
> The remaining gap to Google's 2.1–2.7M number most likely depends on
> undisclosed techniques or genuinely new reversible-algorithm research.
