# Observed bulk-prefix phase pattern

Current strict passing `k` values confirmed through `main.rs`:
- 3, 6, 7, 24, 32, 72, 96

Current strict failing `k` values confirmed through `main.rs`:
- 4, 5, 8, 16, 40, 64, 80, 100, 104, 108, 112, 120, 124, 128

## Immediate pattern observations
- All passing values >= 24 are multiples of **24** or **32**, except none so far.
- Passing set includes:
  - 24 = 3 * 8
  - 32 = 1 * 32
  - 72 = 3 * 24
  - 96 = lcm(24, 32)
- Failing values often but not always lie on nearby multiples of 8.
- Notably:
  - 120 (a multiple of 24) fails,
  - 64 and 128 (powers of two) fail,
  - so neither “multiple of 24” nor “multiple of 32” alone explains the set.

## Working interpretation
The data does **not** look random. It looks like a coherent cancellation pattern
with at least two interacting periodic components.

The simplest candidate is that the specialized prefix injects a hidden phase
contribution whose cancellation depends on the prefix length modulo multiple
small periods, likely involving both:
- a period tied to the bulk step structure itself,
- and another tied to the surrounding scaffold / cleanup composition.

This is consistent with the currently best passing values at 24, 32, 72, 96.
