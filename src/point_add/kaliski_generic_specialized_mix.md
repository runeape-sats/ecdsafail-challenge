# Proposed isolated phase-correct specialization strategy

The direct specialized bulk-prefix body keeps reintroducing phase bugs because
it tries to be a new circuit while still composing with the old measurement-
uncompute scaffold.

A more isolated strategy is:
- keep the **generic Kaliski skeleton** exactly,
- only specialize the arithmetic sub-operations that are phase-neutral, and
- never alter the HMR / feed-forward / cleanup structure.

## Practical meaning
Instead of replacing the whole `kaliski_iteration` body, only replace pieces
like:
- the arithmetic in STEP 4,
- safe `mod_double_no_corr` uses,
- safe width truncations,
- or known-constant controls,

while preserving:
- the same ancilla allocation order,
- the same HMR sequence,
- the same CZ/CCZ feed-forward locations,
- the same step skeleton forward and backward.

This is the lowest-risk path to a truly phase-correct optimization, because it
turns the problem from “invent a new measurement-uncompute protocol” into “keep
the old protocol, specialize only the phase-insensitive arithmetic.”
