//! Classical / qubit-budget prototype for Kim 2026 style unconditional Kaliski.
//!
//! Purpose: validate *fast* whether Kim's unconditional-execution trick can be
//! imported into our scaffold, and if not, what exact statement survives.
//!
//! Key correction vs earlier monologue: the naive claim
//!   "just keep stepping after v=0 and r only doubles"
//! is FALSE under our current 256-bit `U256`-truncated classical model,
//! because the Kim paper explicitly postpones reduction into a 2n-bit `r`.
//! So the right prototype uses a wide `U512`-style accumulator, not `U256`.
//!
//! Current status of this prototype:
//! - `dy_over_dx_reference_sanity` is a live sanity check.
//! - `naive_*` tests are kept as ignored negative results for the old wrong
//!   formulation.
//! - `wide_unconditional_exec_*` are the real tests for whether Kim is still
//!   alive in a widened-r model.

#![cfg(test)]

use alloy_primitives::{U256, U512};

use super::SECP256K1_P;

#[derive(Clone, Debug)]
struct St {
    u: U256,
    v: U256,
    r: U512,
    s: U512,
}

fn u256_to_u512(x: U256) -> U512 {
    U512::from_limbs([
        x.as_limbs()[0],
        x.as_limbs()[1],
        x.as_limbs()[2],
        x.as_limbs()[3],
        0,
        0,
        0,
        0,
    ])
}

fn low_u256(x: U512) -> U256 {
    let limbs = x.as_limbs();
    U256::from_limbs([limbs[0], limbs[1], limbs[2], limbs[3]])
}

fn mod_p_from_u512(x: U512) -> U256 {
    // Slow but fine for tests. Convert via bytes then reduce by repeated fold.
    // Since U512 does not expose a direct mod-U256 helper in our codebase, we
    // use a simple shift/add reduction on bytes.
    let bytes = x.to_le_bytes::<64>();
    let lo = U256::from_le_slice(&bytes[0..32]);
    let hi = U256::from_le_slice(&bytes[32..64]);
    // x = lo + 2^256 * hi ≡ lo + (2^32 + 977) * hi mod p, because
    // 2^256 ≡ 2^32 + 977 mod p for secp256k1.
    let p = SECP256K1_P;
    let c = U256::from(1u64 << 32).add_mod(U256::from(977u64), p);
    lo.add_mod(hi.mul_mod(c, p), p)
}

/// Conditional step in the *current* branch logic, but keeping r,s wide.
fn conditional_step_wide(st: &mut St) -> bool {
    if st.v.is_zero() {
        return false;
    }
    let u = st.u;
    let v = st.v;
    let r = st.r;
    let s = st.s;
    if !u.bit(0) {
        st.u = u >> 1;
        st.v = v;
        st.r = r;
        st.s = s << 1;
    } else if !v.bit(0) {
        st.u = u;
        st.v = v >> 1;
        st.r = r << 1;
        st.s = s;
    } else if u > v {
        st.u = (u.wrapping_sub(v)) >> 1;
        st.v = v;
        st.r = r + s;
        st.s = s << 1;
    } else {
        st.u = u;
        st.v = (v.wrapping_sub(u)) >> 1;
        st.r = r << 1;
        st.s = r + s;
    }
    true
}

/// Unconditional extension after v=0 in the *Kim-style wide-r model*:
/// keep the same round logic, but when v=0 we only apply the residual
/// doubling on r. This is the exact claim we want to test numerically.
fn unconditional_step_wide(st: &mut St) {
    if st.v.is_zero() {
        st.r <<= 1;
        return;
    }
    let _ = conditional_step_wide(st);
}

fn run_conditional_wide(v0: U256, max_steps: usize) -> (St, usize) {
    let mut st = St {
        u: SECP256K1_P,
        v: v0,
        r: U512::ZERO,
        s: U512::from(1u64),
    };
    let mut k = 0usize;
    while k < max_steps && conditional_step_wide(&mut st) {
        k += 1;
    }
    (st, k)
}

fn run_unconditional_wide(v0: U256, rounds: usize) -> St {
    let mut st = St {
        u: SECP256K1_P,
        v: v0,
        r: U512::ZERO,
        s: U512::from(1u64),
    };
    for _ in 0..rounds {
        unconditional_step_wide(&mut st);
    }
    st
}

fn sub_mod(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        (a - b) % p
    } else {
        p - ((b - a) % p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

    fn curve() -> WeierstrassEllipticCurve {
        WeierstrassEllipticCurve {
            modulus: SECP256K1_P,
            a: U256::from(0),
            b: U256::from(7),
            gx: U256::from_str_radix(
                "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
                16,
            )
            .unwrap(),
            gy: U256::from_str_radix(
                "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
                16,
            )
            .unwrap(),
            order: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
                16,
            )
            .unwrap(),
        }
    }

    fn rand_u256(rng: &mut u64) -> U256 {
        let mut limbs = [0u64; 4];
        for l in &mut limbs {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *l = *rng;
        }
        U256::from_limbs(limbs) % SECP256K1_P
    }

    #[test]
    #[ignore = "negative result from earlier wrong narrow-r model"]
    fn naive_unconditional_exec_turns_dynamic_correction_into_fixed_tail() {}

    #[test]
    #[ignore = "negative result from earlier wrong narrow-r model"]
    fn naive_unconditional_exec_keeps_scale_deterministic_at_2n() {}

    #[test]
    #[ignore = "negative result from earlier wrong narrow-r model"]
    fn naive_pair1_pair2_correction_loops_are_exactly_the_dynamic_tail_today() {}

    #[test]
    fn wide_unconditional_exec_tail_matches_fixed_doubling() {
        let mut rng = 0x1234_5678_9abc_def0u64;
        for _ in 0..200 {
            let mut x = rand_u256(&mut rng);
            while x.is_zero() {
                x = rand_u256(&mut rng);
            }
            let (cond, k) = run_conditional_wide(x, 2 * 256);
            assert!(k >= 256 && k <= 511);
            let uncond = run_unconditional_wide(x, 2 * 256);
            let expected_r = cond.r << (2 * 256 - k);
            assert_eq!(
                uncond.r, expected_r,
                "wide unconditional tail is not fixed-count doubling"
            );
        }
    }

    #[test]
    fn wide_unconditional_exec_final_low_word_has_fixed_scale() {
        let mut rng = 0x0ddc0ffee1234567u64;
        let p = SECP256K1_P;
        let two = U256::from(2);
        let scale_2n = two.pow_mod(U256::from(512u64), p);

        for _ in 0..100 {
            let mut x = rand_u256(&mut rng);
            while x.is_zero() {
                x = rand_u256(&mut rng);
            }
            let st = run_unconditional_wide(x, 512);
            let low = mod_p_from_u512(st.r);
            let expect = x.inv_mod(p).unwrap().mul_mod(scale_2n, p);
            // Sign is intentionally NOT asserted yet — the classical branch
            // convention here is not proven to match our quantum-sign choice.
            let expect_neg = sub_mod(U256::ZERO, expect, p);
            assert!(
                low == expect || low == expect_neg,
                "wide unconditional low word is not ±x^-1 * 2^(2n)"
            );
        }
    }

    #[test]
    fn wide_conditional_k_range_is_tight() {
        let mut rng = 0xabcdef0123456789u64;
        let mut min_k = usize::MAX;
        let mut max_k = 0usize;
        let mut total_k = 0usize;
        for _ in 0..200 {
            let mut x = rand_u256(&mut rng);
            while x.is_zero() {
                x = rand_u256(&mut rng);
            }
            let (_st, k) = run_conditional_wide(x, 512);
            min_k = min_k.min(k);
            max_k = max_k.max(k);
            total_k += k;
        }
        let avg_k = total_k as f64 / 200.0;
        eprintln!("wide conditional Kaliski termination k range over 200 samples: [{min_k}, {max_k}], avg={avg_k:.2}");
        assert!(min_k >= 256);
        assert!(max_k <= 511);
        assert!(avg_k > 330.0 && avg_k < 390.0);
    }

    #[test]
    fn dy_over_dx_reference_sanity() {
        let c = curve();
        let (px, py) = c.mul(c.gx, c.gy, U256::from(11u64));
        let (qx, qy) = c.mul(c.gx, c.gy, U256::from(19u64));
        let dx = sub_mod(px, qx, SECP256K1_P);
        let dy = sub_mod(py, qy, SECP256K1_P);
        let lam = dy.mul_mod(dx.inv_mod(SECP256K1_P).unwrap(), SECP256K1_P);
        let (rx, _ry) = c.add(px, py, qx, qy);
        let rx_formula = sub_mod(
            sub_mod(lam.mul_mod(lam, SECP256K1_P), px, SECP256K1_P),
            qx,
            SECP256K1_P,
        );
        assert_eq!(rx, rx_formula);
    }
}
