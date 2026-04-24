//! Classical prototype / qubit-budget scratchpad for Luo-style register sharing.
//!
//! Goal of this file: validate *cheaply* whether a Luo/PZ-style inversion track
//! is even compatible with the user's budget:
//!   - only ~600 qubits over the 512 input-point-coordinate qubits,
//!   - i.e. total target around 1100–1200 for the full point-add, or at least
//!     meaningfully below our current 2716q.
//!
//! We do NOT attempt a reversible circuit here. We only:
//!   1. model the qubit budget implied by Luo's register sharing,
//!   2. compare it to our current Kaliski scaffold,
//!   3. record the minimum structural consequences.
//!
//! Literature anchor: `/tmp/luo_ec_clean.txt`, especially Table 1 and Algorithm 3.

#![cfg(test)]

use alloy_primitives::U256;

use super::SECP256K1_P;

/// Very coarse qubit budget for the current live affine scaffold.
#[derive(Debug, Clone, Copy)]
struct Budget {
    tx_ty: usize,
    inversion_state: usize,
    lambda_and_mul_state: usize,
    classical_bits: usize,
    total: usize,
}

/// Current live build (best stable before the 511 detour):
/// - tx,ty = 2n = 512
/// - Kaliski persistent state = u,v_w,r,s,m_hist,f ≈ 4n + iters + 1
///   with iters ≈ 407/404 → ~1432
/// - live body state + mul transients explain the observed 2716 peak.
fn current_budget_estimate(n: usize, iters: usize) -> Budget {
    let tx_ty = 2 * n;
    let inversion_state = 4 * n + iters + 1;
    // Remaining gap to the observed peak (2716) is dominated by lam,
    // tmp_ext, carries, and a few flags.
    let total = 2716;
    let lambda_and_mul_state = total - tx_ty - inversion_state;
    Budget {
        tx_ty,
        inversion_state,
        lambda_and_mul_state,
        classical_bits: 2 * n, // ox, oy
        total,
    }
}

/// Luo-style inversion state from `/tmp/luo_ec_clean.txt`:
/// Table 1 says inversion can be done in roughly `3n + 4 log2 n + O(1)`
/// qubits *total* for the inversion component.
///
/// For n=256 this is about 3*256 + 4*8 = 800 qubits total, INCLUDING the
/// input/output pair of the inversion itself.
///
/// In our point-add context the inversion input is one n-bit value (dx or
/// similar), and we still need the 2n point coordinates live. So the key
/// number is the non-tx/ty overhead: about n + O(log n), not 4n+iters.
fn luo_inversion_total_qubits(n: usize) -> usize {
    3 * n + 4 * (n.ilog2() as usize)
}

/// Lower-bound point-add budget if we swapped our Kaliski block for a Luo/PZ
/// block *without* changing anything else in the affine scaffold.
fn naive_luo_point_add_budget(n: usize, current_other_peak: usize) -> usize {
    // Keep tx,ty live. Add the Luo inversion block. Keep the rest of the
    // current non-inversion transients as-is.
    let tx_ty = 2 * n;
    let inversion_total = luo_inversion_total_qubits(n);
    // current_other_peak is "everything except tx_ty and current inversion".
    tx_ty + inversion_total + current_other_peak
}

/// Clean arithmetic helper for a tiny classical sanity check.
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

    #[test]
    fn luo_budget_is_qubit_relevant() {
        let n = 256usize;
        let cur = current_budget_estimate(n, 407);
        eprintln!("current budget estimate: {:?}", cur);

        // Everything that's NOT tx/ty or persistent Kaliski state.
        let current_other_peak = cur.lambda_and_mul_state;
        let luo_total = naive_luo_point_add_budget(n, current_other_peak);
        eprintln!("naive Luo swap-in peak estimate: {luo_total}");

        // Current live peak is 2716. Luo Table 1 claims inversion itself is
        // only ~800 qubits total, vs our current inversion state + body share
        // of ~1948 qubits. Even the naive swap should cut a lot.
        assert!(luo_total < cur.total, "Luo-style inversion must reduce peak");

        // User asked for only ~600 qubits over the 512 input point coords,
        // i.e. roughly 1112 total. A naive Luo swap alone will NOT get us
        // there — the rest of the affine scaffold is still too wide.
        assert!(
            luo_total > 1112,
            "If this ever flips, Luo-alone got us into the user budget; revisit immediately"
        );
    }

    #[test]
    fn luo_alone_is_not_sota_but_is_structural() {
        let n = 256usize;
        let cur = current_budget_estimate(n, 407);
        let current_other_peak = cur.lambda_and_mul_state;
        let luo_total = naive_luo_point_add_budget(n, current_other_peak);

        // This is the key conclusion for next actions:
        // Luo's register sharing plausibly drops us from ~2716 to ~2k-ish,
        // but not all the way to ~1100 by itself. So Luo is only worth it if
        // paired with an affine-scaffold collapse (single inversion,
        // different coordinate flow, or killing lam/m_hist simultaneously).
        eprintln!(
            "luo_total={}, current_total={}, saved={} qubits",
            luo_total,
            cur.total,
            cur.total - luo_total
        );
        assert!(cur.total - luo_total >= 500, "Luo should save ~500+ qubits");
    }

    #[test]
    fn dy_py_relation_sanity() {
        // Tiny guard against the kind of algebra drift we had earlier.
        let p = SECP256K1_P;
        let px = U256::from(123u64);
        let py = U256::from(456u64);
        let qx = U256::from(17u64);
        let qy = U256::from(31u64);
        let dx = sub_mod(px, qx, p);
        let dy = sub_mod(py, qy, p);
        assert_eq!(dx, U256::from(106u64));
        assert_eq!(dy, U256::from(425u64));
    }
}
