//! Classical analysis for a possible **hybrid Kaliski-jump** moonshot.
//!
//! Idea: keep the existing Kaliski cleanup machinery (`r,s,m_hist`) but try to
//! batch *local* `(u, v)` updates over a small fixed number of steps `t`, keyed
//! by the low `w` bits of `(u, v)`. If the resulting t-step transition matrices
//! come from a very small family per low-bit class, then a compressed QROM could
//! replace several expensive per-step parity/compare/cswap/sub/halve operations.
//!
//! This file is **classical-only** research infrastructure. It does not affect
//! the quantum circuit.
//!
//! Standard almost-inverse / binary-GCD step on nonnegative integers:
//!
//! ```text
//! if u even:                   (u, v) ← (u/2, v)
//! elif v even:                (u, v) ← (u, v/2)
//! elif u > v:                 (u, v) ← ((u-v)/2, v)
//! else:                       (u, v) ← (u, (v-u)/2)
//! ```
//!
//! Each branch can be represented as a linear map with a shared `1/2` factor:
//!
//! ```text
//! U-even:  (u', v') = (1/2) * [[1,  0], [0, 2]] * (u, v)
//! V-even:  (u', v') = (1/2) * [[2,  0], [0, 1]] * (u, v)
//! U>V:     (u', v') = (1/2) * [[1, -1], [0, 2]] * (u, v)
//! V>U:     (u', v') = (1/2) * [[2,  0], [-1,1]] * (u, v)
//! ```
//!
//! Over `t` steps, we can accumulate an integer matrix `P_t` such that:
//!
//! ```text
//! (u_t, v_t)^T = (1 / 2^t) * P_t * (u_0, v_0)^T
//! ```
//!
//! The research questions here are:
//! 1. How many distinct `P_t` appear along actual secp256k1 Kaliski trajectories?
//! 2. For a fixed low-bit class `(u mod 2^w, v mod 2^w)`, how many distinct
//!    `P_t` values occur? If this is very small, a compressed lookup might work.
//! 3. How big do the entries of `P_t` get in practice (vs. the trivial 2^t bound)?

use std::collections::{BTreeMap, BTreeSet};

use alloy_primitives::U256;
use sha3::digest::{ExtendableOutput, Update, XofReader};

use super::SECP256K1_P;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Mat2 {
    pub a00: i128,
    pub a01: i128,
    pub a10: i128,
    pub a11: i128,
}

impl Mat2 {
    pub const ID: Mat2 = Mat2 { a00: 1, a01: 0, a10: 0, a11: 1 };

    pub fn mul(self, rhs: Mat2) -> Mat2 {
        Mat2 {
            a00: self.a00 * rhs.a00 + self.a01 * rhs.a10,
            a01: self.a00 * rhs.a01 + self.a01 * rhs.a11,
            a10: self.a10 * rhs.a00 + self.a11 * rhs.a10,
            a11: self.a10 * rhs.a01 + self.a11 * rhs.a11,
        }
    }

    pub fn max_abs(&self) -> i128 {
        [self.a00.abs(), self.a01.abs(), self.a10.abs(), self.a11.abs()]
            .into_iter()
            .max()
            .unwrap_or(0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum KCase {
    UEven,
    VEven,
    UGtV,
    VGtU,
}

impl KCase {
    /// Matrix for the (u, v) register update over one Kaliski micro-step:
    ///
    ///   (u', v')^T = (1/2) · M_uv · (u, v)^T
    pub fn uv_matrix(self) -> Mat2 {
        match self {
            KCase::UEven => Mat2 { a00: 1, a01: 0, a10: 0, a11: 2 },
            KCase::VEven => Mat2 { a00: 2, a01: 0, a10: 0, a11: 1 },
            KCase::UGtV  => Mat2 { a00: 1, a01: -1, a10: 0, a11: 2 },
            KCase::VGtU  => Mat2 { a00: 2, a01: 0, a10: -1, a11: 1 },
        }
    }

    /// Matrix for the coefficient-side (r, s) update over one Kaliski step.
    ///
    /// Derived directly from the implemented sequence in `kaliski_iteration`:
    ///
    /// - UEven:  swap(r,s); double r; swap back  =>  (r, s) -> (r, 2s)
    /// - VEven:  double r                        =>  (r, s) -> (2r, s)
    /// - UGtV:   swap; s += r; double r; swap   =>  (r, s) -> (r+s, 2s)
    /// - VGtU:   s += r; double r                =>  (r, s) -> (2r, r+s)
    pub fn rs_matrix(self) -> Mat2 {
        match self {
            KCase::UEven => Mat2 { a00: 1, a01: 0, a10: 0, a11: 2 },
            KCase::VEven => Mat2 { a00: 2, a01: 0, a10: 0, a11: 1 },
            KCase::UGtV  => Mat2 { a00: 1, a01: 1, a10: 0, a11: 2 },
            KCase::VGtU  => Mat2 { a00: 2, a01: 0, a10: 1, a11: 1 },
        }
    }
}
#[inline(always)]
fn kaliski_case(u: U256, v: U256) -> KCase {
    if !u.bit(0) {
        KCase::UEven
    } else if !v.bit(0) {
        KCase::VEven
    } else if u > v {
        KCase::UGtV
    } else {
        KCase::VGtU
    }
}

#[inline(always)]
fn kaliski_step_uv(u: U256, v: U256) -> (U256, U256, KCase) {
    match kaliski_case(u, v) {
        KCase::UEven => (u >> 1, v, KCase::UEven),
        KCase::VEven => (u, v >> 1, KCase::VEven),
        KCase::UGtV  => ((u.wrapping_sub(v)) >> 1, v, KCase::UGtV),
        KCase::VGtU  => (u, (v.wrapping_sub(u)) >> 1, KCase::VGtU),
    }
}

#[derive(Clone, Debug)]
pub struct WindowObs {
    pub low_u: u16,
    pub low_v: u16,
    pub uv_mat: Mat2,
    pub rs_mat: Mat2,
    pub cases: Vec<KCase>,
}

/// Observe a t-step Kaliski window starting from full-width `(u, v)`.
/// Returns `(u_t, v_t, obs)`.
pub fn observe_window(mut u: U256, mut v: U256, w: usize, t: usize) -> (U256, U256, WindowObs) {
    assert!(w <= 16, "low-bit class currently stored as u16");
    let low_mask = if w == 16 {
        U256::from(0xFFFFu64)
    } else {
        (U256::from(1u64) << w).wrapping_sub(U256::from(1u64))
    };
    let low_u = (u & low_mask).to::<u16>();
    let low_v = (v & low_mask).to::<u16>();
    let mut uv_mat = Mat2::ID;
    let mut rs_mat = Mat2::ID;
    let mut cases = Vec::with_capacity(t);
    for _ in 0..t {
        if v.is_zero() { break; }
        let (nu, nv, kc) = kaliski_step_uv(u, v);
        uv_mat = kc.uv_matrix().mul(uv_mat);
        rs_mat = kc.rs_matrix().mul(rs_mat);
        cases.push(kc);
        u = nu;
        v = nv;
    }
    (u, v, WindowObs { low_u, low_v, uv_mat, rs_mat, cases })
}

pub struct Sampler {
    reader: Box<dyn XofReader>,
    p: U256,
}

impl Sampler {
    pub fn new(seed: &[u8], p: U256) -> Self {
        let mut hasher = sha3::Shake128::default();
        hasher.update(seed);
        Self { reader: Box::new(hasher.finalize_xof()), p }
    }

    pub fn next(&mut self) -> U256 {
        loop {
            let mut buf = [0u8; 32];
            self.reader.read(&mut buf);
            let x = U256::from_le_slice(&buf);
            if x < self.p && !x.is_zero() {
                return x;
            }
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct HybridStats {
    pub inputs: usize,
    pub windows: usize,

    pub distinct_global_uv_mats: usize,
    pub distinct_global_rs_mats: usize,

    pub max_uv_entry_abs: i128,
    pub mean_log2_uv_entry_abs: f64,
    pub max_rs_entry_abs: i128,
    pub mean_log2_rs_entry_abs: f64,

    pub low_classes_seen: usize,
    pub mean_uv_mats_per_class: f64,
    pub max_uv_mats_per_class: usize,
    pub singleton_uv_classes: usize,
    pub most_common_uv_class_count: usize,
    pub most_common_uv_class: Option<(u16, u16)>,

    pub mean_rs_mats_per_class: f64,
    pub max_rs_mats_per_class: usize,
    pub singleton_rs_classes: usize,
    pub most_common_rs_class_count: usize,
    pub most_common_rs_class: Option<(u16, u16)>,

    // Joint compression: how many distinct (uv_mat, rs_mat) pairs actually occur.
    pub distinct_joint_pairs: usize,
    pub mean_joint_pairs_per_class: f64,
    pub max_joint_pairs_per_class: usize,
    pub singleton_joint_classes: usize,
    pub most_common_joint_class_count: usize,
    pub most_common_joint_class: Option<(u16, u16)>,
}

/// Sample actual secp256k1 Kaliski trajectories and measure the compressibility
/// of t-step local transition matrices keyed by low-w bits.
pub fn hybrid_kaliski_window_survey(
    seed: &[u8],
    n_inputs: usize,
    w: usize,
    t: usize,
) -> HybridStats {
    let mut sampler = Sampler::new(seed, SECP256K1_P);
    let mut global_uv_mats: BTreeSet<Mat2> = BTreeSet::new();
    let mut global_rs_mats: BTreeSet<Mat2> = BTreeSet::new();
    let mut global_joint_pairs: BTreeSet<(Mat2, Mat2)> = BTreeSet::new();
    let mut by_class_uv: BTreeMap<(u16, u16), BTreeSet<Mat2>> = BTreeMap::new();
    let mut by_class_rs: BTreeMap<(u16, u16), BTreeSet<Mat2>> = BTreeMap::new();
    let mut by_class_joint: BTreeMap<(u16, u16), BTreeSet<(Mat2, Mat2)>> = BTreeMap::new();
    let mut windows = 0usize;
    let mut max_uv_entry_abs = 0i128;
    let mut sum_log2_uv_entry_abs = 0.0f64;
    let mut counted_uv_mats = 0usize;
    let mut max_rs_entry_abs = 0i128;
    let mut sum_log2_rs_entry_abs = 0.0f64;
    let mut counted_rs_mats = 0usize;

    for _ in 0..n_inputs {
        let mut u = SECP256K1_P;
        let mut v = sampler.next();
        for _ in 0..742 {
            if v.is_zero() { break; }
            let (nu, nv, obs) = observe_window(u, v, w, t);
            global_uv_mats.insert(obs.uv_mat);
            global_rs_mats.insert(obs.rs_mat);
            global_joint_pairs.insert((obs.uv_mat, obs.rs_mat));
            by_class_uv.entry((obs.low_u, obs.low_v)).or_default().insert(obs.uv_mat);
            by_class_rs.entry((obs.low_u, obs.low_v)).or_default().insert(obs.rs_mat);
            by_class_joint.entry((obs.low_u, obs.low_v)).or_default().insert((obs.uv_mat, obs.rs_mat));
            let uv_abs = obs.uv_mat.max_abs();
            if uv_abs > max_uv_entry_abs { max_uv_entry_abs = uv_abs; }
            if uv_abs > 0 {
                sum_log2_uv_entry_abs += (uv_abs as f64).log2();
                counted_uv_mats += 1;
            }
            let rs_abs = obs.rs_mat.max_abs();
            if rs_abs > max_rs_entry_abs { max_rs_entry_abs = rs_abs; }
            if rs_abs > 0 {
                sum_log2_rs_entry_abs += (rs_abs as f64).log2();
                counted_rs_mats += 1;
            }
            windows += 1;
            let (u1, v1, _kc) = kaliski_step_uv(u, v);
            u = u1;
            v = v1;
        }
    }

    let low_classes_seen = by_class_uv.len();

    let mut total_uv_mats_per_class = 0usize;
    let mut max_uv_mats_per_class = 0usize;
    let mut singleton_uv_classes = 0usize;
    let mut most_common_uv_class_count = 0usize;
    let mut most_common_uv_class = None;
    for (cls, mats) in &by_class_uv {
        let c = mats.len();
        total_uv_mats_per_class += c;
        if c > max_uv_mats_per_class { max_uv_mats_per_class = c; }
        if c == 1 { singleton_uv_classes += 1; }
        if c > most_common_uv_class_count {
            most_common_uv_class_count = c;
            most_common_uv_class = Some(*cls);
        }
    }

    let mut total_rs_mats_per_class = 0usize;
    let mut max_rs_mats_per_class = 0usize;
    let mut singleton_rs_classes = 0usize;
    let mut most_common_rs_class_count = 0usize;
    let mut most_common_rs_class = None;
    for (cls, mats) in &by_class_rs {
        let c = mats.len();
        total_rs_mats_per_class += c;
        if c > max_rs_mats_per_class { max_rs_mats_per_class = c; }
        if c == 1 { singleton_rs_classes += 1; }
        if c > most_common_rs_class_count {
            most_common_rs_class_count = c;
            most_common_rs_class = Some(*cls);
        }
    }

    let mut total_joint_pairs_per_class = 0usize;
    let mut max_joint_pairs_per_class = 0usize;
    let mut singleton_joint_classes = 0usize;
    let mut most_common_joint_class_count = 0usize;
    let mut most_common_joint_class = None;
    for (cls, mats) in &by_class_joint {
        let c = mats.len();
        total_joint_pairs_per_class += c;
        if c > max_joint_pairs_per_class { max_joint_pairs_per_class = c; }
        if c == 1 { singleton_joint_classes += 1; }
        if c > most_common_joint_class_count {
            most_common_joint_class_count = c;
            most_common_joint_class = Some(*cls);
        }
    }

    HybridStats {
        inputs: n_inputs,
        windows,
        distinct_global_uv_mats: global_uv_mats.len(),
        distinct_global_rs_mats: global_rs_mats.len(),
        max_uv_entry_abs,
        mean_log2_uv_entry_abs: if counted_uv_mats == 0 { 0.0 } else { sum_log2_uv_entry_abs / counted_uv_mats as f64 },
        max_rs_entry_abs,
        mean_log2_rs_entry_abs: if counted_rs_mats == 0 { 0.0 } else { sum_log2_rs_entry_abs / counted_rs_mats as f64 },
        low_classes_seen,
        mean_uv_mats_per_class: if low_classes_seen == 0 { 0.0 } else { total_uv_mats_per_class as f64 / low_classes_seen as f64 },
        max_uv_mats_per_class,
        singleton_uv_classes,
        most_common_uv_class_count,
        most_common_uv_class,
        mean_rs_mats_per_class: if low_classes_seen == 0 { 0.0 } else { total_rs_mats_per_class as f64 / low_classes_seen as f64 },
        max_rs_mats_per_class,
        singleton_rs_classes,
        most_common_rs_class_count,
        most_common_rs_class,
        distinct_joint_pairs: global_joint_pairs.len(),
        mean_joint_pairs_per_class: if low_classes_seen == 0 { 0.0 } else { total_joint_pairs_per_class as f64 / low_classes_seen as f64 },
        max_joint_pairs_per_class,
        singleton_joint_classes,
        most_common_joint_class_count,
        most_common_joint_class,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_window_smoke() {
        let u = SECP256K1_P;
        let v = U256::from(123456789u64);
        let (_u2, _v2, obs) = observe_window(u, v, 8, 4);
        assert!(obs.cases.len() >= 1);
        assert!(obs.uv_mat.max_abs() >= 1);
        assert!(obs.rs_mat.max_abs() >= 1);
    }

    #[test]
    fn hybrid_kaliski_window_survey_test() {
        for &(w, t) in &[(6usize, 4usize), (8usize, 4usize), (8usize, 6usize)] {
            let s = hybrid_kaliski_window_survey(b"hybrid-kaliski-window-seed-v1", 10_000, w, t);
            eprintln!("=== hybrid Kaliski window survey (w={}, t={}) ===", w, t);
            eprintln!("inputs                  : {}", s.inputs);
            eprintln!("windows                 : {}", s.windows);
            eprintln!("distinct global uv mats : {}", s.distinct_global_uv_mats);
            eprintln!("distinct global rs mats : {}", s.distinct_global_rs_mats);
            eprintln!("max |uv entry|          : {}", s.max_uv_entry_abs);
            eprintln!("mean log2 |uv entry|    : {:.3}", s.mean_log2_uv_entry_abs);
            eprintln!("max |rs entry|          : {}", s.max_rs_entry_abs);
            eprintln!("mean log2 |rs entry|    : {:.3}", s.mean_log2_rs_entry_abs);
            eprintln!("classes seen            : {}", s.low_classes_seen);
            eprintln!("mean uv mats/class      : {:.3}", s.mean_uv_mats_per_class);
            eprintln!("max uv mats/class       : {}", s.max_uv_mats_per_class);
            eprintln!("singleton uv classes    : {}", s.singleton_uv_classes);
            eprintln!("most common uv class ct : {}", s.most_common_uv_class_count);
            if let Some((ucls, vcls)) = s.most_common_uv_class {
                eprintln!("most common uv class    : (u_low={}, v_low={})", ucls, vcls);
            }
            eprintln!("mean rs mats/class      : {:.3}", s.mean_rs_mats_per_class);
            eprintln!("max rs mats/class       : {}", s.max_rs_mats_per_class);
            eprintln!("singleton rs classes    : {}", s.singleton_rs_classes);
            eprintln!("most common rs class ct : {}", s.most_common_rs_class_count);
            if let Some((ucls, vcls)) = s.most_common_rs_class {
                eprintln!("most common rs class    : (u_low={}, v_low={})", ucls, vcls);
            }
            eprintln!("distinct joint pairs    : {}", s.distinct_joint_pairs);
            eprintln!("mean joint/class        : {:.3}", s.mean_joint_pairs_per_class);
            eprintln!("max joint/class         : {}", s.max_joint_pairs_per_class);
            eprintln!("singleton joint classes : {}", s.singleton_joint_classes);
            eprintln!("most common joint ct    : {}", s.most_common_joint_class_count);
            if let Some((ucls, vcls)) = s.most_common_joint_class {
                eprintln!("most common joint class : (u_low={}, v_low={})", ucls, vcls);
            }
            eprintln!("===============================================");
            assert!(s.windows > 0);
            assert!(s.distinct_global_uv_mats >= 1);
            assert!(s.distinct_global_rs_mats >= 1);
            assert!(s.distinct_joint_pairs >= 1);
        }
    }
}
