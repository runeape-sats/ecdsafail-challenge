//! Diagnostic equivalence checks for experimental Kaliski bulk primitives.
//!
//! This is not a performance profiler. Its job is to answer a sharper question:
//!
//! > Does an experimental forward bulk primitive produce the *same persistent
//! > Kaliski state* as the generic step on actual reachable secp256k1 states?
//!
//! That is the right tool for integrating something real into the live circuit,
//! because the current blocker is backward/history compatibility, not lack of
//! raw cost numbers.

use alloy_primitives::U256;
use sha3::digest::{ExtendableOutput, Update};

use crate::sim::Simulator;

use super::kaliski_jump::Sampler;
use super::test_timeout::{check_deadline, two_min_deadline};
use super::{
    kaliski_iteration, kaliski_iteration_bulk_prefix3, B, N, Op, QubitId, SECP256K1_P,
};

#[derive(Clone)]
struct StepCircuit {
    ops: Vec<Op>,
    num_qubits: usize,
    num_bits: usize,
    u: Vec<QubitId>,
    v: Vec<QubitId>,
    r: Vec<QubitId>,
    s: Vec<QubitId>,
    m: QubitId,
    f: QubitId,
}

fn build_generic_step(iter_idx: usize) -> StepCircuit {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    let f = b.alloc_qubit();
    kaliski_iteration(&mut b, SECP256K1_P, &u, &v, &r, &s, m, f, iter_idx);
    StepCircuit {
        ops: b.ops,
        num_qubits: b.next_qubit as usize,
        num_bits: b.next_bit as usize,
        u, v, r, s, m, f,
    }
}

fn build_special_step(iter_idx: usize) -> StepCircuit {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    let f = b.alloc_qubit();
    kaliski_iteration_bulk_prefix3(&mut b, &u, &v, &r, &s, m, iter_idx);
    StepCircuit {
        ops: b.ops,
        num_qubits: b.next_qubit as usize,
        num_bits: b.next_bit as usize,
        u, v, r, s, m, f,
    }
}

fn build_special_three_steps() -> StepCircuit {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let f = b.alloc_qubit();
    let m_hist = b.alloc_qubits(3);
    for i in 0..3 {
        kaliski_iteration_bulk_prefix3(&mut b, &u, &v, &r, &s, m_hist[i], i);
    }
    StepCircuit {
        ops: b.ops,
        num_qubits: b.next_qubit as usize,
        num_bits: b.next_bit as usize,
        u, v, r, s, m: m_hist[2], f,
    }
}

fn set_slice<R: sha3::digest::XofReader>(sim: &mut Simulator<R>, qs: &[QubitId], val: U256) {
    for (i, &q) in qs.iter().enumerate() {
        if val.bit(i) {
            *sim.qubit_mut(q) |= 1;
        } else {
            *sim.qubit_mut(q) &= !1;
        }
    }
}

fn get_slice<R: sha3::digest::XofReader>(sim: &Simulator<R>, qs: &[QubitId]) -> U256 {
    let mut out = U256::ZERO;
    for (i, &q) in qs.iter().enumerate() {
        out.set_bit(i, (sim.qubit(q) & 1) != 0);
    }
    out
}

fn run_step_circuit(c: &StepCircuit, u0: U256, v0: U256, r0: U256, s0: U256, m0: bool, f0: bool)
    -> (U256, U256, U256, U256, bool, bool)
{
    let mut hasher = sha3::Shake128::default();
    hasher.update(b"kaliski-equivalence-seed-v1");
    let mut xof = hasher.finalize_xof();
    let mut sim = Simulator::new(c.num_qubits, c.num_bits, &mut xof);
    set_slice(&mut sim, &c.u, u0);
    set_slice(&mut sim, &c.v, v0);
    set_slice(&mut sim, &c.r, r0);
    set_slice(&mut sim, &c.s, s0);
    if m0 { *sim.qubit_mut(c.m) |= 1; }
    if f0 { *sim.qubit_mut(c.f) |= 1; }
    sim.apply(&c.ops);
    (
        get_slice(&sim, &c.u),
        get_slice(&sim, &c.v),
        get_slice(&sim, &c.r),
        get_slice(&sim, &c.s),
        (sim.qubit(c.m) & 1) != 0,
        (sim.qubit(c.f) & 1) != 0,
    )
}

#[derive(Clone, Debug)]
struct ClassicalState {
    u: U256,
    v: U256,
    r: U256,
    s: U256,
}

fn classical_step(st: &mut ClassicalState) -> bool {
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
        st.r = r.wrapping_add(s);
        st.s = s << 1;
    } else {
        st.u = u;
        st.v = (v.wrapping_sub(u)) >> 1;
        st.r = r << 1;
        st.s = r.wrapping_add(s);
    }
    true
}

fn classical_m_bit(u: U256, v: U256) -> bool {
    if !u.bit(0) {
        false
    } else if !v.bit(0) {
        true
    } else {
        u > v
    }
}

fn reachable_state(v0: U256, n_steps: usize) -> ClassicalState {
    let mut st = ClassicalState { u: SECP256K1_P, v: v0, r: U256::ZERO, s: U256::from(1) };
    for _ in 0..n_steps {
        let ok = classical_step(&mut st);
        assert!(ok);
    }
    st
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bulk_prefix3_step_matches_generic_on_reachable_states() {
        let deadline = two_min_deadline();
        let generic = [build_generic_step(0), build_generic_step(1), build_generic_step(2)];
        let special = [build_special_step(0), build_special_step(1), build_special_step(2)];
        let mut sampler = Sampler::new(b"kaliski-equiv-sampler-v1", SECP256K1_P);

        for sample_idx in 0..512usize {
            if (sample_idx & 63) == 0 { check_deadline(deadline, "kaliski_equiv::bulk_prefix3_step_matches_generic_on_reachable_states"); }
            let v0 = sampler.next();
            for iter_idx in 0..3 {
                let st = reachable_state(v0, iter_idx);
                let g = run_step_circuit(&generic[iter_idx], st.u, st.v, st.r, st.s, false, true);
                let s = run_step_circuit(&special[iter_idx], st.u, st.v, st.r, st.s, false, true);
                assert_eq!(g, s, "specialized step mismatch at sample {} iter {}\nstate_before={:?}\ngeneric={:?}\nspecial={:?}", sample_idx, iter_idx, st, g, s);
            }
        }
    }

    #[test]
    fn bulk_prefix3_step_matches_classical_transition() {
        let deadline = two_min_deadline();
        let generic = [build_generic_step(0), build_generic_step(1), build_generic_step(2)];
        let mut sampler = Sampler::new(b"kaliski-equiv-sampler-v2", SECP256K1_P);

        for sample_idx in 0..512usize {
            if (sample_idx & 63) == 0 { check_deadline(deadline, "kaliski_equiv::bulk_prefix3_step_matches_classical_transition"); }
            let v0 = sampler.next();
            for iter_idx in 0..3 {
                let st = reachable_state(v0, iter_idx);
                let mut exp = st.clone();
                classical_step(&mut exp);
                let out = run_step_circuit(&generic[iter_idx], st.u, st.v, st.r, st.s, false, true);
                assert_eq!(out.0, exp.u, "u mismatch at sample {} iter {}", sample_idx, iter_idx);
                assert_eq!(out.1, exp.v, "v mismatch at sample {} iter {}", sample_idx, iter_idx);
                assert_eq!(out.2, exp.r, "r mismatch at sample {} iter {}", sample_idx, iter_idx);
                assert_eq!(out.3, exp.s, "s mismatch at sample {} iter {}", sample_idx, iter_idx);
                assert_eq!(out.4, classical_m_bit(st.u, st.v), "m mismatch at sample {} iter {}", sample_idx, iter_idx);
                assert!(out.5, "f should remain 1 through first 3 iterations; sample {} iter {}", sample_idx, iter_idx);
            }
        }
    }

    #[test]
    fn bulk_prefix3_three_step_sequence_matches_generic_forward_state() {
        let deadline = two_min_deadline();
        let mut sampler = Sampler::new(b"kaliski-equiv-sampler-v3", SECP256K1_P);

        let mut b = B::new();
        let u = b.alloc_qubits(N);
        let v = b.alloc_qubits(N);
        let r = b.alloc_qubits(N);
        let s = b.alloc_qubits(N);
        let f = b.alloc_qubit();
        let m_hist = b.alloc_qubits(3);
        for i in 0..3 {
            kaliski_iteration(&mut b, SECP256K1_P, &u, &v, &r, &s, m_hist[i], f, i);
        }
        let generic3 = StepCircuit {
            ops: b.ops,
            num_qubits: b.next_qubit as usize,
            num_bits: b.next_bit as usize,
            u, v, r, s, m: m_hist[2], f,
        };
        let special3 = build_special_three_steps();

        for sample_idx in 0..512usize {
            if (sample_idx & 63) == 0 { check_deadline(deadline, "kaliski_equiv::bulk_prefix3_three_step_sequence_matches_generic_forward_state"); }
            let v0 = sampler.next();
            let g = run_step_circuit(&generic3, SECP256K1_P, v0, U256::ZERO, U256::from(1), false, true);
            let s = run_step_circuit(&special3, SECP256K1_P, v0, U256::ZERO, U256::from(1), false, true);
            assert_eq!(g, s, "three-step mismatch at sample {}\ngeneric={:?}\nspecial={:?}", sample_idx, g, s);
        }
    }

    #[test]
    fn bulk_prefix3_backward_matches_generic_three_step_inverse() {
        let deadline = two_min_deadline();
        let mut sampler = Sampler::new(b"kaliski-equiv-sampler-v4", SECP256K1_P);

        let mut bg = B::new();
        let ug = bg.alloc_qubits(N);
        let vg = bg.alloc_qubits(N);
        let rg = bg.alloc_qubits(N);
        let sg = bg.alloc_qubits(N);
        let fg = bg.alloc_qubit();
        let mhg = bg.alloc_qubits(3);
        for i in 0..3 {
            kaliski_iteration(&mut bg, SECP256K1_P, &ug, &vg, &rg, &sg, mhg[i], fg, i);
        }
        for i in (0..3).rev() {
            super::super::kaliski_iteration_backward(&mut bg, SECP256K1_P, &ug, &vg, &rg, &sg, mhg[i], fg, i);
        }
        let generic_fb = StepCircuit {
            ops: bg.ops,
            num_qubits: bg.next_qubit as usize,
            num_bits: bg.next_bit as usize,
            u: ug, v: vg, r: rg, s: sg, m: mhg[2], f: fg,
        };

        let mut bs = B::new();
        let us = bs.alloc_qubits(N);
        let vs = bs.alloc_qubits(N);
        let rs = bs.alloc_qubits(N);
        let ss = bs.alloc_qubits(N);
        let fs = bs.alloc_qubit();
        let mhs = bs.alloc_qubits(3);
        for i in 0..3 {
            kaliski_iteration_bulk_prefix3(&mut bs, &us, &vs, &rs, &ss, mhs[i], i);
        }
        for i in (0..3).rev() {
            super::super::kaliski_iteration_bulk_prefix3_backward(&mut bs, &us, &vs, &rs, &ss, mhs[i], i);
        }
        let special_fb = StepCircuit {
            ops: bs.ops,
            num_qubits: bs.next_qubit as usize,
            num_bits: bs.next_bit as usize,
            u: us, v: vs, r: rs, s: ss, m: mhs[2], f: fs,
        };

        for sample_idx in 0..256usize {
            if (sample_idx & 63) == 0 { check_deadline(deadline, "kaliski_equiv::bulk_prefix3_backward_matches_generic_three_step_inverse"); }
            let v0 = sampler.next();
            let g = run_step_circuit(&generic_fb, SECP256K1_P, v0, U256::ZERO, U256::from(1), false, true);
            let s = run_step_circuit(&special_fb, SECP256K1_P, v0, U256::ZERO, U256::from(1), false, true);
            assert_eq!(g, s, "forward+backward mismatch at sample {}\ngeneric={:?}\nspecial={:?}", sample_idx, g, s);
        }
    }
}
