//! Small targeted phase-map experiments for the specialized bulk-prefix step.
//!
//! Goal: detect whether the phase/nonphase pass pattern correlates with simple
//! arithmetic structure in the prefix length `k`.
//!
//! This stays tiny: no harness edits, no broad loops, just a few helper
//! classifiers that can guide further hard implementation work.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrefixPhaseRow {
    pub k: usize,
    pub strict_pass: bool,
}

pub const OBSERVED_ROWS: &[PrefixPhaseRow] = &[
    PrefixPhaseRow { k: 3, strict_pass: true },
    PrefixPhaseRow { k: 4, strict_pass: false },
    PrefixPhaseRow { k: 5, strict_pass: false },
    PrefixPhaseRow { k: 6, strict_pass: true },
    PrefixPhaseRow { k: 7, strict_pass: true },
    PrefixPhaseRow { k: 8, strict_pass: false },
    PrefixPhaseRow { k: 16, strict_pass: false },
    PrefixPhaseRow { k: 24, strict_pass: true },
    PrefixPhaseRow { k: 32, strict_pass: true },
    PrefixPhaseRow { k: 40, strict_pass: false },
    PrefixPhaseRow { k: 64, strict_pass: false },
    PrefixPhaseRow { k: 72, strict_pass: true },
    PrefixPhaseRow { k: 80, strict_pass: false },
    PrefixPhaseRow { k: 96, strict_pass: true },
    PrefixPhaseRow { k: 100, strict_pass: false },
    PrefixPhaseRow { k: 104, strict_pass: false },
    PrefixPhaseRow { k: 108, strict_pass: false },
    PrefixPhaseRow { k: 112, strict_pass: false },
    PrefixPhaseRow { k: 120, strict_pass: false },
    PrefixPhaseRow { k: 124, strict_pass: false },
    PrefixPhaseRow { k: 128, strict_pass: false },
];

fn v2(mut x: usize) -> usize {
    let mut c = 0;
    while x > 0 && (x & 1) == 0 {
        x >>= 1;
        c += 1;
    }
    c
}

fn is_multiple_of_24(x: usize) -> bool { x % 24 == 0 }
fn is_multiple_of_8(x: usize) -> bool { x % 8 == 0 }
fn is_multiple_of_16(x: usize) -> bool { x % 16 == 0 }
fn is_multiple_of_32(x: usize) -> bool { x % 32 == 0 }
fn is_multiple_of_3(x: usize) -> bool { x % 3 == 0 }
fn is_multiple_of_6(x: usize) -> bool { x % 6 == 0 }
fn is_multiple_of_12(x: usize) -> bool { x % 12 == 0 }

#[derive(Debug, Clone)]
pub struct PhaseMapSummary {
    pub passing: Vec<usize>,
    pub failing: Vec<usize>,
    pub passing_v2: Vec<usize>,
    pub failing_v2: Vec<usize>,
    pub passing_mod8: Vec<usize>,
    pub failing_mod8: Vec<usize>,
    pub passing_mult24: Vec<usize>,
    pub failing_mult24: Vec<usize>,
}

pub fn summarize_phase_map() -> PhaseMapSummary {
    let passing: Vec<usize> = OBSERVED_ROWS.iter().filter(|r| r.strict_pass).map(|r| r.k).collect();
    let failing: Vec<usize> = OBSERVED_ROWS.iter().filter(|r| !r.strict_pass).map(|r| r.k).collect();
    PhaseMapSummary {
        passing_v2: passing.iter().map(|&k| v2(k)).collect(),
        failing_v2: failing.iter().map(|&k| v2(k)).collect(),
        passing_mod8: passing.iter().map(|&k| k % 8).collect(),
        failing_mod8: failing.iter().map(|&k| k % 8).collect(),
        passing_mult24: passing.iter().copied().filter(|&k| is_multiple_of_24(k)).collect(),
        failing_mult24: failing.iter().copied().filter(|&k| is_multiple_of_24(k)).collect(),
        passing,
        failing,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inspect_phase_map_patterns() {
        let s = summarize_phase_map();
        eprintln!("=== bulk-prefix phase map patterns ===");
        eprintln!("passing ks        : {:?}", s.passing);
        eprintln!("failing ks        : {:?}", s.failing);
        eprintln!("passing v2        : {:?}", s.passing_v2);
        eprintln!("failing v2        : {:?}", s.failing_v2);
        eprintln!("passing k mod 8   : {:?}", s.passing_mod8);
        eprintln!("failing k mod 8   : {:?}", s.failing_mod8);
        eprintln!("passing mult of24 : {:?}", s.passing_mult24);
        eprintln!("failing mult of24 : {:?}", s.failing_mult24);
        eprintln!("helper flags:");
        for r in OBSERVED_ROWS {
            eprintln!(
                "k={:<3} pass={} v2={} mod8={} m3={} m6={} m12={} m16={} m24={} m32={} m8={}",
                r.k,
                r.strict_pass,
                v2(r.k),
                r.k % 8,
                is_multiple_of_3(r.k),
                is_multiple_of_6(r.k),
                is_multiple_of_12(r.k),
                is_multiple_of_16(r.k),
                is_multiple_of_24(r.k),
                is_multiple_of_32(r.k),
                is_multiple_of_8(r.k),
            );
        }
        eprintln!("=====================================");
        assert!(s.passing.contains(&96));
        assert!(s.failing.contains(&128));
    }
}
