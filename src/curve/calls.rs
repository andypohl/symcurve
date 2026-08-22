//! Nucleosome calling from symmetry scores, and greedy selection of non-overlapping calls.

use std::collections::BTreeSet;

/// Half a nucleosome's footprint: a call spans `dyad - 73 ..= dyad + 73`.
pub const NUCLEOSOME_HALF_WIDTH: usize = 73;

/// A nucleosome's footprint in bases, `2 * NUCLEOSOME_HALF_WIDTH + 1`.
pub const NUCLEOSOME_SIZE: usize = 2 * NUCLEOSOME_HALF_WIDTH + 1;

/// Scores at or above this are reported clamped to it.
///
/// The symmetry stage substitutes this value where a dyad's symmetry component came out
/// exactly zero, which is a saturated score rather than a measured one.
pub const SATURATED_SCORE: f64 = 100.0;

/// The parameters of nucleosome calling.
#[derive(Debug, Clone, Copy)]
pub struct CallParams {
    /// Half the called footprint, on each side of the dyad.
    pub half_width: usize,
    /// Minimum gap required between two accepted calls, beyond the footprint itself.
    pub spacer: usize,
}

impl Default for CallParams {
    fn default() -> Self {
        Self {
            half_width: NUCLEOSOME_HALF_WIDTH,
            spacer: 30,
        }
    }
}

impl CallParams {
    /// How far apart two accepted dyads must be.
    ///
    /// The reference adds the footprint and the spacer, so with the defaults a call
    /// excludes anything within 177 bases of it on either side.
    pub fn exclusion(&self) -> usize {
        self.spacer + 2 * self.half_width + 1
    }
}

/// One called nucleosome, positioned by its dyad.
///
/// `dyad` is a zero-based index into the record, matching how the reference indexes its
/// arrays. The call covers `dyad - half_width ..= dyad + half_width`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NucleosomeCall {
    pub dyad: usize,
    pub score: f64,
}

impl NucleosomeCall {
    /// The score as reported, clamped at the saturated value.
    pub fn reported_score(&self) -> f64 {
        self.score.min(SATURATED_SCORE)
    }
}

/// Turn symmetry scores into overlapping nucleosome calls.
///
/// `scores` are `(dyad, score)` pairs with zero-based dyads, in ascending order. A dyad is
/// called when its score is above zero and its whole footprint fits inside the record.
///
/// Following the reference, the bounds are strict on both sides: `dyad - half_width` must
/// be greater than zero, not merely non-negative, and `dyad + half_width` must be less
/// than the record length rather than within it. That drops one otherwise-callable dyad at
/// each end, and is reproduced so the positions match.
pub fn call_nucleosomes(
    scores: &[(usize, f64)],
    record_len: usize,
    params: &CallParams,
) -> Vec<NucleosomeCall> {
    scores
        .iter()
        .filter(|(dyad, score)| {
            *score > 0.0 && *dyad > params.half_width && dyad + params.half_width < record_len
        })
        .map(|&(dyad, score)| NucleosomeCall { dyad, score })
        .collect()
}

/// Select non-overlapping calls, taking the highest scoring first.
///
/// Candidates are considered in descending score order, and one is accepted only if no
/// already-accepted dyad lies within [`CallParams::exclusion`] of it. The result is
/// returned in ascending dyad order, as the reference prints it.
///
/// The reference scans every accepted call for every candidate, which is quadratic and
/// dominates its runtime at chromosome scale. Accepted dyads are kept in a sorted set here
/// and the exclusion window is a range query over it, which is the same rule in
/// `O(n log n)`.
///
/// Two behaviours of the reference are reproduced deliberately. Its accepted-position
/// array is initialised holding a single zero, and its scan includes that element, so a
/// phantom call at position 0 rejects every candidate within the exclusion window of it;
/// no dyad at or below `exclusion` can ever be accepted. And where scores tie, the
/// reference's order comes from Perl hash iteration, which is randomised per process, so
/// its choice among equal scores is not reproducible even against itself; ties are broken
/// here by ascending dyad so that this implementation at least is deterministic.
pub fn greedy_non_overlapping(
    calls: &[NucleosomeCall],
    params: &CallParams,
) -> Vec<NucleosomeCall> {
    let exclusion = params.exclusion();
    let mut order: Vec<&NucleosomeCall> = calls.iter().collect();
    order.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.dyad.cmp(&b.dyad))
    });

    // Seeded with the reference's phantom position 0.
    let mut accepted: BTreeSet<usize> = BTreeSet::from([0]);
    let mut chosen: Vec<NucleosomeCall> = Vec::new();

    for call in order {
        let low = call.dyad.saturating_sub(exclusion);
        let high = call.dyad.saturating_add(exclusion);
        if accepted.range(low..=high).next().is_none() {
            accepted.insert(call.dyad);
            chosen.push(*call);
        }
    }

    chosen.sort_by_key(|c| c.dyad);
    chosen
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A direct transcription of the reference's GREEDYPOS, kept quadratic so it reads
    /// against the Perl and can serve as an oracle for the fast version.
    ///
    /// ```perl
    /// my @position = 0;
    /// foreach $dyad (sort { $calls{$b} <=> $calls{$a} } keys %calls) {
    ///     $index = 0;
    ///     for (my $i = 0; $i <= $no; $i++) {
    ///         if (($pos >= $position[$i]-$spacer-$size) and ($pos <= $position[$i]+$spacer+$size)) { $index = 1; }
    ///     }
    ///     if ($index == 0) { $no++; $position[$no] = $pos; }
    /// }
    /// ```
    fn perl_greedy(calls: &[NucleosomeCall], params: &CallParams) -> Vec<NucleosomeCall> {
        let exclusion = params.exclusion() as i64;
        let mut order: Vec<&NucleosomeCall> = calls.iter().collect();
        order.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.dyad.cmp(&b.dyad))
        });
        let mut position: Vec<i64> = vec![0]; // `my @position = 0;`
        let mut chosen = Vec::new();
        for call in order {
            let pos = call.dyad as i64;
            let mut index = false;
            for &p in &position {
                if pos >= p - exclusion && pos <= p + exclusion {
                    index = true;
                }
            }
            if !index {
                position.push(pos);
                chosen.push(*call);
            }
        }
        chosen.sort_by_key(|c| c.dyad);
        chosen
    }

    fn synthetic_calls(n: usize, span: usize, seed: u64) -> Vec<NucleosomeCall> {
        let mut x = seed;
        let mut out: Vec<NucleosomeCall> = (0..n)
            .map(|i| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                NucleosomeCall {
                    dyad: (x as usize) % span,
                    // Deliberately coarse so scores tie often.
                    score: ((x >> 20) % 50) as f64 / 10.0 + 0.1 + i as f64 * 0.0,
                }
            })
            .collect();
        out.sort_by_key(|c| c.dyad);
        out.dedup_by_key(|c| c.dyad);
        out
    }

    #[test]
    fn test_greedy_matches_the_reference_transcription() {
        for (n, span, seed) in [
            (50usize, 2000usize, 1u64),
            (200, 10_000, 2),
            (500, 20_000, 3),
            (1000, 5_000, 4), // dense: most candidates rejected
            (20, 100, 5),     // everything inside the phantom's exclusion
        ] {
            let calls = synthetic_calls(n, span, seed);
            for spacer in [0usize, 30, 200] {
                let params = CallParams {
                    half_width: NUCLEOSOME_HALF_WIDTH,
                    spacer,
                };
                let fast = greedy_non_overlapping(&calls, &params);
                let slow = perl_greedy(&calls, &params);
                assert_eq!(fast, slow, "n={n} span={span} spacer={spacer}");
            }
        }
    }

    #[test]
    fn test_accepted_calls_are_far_enough_apart() {
        let params = CallParams::default();
        let calls = synthetic_calls(800, 40_000, 9);
        let chosen = greedy_non_overlapping(&calls, &params);
        assert!(chosen.len() > 10);
        for pair in chosen.windows(2) {
            assert!(
                pair[1].dyad - pair[0].dyad > params.exclusion(),
                "{} and {} are too close",
                pair[0].dyad,
                pair[1].dyad
            );
        }
    }

    #[test]
    fn test_phantom_position_zero_blocks_the_start() {
        // The reference's accepted array starts holding a zero and its scan includes it,
        // so nothing within the exclusion window of position 0 can be accepted.
        let params = CallParams::default();
        let exclusion = params.exclusion();
        let calls = vec![
            NucleosomeCall {
                dyad: 100,
                score: 99.0,
            }, // inside the phantom's window
            NucleosomeCall {
                dyad: exclusion,
                score: 98.0,
            }, // exactly on the boundary
            NucleosomeCall {
                dyad: exclusion + 1,
                score: 1.0,
            }, // just outside
        ];
        let chosen = greedy_non_overlapping(&calls, &params);
        assert_eq!(chosen.len(), 1);
        assert_eq!(chosen[0].dyad, exclusion + 1);
        assert_eq!(chosen, perl_greedy(&calls, &params));
    }

    #[test]
    fn test_higher_scores_win() {
        let params = CallParams::default();
        let base = 10_000usize;
        let calls = vec![
            NucleosomeCall {
                dyad: base,
                score: 1.0,
            },
            NucleosomeCall {
                dyad: base + 10,
                score: 5.0,
            },
            NucleosomeCall {
                dyad: base + 20,
                score: 3.0,
            },
        ];
        let chosen = greedy_non_overlapping(&calls, &params);
        assert_eq!(chosen.len(), 1);
        assert_eq!(chosen[0].score, 5.0, "the strongest candidate should win");
    }

    #[test]
    fn test_calling_bounds_follow_the_reference() {
        let params = CallParams::default();
        let hw = params.half_width;
        let len = 1000usize;
        let scores = vec![
            (hw, 1.0),           // dyad - hw == 0, rejected: the test is strict
            (hw + 1, 1.0),       // first callable
            (len - hw - 1, 1.0), // last callable
            (len - hw, 1.0),     // dyad + hw == len, rejected
            (500, 0.0),          // zero score is never called
        ];
        let calls = call_nucleosomes(&scores, len, &params);
        let dyads: Vec<usize> = calls.iter().map(|c| c.dyad).collect();
        assert_eq!(dyads, vec![hw + 1, len - hw - 1]);
    }

    #[test]
    fn test_saturated_scores_are_reported_clamped() {
        let call = NucleosomeCall {
            dyad: 500,
            score: SATURATED_SCORE,
        };
        assert_eq!(call.reported_score(), 100.0);
        let call = NucleosomeCall {
            dyad: 500,
            score: 250.0,
        };
        assert_eq!(call.reported_score(), 100.0);
        let call = NucleosomeCall {
            dyad: 500,
            score: 0.25,
        };
        assert_eq!(call.reported_score(), 0.25);
    }
}
