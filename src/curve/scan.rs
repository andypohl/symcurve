//! Scoring whole records, piece by piece, across threads.
//!
//! [`crate::fasta::split_seq_by_gaps`] breaks a record at unscoreable bases, and the
//! resulting pieces are independent: curvature is computed from a sliding window that
//! never spans a gap, so no piece needs anything from its neighbours. That makes the
//! pieces the natural unit of parallelism.

use rayon::prelude::*;

use crate::curve::iters::CurveIter;
use crate::curve::matrix::RollType;
use crate::fasta::RecordPiece;

/// The curvature scores for one piece, with the positions they belong to.
///
/// `start` is the 1-based position of the piece's first base in the original record.
/// `curve_start` is the 1-based position that `curves[0]` scores, which is further in by
/// the pipeline's lead-in: the windows need context on both sides, so the first several
/// bases of a piece get no score.
#[derive(Debug, Clone)]
pub struct PieceCurves {
    pub start: usize,
    pub curve_start: usize,
    pub curves: Vec<f64>,
}

/// The parameters that define a curvature calculation.
///
/// Grouped into one struct because they travel together through every layer, and a
/// function taking five bare numbers is easy to call wrongly.
#[derive(Debug, Clone, Copy)]
pub struct CurveParams {
    pub roll_type: RollType,
    /// Half the rolling-mean window, minus one: the window is `2 * step_b + 1`.
    pub step_b: usize,
    /// Distance from the midpoint base to each side of the curve window.
    pub step_c: usize,
    pub curve_scale: f64,
}

impl CurveParams {
    /// How many bases at each end of a piece receive no score.
    ///
    /// The rolling mean consumes `step_b` on each side and the distance window a further
    /// `step_c`, plus one for the triplet. A piece of length `n` therefore yields
    /// `n - 2 * lead_in()` scores, the first of which belongs to the base at offset
    /// `lead_in()`. This matches the reference implementation, which indexes curvature
    /// from `curvstep + stepone`.
    pub fn lead_in(&self) -> usize {
        self.step_b + self.step_c + 1
    }
}

/// The default number of scores one chunk produces.
///
/// Large enough that the per-chunk lead-in overhead is negligible, small enough that a
/// handful in flight stay in cache and bound memory. Overridden from the memory budget.
pub const DEFAULT_CHUNK_SCORES: usize = 1 << 20;

/// One unit of independently scoreable work within a piece.
///
/// Chunks can be scored with fresh state and their results concatenated, because
/// curvature is a distance between local averages: a chunk that starts partway into a
/// piece has the wrong starting coordinate and the wrong accumulated twist, but the first
/// is a translation of the traced path and the second a rotation of it, and neither
/// changes the distances the scores are made of. Each chunk therefore reads `lead_in`
/// bases of context beyond its output range at each end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chunk {
    /// Offset into the piece where this chunk starts reading.
    pub in_start: usize,
    /// Offset into the piece where this chunk stops reading, exclusive.
    pub in_end: usize,
    /// Index into the piece's scores of this chunk's first score.
    pub out_start: usize,
}

/// Divide a piece of `piece_len` bases into chunks producing `chunk_scores` scores each.
///
/// Returns no chunks when the piece is too short to produce any score at all.
pub fn chunks(piece_len: usize, params: &CurveParams, chunk_scores: usize) -> Vec<Chunk> {
    let lead = params.lead_in();
    let out_len = piece_len.saturating_sub(2 * lead);
    let chunk_scores = chunk_scores.max(1);
    (0..out_len)
        .step_by(chunk_scores)
        .map(|out_start| {
            let out_end = (out_start + chunk_scores).min(out_len);
            Chunk {
                // Scoring bases [p, q + 2*lead) yields exactly the scores [p, q).
                in_start: out_start,
                in_end: out_end + 2 * lead,
                out_start,
            }
        })
        .collect()
}

/// Score a run of bases with fresh state.
pub fn score_bases(bases: &[u8], params: &CurveParams) -> Vec<f64> {
    CurveIter::new(
        bases.iter().copied(),
        params.roll_type,
        params.step_b,
        params.step_c,
        params.curve_scale,
    )
    .collect()
}

/// Score a single piece.
pub fn score_piece(piece: &RecordPiece, params: &CurveParams) -> PieceCurves {
    score_piece_chunked(piece, params, DEFAULT_CHUNK_SCORES)
}

/// Score a single piece, dividing it into chunks of the given size.
///
/// A record is often one enormous piece, so splitting only by piece leaves a chromosome
/// on a single thread. Chunking within the piece both spreads that work and bounds how
/// much of it is in memory at once.
pub fn score_piece_chunked(
    piece: &RecordPiece,
    params: &CurveParams,
    chunk_scores: usize,
) -> PieceCurves {
    // `bases` borrows straight out of the shared record, so no copy is made here.
    let bases = piece.bases();
    // collect() on an indexed parallel iterator preserves order, so the per-chunk score
    // runs concatenate back into the piece's scores in position order.
    let per_chunk: Vec<Vec<f64>> = chunks(bases.len(), params, chunk_scores)
        .into_par_iter()
        .map(|chunk| score_bases(&bases[chunk.in_start..chunk.in_end], params))
        .collect();
    let curves: Vec<f64> = per_chunk.concat();
    let start = usize::from(piece.start);
    PieceCurves {
        start,
        curve_start: start + params.lead_in(),
        curves,
    }
}

/// Score every piece, in parallel, preserving input order.
///
/// Pieces vary enormously in length -- a record may hold one piece covering most of a
/// chromosome alongside many short ones -- so this splits by piece and lets rayon's work
/// stealing handle the imbalance, rather than dividing the work evenly up front.
pub fn score_pieces(pieces: &[RecordPiece], params: &CurveParams) -> Vec<PieceCurves> {
    pieces
        .par_iter()
        .map(|piece| score_piece(piece, params))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::split_seq_by_gaps;
    use approx::assert_relative_eq;

    fn params() -> CurveParams {
        CurveParams {
            roll_type: RollType::Simple,
            step_b: 5,
            step_c: 15,
            curve_scale: 0.33335,
        }
    }

    fn record_of(seq: &[u8]) -> noodles_fasta::Record {
        let mut src = b">chr42\n".to_vec();
        src.extend_from_slice(seq);
        src.push(b'\n');
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        reader.records().next().unwrap().unwrap()
    }

    /// Scoring pieces across threads requires them to be Send and Sync. This is what
    /// Arc buys over Rc, which is not Send; asserting it here states the requirement
    /// rather than leaving it implied by a par_iter call elsewhere.
    #[test]
    fn test_record_piece_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<RecordPiece>();
        assert_send_sync::<PieceCurves>();
        assert_send_sync::<CurveParams>();
    }

    #[test]
    fn test_parallel_matches_serial() {
        // Several pieces of differing lengths, separated by gaps.
        let unit = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let mut seq = Vec::new();
        for reps in [3usize, 1, 5, 2] {
            for _ in 0..reps {
                seq.extend_from_slice(unit);
            }
            seq.extend_from_slice(b"NNNN");
        }
        let pieces = split_seq_by_gaps(record_of(&seq));
        assert_eq!(pieces.len(), 4);

        let parallel = score_pieces(&pieces, &params());
        let serial: Vec<PieceCurves> = pieces.iter().map(|p| score_piece(p, &params())).collect();

        assert_eq!(parallel.len(), serial.len());
        for (par, ser) in parallel.iter().zip(&serial) {
            // Order must be preserved, not just the multiset of results.
            assert_eq!(par.start, ser.start);
            assert_eq!(par.curves.len(), ser.curves.len());
            for (a, b) in par.curves.iter().zip(&ser.curves) {
                assert_relative_eq!(a, b, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn test_piece_start_positions_are_reported() {
        let seq = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATCNNNNCCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let pieces = split_seq_by_gaps(record_of(seq));
        let scored = score_pieces(&pieces, &params());
        assert_eq!(scored.len(), 2);
        assert_eq!(scored[0].start, 1);
        assert_eq!(scored[1].start, 55);
    }

    #[test]
    fn test_lead_in_matches_the_reference_convention() {
        // The Perl reference scores indices curvstep+stepone .. len-curvstep-stepone,
        // with stepone = step_b + 1 and curvstep = step_c. Check both the count and the
        // position of the first score, since an off-by-one here shifts every value in
        // the output file relative to the genome.
        let p = params();
        assert_eq!(p.lead_in(), 21); // step_b 5 + step_c 15 + 1, i.e. stepone 6 + curvstep 15
        let unit = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let mut seq = Vec::new();
        for _ in 0..4 {
            seq.extend_from_slice(unit);
        }
        let n = seq.len();
        let pieces = split_seq_by_gaps(record_of(&seq));
        let scored = score_pieces(&pieces, &p);
        assert_eq!(scored.len(), 1);
        assert_eq!(scored[0].curves.len(), n - 2 * p.lead_in());
        assert_eq!(scored[0].start, 1);
        assert_eq!(scored[0].curve_start, 1 + p.lead_in());
    }

    #[test]
    fn test_curve_start_is_offset_from_the_piece_not_the_record() {
        let seq = b"NNNNNCCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATCCCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let pieces = split_seq_by_gaps(record_of(seq));
        assert_eq!(pieces.len(), 1);
        let scored = score_pieces(&pieces, &params());
        assert_eq!(scored[0].start, 6); // after the 5 leading Ns, 1-based
        assert_eq!(scored[0].curve_start, 6 + params().lead_in());
    }

    #[test]
    fn test_chunk_ranges_tile_the_output_exactly() {
        let p = params();
        let lead = p.lead_in();
        let piece_len = 1000usize;
        let out_len = piece_len - 2 * lead;
        let cs = 100usize;
        let cs_list = chunks(piece_len, &p, cs);
        assert_eq!(cs_list.len(), out_len.div_ceil(cs));
        // Chunks must tile the output range with no gap and no overlap.
        let mut expected_out = 0usize;
        for c in &cs_list {
            assert_eq!(c.out_start, expected_out);
            // Scoring bases [in_start, in_end) yields in_end - in_start - 2*lead scores.
            let produced = c.in_end - c.in_start - 2 * lead;
            expected_out += produced;
            assert!(c.in_end <= piece_len, "chunk reads past the piece");
        }
        assert_eq!(
            expected_out, out_len,
            "chunks do not cover the output exactly"
        );
    }

    #[test]
    fn test_short_piece_produces_no_chunks() {
        let p = params();
        for len in [0usize, 1, 2 * p.lead_in(), 2 * p.lead_in() + 1] {
            let got = chunks(len, &p, 100);
            let expect_scores = len.saturating_sub(2 * p.lead_in());
            assert_eq!(got.is_empty(), expect_scores == 0, "len {len}");
        }
    }

    #[test]
    fn test_chunked_scoring_matches_unchunked() {
        // The whole point of chunking is that it changes nothing about the answer.
        // Chunk sizes are chosen to straddle the boundaries: one chunk, exact multiples,
        // and sizes that leave a short final chunk.
        let unit = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let mut seq = Vec::new();
        for _ in 0..40 {
            seq.extend_from_slice(unit);
        }
        let pieces = split_seq_by_gaps(record_of(&seq));
        let piece = &pieces[0];
        let p = params();

        let reference = score_bases(piece.bases(), &p);
        assert!(reference.len() > 1000);

        for chunk_scores in [
            1usize,
            7,
            64,
            500,
            1000,
            reference.len(),
            reference.len() * 2,
        ] {
            let got = score_piece_chunked(piece, &p, chunk_scores);
            assert_eq!(
                got.curves.len(),
                reference.len(),
                "length differs at chunk_scores {chunk_scores}"
            );
            for (i, (a, b)) in reference.iter().zip(&got.curves).enumerate() {
                // Not bit-identical: a chunk accumulates twist over a shorter run, so the
                // rounding differs. The values are mathematically the same.
                assert_relative_eq!(a, b, epsilon = 1e-9, max_relative = 1e-9);
                let _ = i;
            }
        }
    }

    #[test]
    fn test_empty_input_is_not_an_error() {
        let scored = score_pieces(&[], &params());
        assert!(scored.is_empty());
    }

    #[test]
    fn test_piece_shorter_than_the_window_yields_nothing() {
        // A piece too short to fill the windows must produce no scores rather than panic.
        let pieces = split_seq_by_gaps(record_of(b"ACGTACGT"));
        assert_eq!(pieces.len(), 1);
        let scored = score_pieces(&pieces, &params());
        assert_eq!(scored.len(), 1);
        assert!(scored[0].curves.is_empty());
    }
}
