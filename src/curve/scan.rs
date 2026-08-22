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

/// The curvature scores for one piece, with the position the piece started at.
///
/// `start` is the 1-based position of the piece's first base in the original record, so
/// callers can map scores back onto the record without consulting the piece again.
#[derive(Debug, Clone)]
pub struct PieceCurves {
    pub start: usize,
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

/// Score a single piece.
pub fn score_piece(piece: &RecordPiece, params: &CurveParams) -> PieceCurves {
    // `bases` borrows straight out of the shared record, so no copy is made here.
    let curves: Vec<f64> = CurveIter::new(
        piece.bases().iter().copied(),
        params.roll_type,
        params.step_b,
        params.step_c,
        params.curve_scale,
    )
    .collect();
    PieceCurves {
        start: usize::from(piece.start),
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
