//! Functions for working with FASTA files.

use std::sync::Arc;

use noodles_core::Position;
use noodles_fasta::Record;
use noodles_fasta::record::Sequence;

/// One Record will be split into multiple RecordPieces.
/// The original Record is kept as an Arc so that each of the
/// RecordPieces can share the same ownership. Arc rather than Rc because the
/// pieces are scored in parallel, and Rc is not Send.
pub struct RecordPiece {
    pub record: Arc<Record>,
    pub start: Position,
    pub end: Position,
}

impl RecordPiece {
    fn new(record: Arc<Record>, start: Position, end: Position) -> Self {
        Self { record, start, end }
    }

    /// Get the sequence of the RecordPiece by slicing into the original Record.
    ///
    /// This copies the slice out. Prefer [`RecordPiece::bases`] when the bytes are only
    /// going to be read, which is the case on the scoring path.
    pub fn sequence(&self) -> Sequence {
        self.record.sequence().slice(self.start..=self.end).unwrap()
    }

    /// Borrow this piece's bases directly out of the shared Record.
    ///
    /// Unlike [`RecordPiece::sequence`] this allocates nothing, so scoring a piece does
    /// not begin by copying it. `start` and `end` are 1-based and inclusive.
    pub fn bases(&self) -> &[u8] {
        let start = usize::from(self.start) - 1;
        let end = usize::from(self.end);
        &self.record.sequence().as_ref()[start..end]
    }
}

/// Returns true for any base that cannot be scored and must therefore break the sequence.
///
/// A, C, G and T are scoreable in either case: RepeatMasker lowercases repetitive regions,
/// but soft-masking is an annotation rather than missing data, so `acgt` is ordinary
/// sequence. Everything else -- `N` and the IUPAC ambiguity codes (R, Y, S, W, K, M, B,
/// D, H, V) -- represents a base that is genuinely unknown.
fn is_gap(base: u8) -> bool {
    !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

#[allow(dead_code)]
/// Given a record, split the sequence at runs of unscoreable bases.
///
/// Returns a vector of pieces, each covering a stretch of sequence that contains only
/// A, C, G and T (in either case). Runs of `N` or IUPAC ambiguity codes are dropped, and
/// the sequence is split there. Each piece records its own start-end position in the
/// original record, the positions being 1-based.
///
/// Splitting matters because curvature is computed from a sliding window over a running
/// sum of coordinates. Merely deleting the unknown bases would let a window span the gap
/// and derive a value from bases that are far apart in the real sequence.
///
/// Input:
/// ```text
/// >chr42
/// ATGCATGC
/// NNNNATGC
/// A
/// ```
///
/// Output:
/// ```text
/// >chr42 1-8
/// ATGCATGC
/// >chr42 13-17
/// ATGCA
/// ```
pub fn split_seq_by_gaps(record: Record) -> Vec<RecordPiece> {
    // Move the record into a single Arc up front. Every piece then clones this
    // one handle, so they all point at the same allocation. Calling Arc::new
    // per piece would instead allocate a fresh box holding a full copy of the
    // sequence, which is what the shared ownership here is meant to avoid.
    let record = Arc::new(record);
    let mut records = Vec::new();
    let n = record.sequence().len();
    let seq = record.sequence().as_ref();
    let mut pos = 0;
    // classic two-pointer approach is tried-and-true
    // but might not be the most idiomatic Rust
    while pos < n {
        while (pos < n) && is_gap(seq[pos]) {
            pos += 1;
        }
        let left = pos;
        while (pos < n) && !is_gap(seq[pos]) {
            pos += 1;
        }
        let right = pos;
        if left < right {
            // Position is 1-based so add 1 to left
            let start = Position::try_from(left + 1).unwrap();
            let end = Position::try_from(right).unwrap();
            let piece = RecordPiece::new(Arc::clone(&record), start, end);
            records.push(piece);
        }
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::matrix;
    use approx::assert_relative_eq;

    #[test]
    fn test_read_fasta() {
        // two sequences
        let src = b">sq0\nACGT\n>sq1\nN\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let first_rec = reader.records().next().unwrap().unwrap();
        let second_rec = reader.records().next().unwrap().unwrap();
        assert_eq!(first_rec.name(), b"sq0");
        assert_eq!(second_rec.name(), b"sq1");
        let start = Position::try_from(2).unwrap();
        let end = Position::try_from(3).unwrap();
        assert_eq!(
            first_rec.sequence().slice(start..=end).unwrap().as_ref(),
            b"CG".to_vec()
        );
    }

    #[test]
    fn test_windows() {
        let seq = b"ACGTACGTACGTACGTACGT";
        let mut seq_bytes = seq.windows(3);
        assert_eq!(seq_bytes.next().unwrap(), b"ACG");
        assert_eq!(seq_bytes.next().unwrap(), b"CGT");
        assert_eq!(seq_bytes.next().unwrap(), b"GTA");
    }

    #[test]
    fn test_splitting() {
        let src = b">chr42\nATGCATGCNNNNATGCA\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let split_records: Vec<_> = reader
            .records()
            .flat_map(|rec| split_seq_by_gaps(rec.unwrap()))
            .collect();
        assert_eq!(split_records.len(), 2);
        assert_eq!(split_records[0].sequence().as_ref(), b"ATGCATGC".to_vec());
        assert_eq!(split_records[1].sequence().as_ref(), b"ATGCA".to_vec());
        assert_eq!(split_records[1].sequence().as_ref(), b"ATGCA".to_vec());
        assert_eq!(usize::from(split_records[1].start), 13);
        assert_eq!(usize::from(split_records[1].end), 17);
    }

    #[test]
    fn test_pieces_share_one_record() {
        // Four pieces separated by runs of Ns. All of them must point at the
        // same underlying Record rather than each holding their own copy.
        let src = b">chr42\nACGTNNACGTNNACGTNNACGT\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let record = reader.records().next().unwrap().unwrap();
        let pieces = split_seq_by_gaps(record);
        assert_eq!(pieces.len(), 4);
        // One allocation, one handle per piece.
        assert_eq!(Arc::strong_count(&pieces[0].record), 4);
        for piece in &pieces[1..] {
            assert!(Arc::ptr_eq(&pieces[0].record, &piece.record));
        }
        // Sharing the record must not disturb the slicing.
        assert_eq!(pieces[0].sequence().as_ref(), b"ACGT".to_vec());
        assert_eq!(pieces[3].sequence().as_ref(), b"ACGT".to_vec());
        assert_eq!(usize::from(pieces[3].start), 19);
        assert_eq!(usize::from(pieces[3].end), 22);
    }

    #[test]
    fn test_softmasked_bases_are_kept() {
        // Lowercase acgt is soft-masked repeat sequence, not missing data, so it must
        // flow through as ordinary sequence rather than splitting the record.
        let src = b">chr42\nACGTacgtACGT\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let record = reader.records().next().unwrap().unwrap();
        let pieces = split_seq_by_gaps(record);
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].sequence().as_ref(), b"ACGTacgtACGT".to_vec());
        assert_eq!(usize::from(pieces[0].start), 1);
        assert_eq!(usize::from(pieces[0].end), 12);
    }

    #[test]
    fn test_softmasked_bases_score_as_their_uppercase_form() {
        // The lookup must agree across cases, or soft-masked regions would yield
        // different curvature than the same sequence unmasked.
        for (upper, lower) in [(b"ACG", b"acg"), (b"TTT", b"ttt"), (b"CCA", b"cca")] {
            assert_relative_eq!(
                matrix::matrix_lookup(upper, &matrix::ROLL_SIMPLE).unwrap(),
                matrix::matrix_lookup(lower, &matrix::ROLL_SIMPLE).unwrap(),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn test_ambiguity_codes_split_like_n() {
        // R and Y are genuinely unknown bases. Before this they reached the matrix
        // lookup and panicked; now they gap-split exactly as N does.
        let src = b">chr42\nACGTRYACGT\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let record = reader.records().next().unwrap().unwrap();
        let pieces = split_seq_by_gaps(record);
        assert_eq!(pieces.len(), 2);
        assert_eq!(pieces[0].sequence().as_ref(), b"ACGT".to_vec());
        assert_eq!(pieces[1].sequence().as_ref(), b"ACGT".to_vec());
        assert_eq!(usize::from(pieces[1].start), 7);
        assert_eq!(usize::from(pieces[1].end), 10);
    }

    #[test]
    fn test_every_iupac_code_is_a_gap() {
        for &code in b"NRYSWKMBDHVnryswkmbdhv" {
            assert!(is_gap(code), "{:?} should be a gap", code as char);
        }
        for &code in b"ACGTacgt" {
            assert!(!is_gap(code), "{:?} should not be a gap", code as char);
        }
    }

    #[test]
    fn test_mixed_gaps_and_softmasking() {
        // A realistic shape: soft-masked repeat, an assembly gap, an ambiguity code.
        let src = b">chr42\nACGTacgtNNNNacgtRACGT\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let record = reader.records().next().unwrap().unwrap();
        let pieces = split_seq_by_gaps(record);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0].sequence().as_ref(), b"ACGTacgt".to_vec());
        assert_eq!(pieces[1].sequence().as_ref(), b"acgt".to_vec());
        assert_eq!(pieces[2].sequence().as_ref(), b"ACGT".to_vec());
    }

    #[test]
    fn test_lookup_rejects_unknown_base_distinctly() {
        // The two failure modes used to be conflated under "must be of length 3".
        let bad_base = matrix::matrix_lookup(b"AAN", &matrix::ROLL_SIMPLE).unwrap_err();
        assert!(
            bad_base.to_string().contains("unrecognized nucleotide"),
            "got: {bad_base}"
        );
        let bad_len = matrix::matrix_lookup(b"AA", &matrix::ROLL_SIMPLE).unwrap_err();
        assert!(bad_len.to_string().contains("length 3"), "got: {bad_len}");
    }

    #[test]
    fn test_splitting_empty() {
        let src = b">chr42\n\n";
        let mut reader = noodles_fasta::io::Reader::new(&src[..]);
        let split_records: Vec<_> = reader
            .records()
            .flat_map(|rec| split_seq_by_gaps(rec.unwrap()))
            .collect();
        assert_eq!(split_records.len(), 0);
    }
}
