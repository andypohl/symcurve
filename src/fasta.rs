//! Functions for working with FASTA files.

use noodles_core::Position;
use noodles_fasta::record::Definition;
use noodles_fasta::{self, Record};

#[allow(dead_code)]
/// Given a record, split the sequence by runs of Ns.
/// Returns a vector of records, each with a sequence that does not contain any Ns.
/// The description of each record is set to the start-end position of the sequence,
/// the positions being 1-based.
pub fn split_seq_by_n(record: Record) -> Vec<Record> {
    let mut records = Vec::new();
    let n = record.sequence().len();
    let seq = record.sequence().as_ref();
    let mut pos = 0;
    // classic two-pointer approach is tried-and-true
    // but might not be the most idiomatic Rust
    while pos < n {
        while (pos < n) && (seq[pos] == b'N') {
            pos += 1;
        }
        let left = pos;
        while (pos < n) && (seq[pos] != b'N') {
            pos += 1;
        }
        let right = pos;
        if left < right {
            // Position is 1-based so add 1 to left
            let start = Position::try_from(left + 1).unwrap();
            let end = Position::try_from(right).unwrap();
            let subseq = record.sequence().slice(start..=end).unwrap();
            let info = format!("{}-{}", start, end).as_bytes().to_vec();
            let def = Definition::new(record.name(), Some(info));
            let rec = Record::new(def, subseq);
            records.push(rec);
        }
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_fasta() {
        // two sequences
        let src = b">sq0\nACGT\n>sq1\nN\n";
        let mut reader = noodles_fasta::Reader::new(&src[..]);
        let first_rec = reader.records().next().unwrap().unwrap();
        let second_rec = reader.records().next().unwrap().unwrap();
        assert_eq!(first_rec.definition().name(), b"sq0");
        assert_eq!(second_rec.definition().name(), b"sq1");
        let start = Position::try_from(2).unwrap();
        let end = Position::try_from(3).unwrap();
        assert_eq!(
            first_rec.sequence().slice(start..=end).unwrap().as_ref(),
            b"CG".to_vec()
        );
    }

    #[test]
    fn test_windows() {
        let seq: &str = "ACGTACGTACGTACGTACGT";
        let mut seq_bytes = seq.as_bytes().windows(3);
        assert_eq!(seq_bytes.next().unwrap(), b"ACG");
        assert_eq!(seq_bytes.next().unwrap(), b"CGT");
        assert_eq!(seq_bytes.next().unwrap(), b"GTA");
    }

    #[test]
    fn test_splitting() {
        let src = b">chr42\nATGCATGCNNNNATGCA\n";
        let mut reader = noodles_fasta::Reader::new(&src[..]);
        let split_records: Vec<_> = reader
            .records()
            .flat_map(|rec| split_seq_by_n(rec.unwrap()))
            .collect();
        assert_eq!(split_records.len(), 2);
        //println!("{:?}", split_records[1]);
        assert_eq!(split_records[0].sequence().as_ref(), b"ATGCATGC".to_vec());
        assert_eq!(split_records[1].sequence().as_ref(), b"ATGCA".to_vec());
        // check the descriptions where the ranges are stashed
        assert_eq!(
            split_records[0].definition().description().unwrap(),
            b"1-8".to_vec()
        );
        assert_eq!(
            split_records[1].definition().description().unwrap(),
            b"13-17".to_vec()
        );
    }

    #[test]
    fn test_splitting_empty() {
        let src = b">chr42\n\n";
        let mut reader = noodles_fasta::Reader::new(&src[..]);
        let split_records: Vec<_> = reader
            .records()
            .flat_map(|rec| split_seq_by_n(rec.unwrap()))
            .collect();
        assert_eq!(split_records.len(), 0);
    }
}
