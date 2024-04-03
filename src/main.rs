fn main() {
    println!("Hello, world!");
}

#[cfg(test)]
mod tests {
    use noodles_fasta as fasta;

    #[test]
    fn test_read_fasta() {
        // two sequences
        let src = b">sq0\nACGT\n>sq1\nN\n";
        let mut reader = fasta::Reader::new(&src[..]);
        let first_rec = reader.records().next().unwrap().unwrap();
        let second_rec = reader.records().next().unwrap().unwrap();
        assert_eq!(first_rec.definition().name(), b"sq0");
        assert_eq!(second_rec.definition().name(), b"sq1");
        assert_eq!(first_rec.sequence().as_ref(), b"ACGT".to_vec());
    }
}
