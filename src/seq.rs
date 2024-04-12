struct NucIterator<I>
where
    I: Iterator<Item = char>,
{
    inner: I,
}

impl<I> Iterator for NucIterator<I>
where
    I: Iterator<Item = char>,
{
    type Item = i32;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(c) = self.inner.next() {
            match c {
                'A' => return Some(0),
                'T' => return Some(1),
                'G' => return Some(2),
                'C' => return Some(3),
                _ => return Some(0),
            }
        }
        None
    }
}

trait IndexConverter: Iterator<Item = char> + Sized {
    fn to_nuc_index(self) -> NucIterator<Self> {
        NucIterator { inner: self }
    }
}

impl<T> IndexConverter for T where T: Iterator<Item = char> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_char_to_index() {
        let ixs: Vec<i32> = "ATGC".chars().to_nuc_index().collect();
        assert_eq!(ixs, vec![0, 1, 2, 3]);
        let ixs: Vec<i32> = "ATNGC".chars().to_nuc_index().collect();
        assert_eq!(ixs, vec![0, 1, 0, 2, 3]);
    }
}
