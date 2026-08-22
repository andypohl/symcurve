//! This module contains some constants/matrices for curvature calculation.
use std::fmt;

/// The number of nucleotides in a triplet, which is also the number of dimensions in the
/// nucleotide matrices used for triplet -> value lookup.
pub const TRIPLET_SIZE: usize = 3;

/// A type alias for a 3D matrix sized 4x4x4 of f64 values. The first dimension is the
/// first nucleotide in a triplet, the second dimension is the second nucleotide in a triplet,
/// and the third dimension is the third nucleotide in a triplet.
pub type NucMatrix = [[[f64; 4]; 4]; 4];

/// The TWIST matrix is used to calculate the twist angle in three nucleotides of DNA.
/// The values are all 0.598647428 for all combinations of nucleotide triplets.
pub const TWIST: NucMatrix = [[[0.598647428; 4]; 4]; 4];

/// The TILT matrix is not really used in the current implementation, but is included here
/// for completeness. The values are all 0.0 for all combinations of nucleotide triplets.
pub const TILT: NucMatrix = [[[0.0; 4]; 4]; 4];

/// The "activated" version of the ROLL matrix is used to calculate the roll angle in three
/// nucleotides of DNA. This matrix differs from the simple matrix in that the angle
/// values represent a more activated state of the nucleosomes bound to the DNA.
///
/// Selected by [`RollType::Active`].
pub const ROLL_ACTIVE: NucMatrix = [
    [
        [0.0633, 0.3500, 4.6709, 2.64115],
        [6.2734, 0.3500, 7.7171, 4.44325],
        [4.8884, 3.9232, 5.0523, 6.8829],
        [5.4903, 3.9232, 5.3055, 5.3055],
    ],
    [
        [4.6709, 6.2734, 5.00295, 5.0673],
        [4.6709, 0.0633, 4.7618, 4.0633],
        [7.7000, 5.4903, 3.05865, 6.75525],
        [7.7000, 4.8884, 7.07195, 4.9907],
    ],
    [
        [4.0633, 4.44325, 5.9806, 5.51645],
        [5.0673, 2.64115, 6.62555, 5.51645],
        [4.9907, 5.3055, 5.89135, 9.0823],
        [6.75525, 6.8829, 5.89135, 9.0823],
    ],
    [
        [4.7618, 7.7171, 6.8996, 6.62555],
        [5.00295, 4.6709, 6.8996, 5.9806],
        [7.07195, 5.3055, 3.869, 5.9000],
        [3.05865, 5.0523, 3.869, 5.827],
    ],
];

/// The simple version of the ROLL matrix is used to calculate the roll angle in three nucleotides
/// of DNA.
///
/// Selected by [`RollType::Simple`].
pub const ROLL_SIMPLE: NucMatrix = [
    [
        [0.1, 0.0, 4.2, 1.6],
        [9.7, 0.0, 8.7, 3.6],
        [6.5, 2.0, 4.7, 6.3],
        [5.8, 2.0, 5.2, 5.2],
    ],
    [
        [7.3, 9.7, 7.8, 6.4],
        [7.3, 0.1, 6.2, 5.1],
        [10.0, 5.8, 0.7, 7.5],
        [10.0, 6.5, 5.8, 6.2],
    ],
    [
        [5.1, 3.6, 6.6, 5.6],
        [6.4, 1.6, 6.8, 5.6],
        [6.2, 5.2, 5.7, 8.2],
        [7.5, 6.3, 4.3, 8.2],
    ],
    [
        [6.2, 8.7, 9.6, 6.8],
        [7.8, 4.2, 9.6, 6.6],
        [5.8, 5.2, 3.0, 4.3],
        [0.7, 4.7, 3.0, 5.7],
    ],
];

/// Why a matrix lookup could not be performed.
///
/// An enum rather than a message string: the two cases are distinct, callers can tell
/// them apart, and neither needs an allocation to report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixLookupError {
    /// The slice handed in was not exactly three bases long.
    WrongLength(usize),
    /// A base was not A, C, G or T in either case.
    UnknownBase(u8),
}

impl fmt::Display for MatrixLookupError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongLength(n) => write!(f, "triplet must be of length 3, got {n}"),
            Self::UnknownBase(base) => {
                write!(f, "unrecognized nucleotide {:?}", *base as char)
            }
        }
    }
}

impl std::error::Error for MatrixLookupError {}

/// Which ROLL matrix a curvature calculation should use.
#[derive(Debug, Clone, Copy)]
pub enum RollType {
    Simple,
    Active,
}

/// Sentinel stored in `NUC_TABLE` for any byte that is not a scoreable base.
const NOT_A_BASE: u8 = u8::MAX;

/// Direct lookup from an ASCII byte to its index in a `NucMatrix`.
///
/// Both cases are populated, so soft-masked sequence (RepeatMasker lowercases
/// repetitive regions) is treated as ordinary sequence rather than as unknown bases.
/// Everything else, including `N` and the IUPAC ambiguity codes, maps to `NOT_A_BASE`;
/// callers are expected to have split those out already.
///
/// This runs three times per base over an entire genome, so it is a flat table rather
/// than a case conversion followed by a match.
const NUC_TABLE: [u8; 256] = {
    let mut table = [NOT_A_BASE; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'T' as usize] = 1;
    table[b't' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'C' as usize] = 3;
    table[b'c' as usize] = 3;
    table
};

/// Maps a nucleotide to its index in a `NucMatrix`, in either case.
///
/// Returns `None` for anything that is not A, C, G, or T.
#[inline]
pub fn nuc_index(base: u8) -> Option<usize> {
    match NUC_TABLE[base as usize] {
        NOT_A_BASE => None,
        ix => Some(ix as usize),
    }
}

/// Decodes a triplet into its three `NucMatrix` indices.
///
/// Several matrices are consulted for the same triplet on every base, so the ASCII
/// decode is done once here and the resulting indices reused, rather than re-decoding
/// per matrix. The result is a fixed-size array, so no allocation is involved.
///
/// # Errors
///
/// Returns a `MatrixLookupError` naming the offending base if any of the three is not
/// A, C, G, or T.
pub fn triplet_indices(
    triplet: &[u8; TRIPLET_SIZE],
) -> Result<[usize; TRIPLET_SIZE], MatrixLookupError> {
    let mut ixs = [0usize; TRIPLET_SIZE];
    for (slot, &base) in ixs.iter_mut().zip(triplet.iter()) {
        *slot = nuc_index(base).ok_or(MatrixLookupError::UnknownBase(base))?;
    }
    Ok(ixs)
}

/// Reads a value out of a matrix using indices already decoded by `triplet_indices`.
#[inline]
pub fn lookup_by_index(ixs: &[usize; TRIPLET_SIZE], matrix: &NucMatrix) -> f64 {
    matrix[ixs[0]][ixs[1]][ixs[2]]
}

/// Looks up a value in a nucleotide matrix based on a triplet of nucleotides.
///
/// This is the convenience form taking an arbitrary slice. The hot path decodes once with
/// `triplet_indices` and then calls `lookup_by_index` per matrix instead.
///
/// # Arguments
///
/// * `triplet` - A slice of u8 representing a triplet of nucleotides. Each u8 should be the ASCII
///   value of 'A', 'C', 'G', or 'T', upper or lower case.
/// * `matrix` - A reference to a `NucMatrix` to look up the value in.
///
/// # Returns
///
/// If the triplet is of length 3 and every base is recognized, this function returns a `Result`
/// containing the value at the corresponding position in the matrix. Otherwise it returns a
/// `Result` containing a `MatrixLookupError`.
///
/// # Errors
///
/// Returns a `MatrixLookupError` if the triplet is not of length 3, or if it contains a base
/// that is not A, C, G, or T. These are reported as distinct errors rather than being conflated.
pub fn matrix_lookup(triplet: &[u8], matrix: &NucMatrix) -> Result<f64, MatrixLookupError> {
    let triplet: &[u8; TRIPLET_SIZE] = triplet
        .try_into()
        .map_err(|_| MatrixLookupError::WrongLength(triplet.len()))?;
    Ok(lookup_by_index(&triplet_indices(triplet)?, matrix))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_spot_check_indexing() {
        assert_relative_eq!(TWIST[0][0][0], 0.598647428, epsilon = 1e-4);
        assert_relative_eq!(TWIST[1][1][1], 0.598647428, epsilon = 1e-4);
        assert_relative_eq!(ROLL_ACTIVE[1][2][0], 7.7, epsilon = 1e-4);
        assert_relative_eq!(
            matrix_lookup(b"AAA", &TWIST).unwrap(),
            0.598647428,
            epsilon = 1e-4
        );
        assert_relative_eq!(
            matrix_lookup(b"CCC", &TWIST).unwrap(),
            0.598647428,
            epsilon = 1e-4
        );
        assert_relative_eq!(
            matrix_lookup(b"CCA", &ROLL_SIMPLE).unwrap(),
            0.7,
            epsilon = 1e-4
        );
        assert!(matrix_lookup(b"AA", &ROLL_ACTIVE).is_err());
        assert!(matrix_lookup(b"AAAA", &ROLL_ACTIVE).is_err());
        assert!(matrix_lookup(b"AAN", &ROLL_ACTIVE).is_err());
    }

    #[test]
    fn test_nuc_index_over_the_whole_byte_range() {
        // The decode table is written out by hand, so check every possible byte
        // rather than a sample: a wrong or missing entry is otherwise easy to miss.
        for byte in 0u8..=255 {
            let expected = match byte {
                b'A' | b'a' => Some(0),
                b'T' | b't' => Some(1),
                b'G' | b'g' => Some(2),
                b'C' | b'c' => Some(3),
                _ => None,
            };
            assert_eq!(
                nuc_index(byte),
                expected,
                "byte {byte:?} ({:?})",
                byte as char
            );
        }
    }

    #[test]
    fn test_triplet_indices_matches_per_base_decoding() {
        let triplet = *b"CgA";
        assert_eq!(triplet_indices(&triplet).unwrap(), [3, 2, 0]);
        // Decoding once must agree with indexing the matrix the long way.
        assert_relative_eq!(
            lookup_by_index(&triplet_indices(&triplet).unwrap(), &ROLL_SIMPLE),
            matrix_lookup(&triplet, &ROLL_SIMPLE).unwrap(),
            epsilon = 1e-12
        );
        assert!(triplet_indices(b"CgN").is_err());
    }

    #[test]
    fn test_matrix_lookup_error_display() {
        assert_eq!(
            MatrixLookupError::WrongLength(2).to_string(),
            "triplet must be of length 3, got 2"
        );
        assert_eq!(
            MatrixLookupError::UnknownBase(b'N').to_string(),
            "unrecognized nucleotide 'N'"
        );
    }

    #[test]
    fn test_matrix_lookup_error_is_a_std_error() {
        // So it composes with `?` into Box<dyn Error>, anyhow and the like.
        fn takes_error<E: std::error::Error>(_: E) {}
        takes_error(MatrixLookupError::UnknownBase(b'N'));
        let boxed: Box<dyn std::error::Error> = MatrixLookupError::WrongLength(4).into();
        assert!(boxed.to_string().contains("length 3"));
    }
}
