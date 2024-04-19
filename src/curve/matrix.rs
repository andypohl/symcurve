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

/// The simple version fo the ROLL matrix is used to calculate the roll angle in three nucleotides
/// of DNA.
pub const ROLL_SIMPLE: NucMatrix = [
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

/// The "activated" version of the ROLL matrix is used to calculate the roll angle in three
/// nucleotides of DNA. This matrix differs from the simple matrix in that the angle
/// values represent a more activated state of the nucleosomes bound to the DNA.
pub const ROLL_ACTIVE: NucMatrix = [
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

#[derive(Debug)]
pub struct MatrixLookupError {
    details: String,
}

impl fmt::Display for MatrixLookupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Error: {}", self.details)
    }
}

#[derive(Debug, Clone)]
pub(crate) enum RollType {
    Simple,
    Active,
}

/// Looks up a value in a nucleotide matrix based on a triplet of nucleotides.
///
/// This function takes a triplet of nucleotides and a nucleotide matrix, and returns the value
/// at the corresponding position in the matrix. The triplet is expected to contain the ASCII
/// values of 'A', 'C', 'G', or 'T'.  
///
/// # Arguments
///
/// * `triplet` - A slice of u8 representing a triplet of nucleotides. Each u8 should be the ASCII
/// value of 'A', 'C', 'G', or 'T'.
/// * `matrix` - A reference to a `NucMatrix` to look up the value in.
///
/// # Returns
///
/// If the triplet is valid and of length 3, this function returns a `Result` containing the value
/// at the corresponding position in the matrix. If the triplet is not valid or not of length 3,
/// it returns a `Result` containing a `MatrixLookupError`.
///
/// # Errors
///
/// Returns a `MatrixLookupError` if the triplet is not of length 3.  An unrecognized nucleotide
/// will cause this error because the triplet will not be of length 3.
pub(crate) fn matrix_lookup(triplet: &[u8], matrix: &NucMatrix) -> Result<f64, MatrixLookupError> {
    let ixs: Vec<usize> = triplet
        .iter()
        .map(|&x| match x {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        })
        .flatten()
        .collect();
    if ixs.len() != 3 {
        return Err(MatrixLookupError {
            details: "triplet must be of length 3".to_string(),
        });
    }
    Ok(matrix[ixs[0]][ixs[1]][ixs[2]])
}

#[cfg(test)]
mod tests {
    extern crate approx;
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_spot_check_indexing() {
        assert_relative_eq!(TWIST[0][0][0], 0.598647428);
        assert_relative_eq!(TWIST[1][1][1], 0.598647428);
        assert_relative_eq!(ROLL_ACTIVE[1][2][0], 10.0);
        assert_relative_eq!(matrix_lookup(b"AAA", &TWIST).unwrap(), 0.598647428);
        assert_relative_eq!(matrix_lookup(b"CCC", &TWIST).unwrap(), 0.598647428);
        assert!(matrix_lookup(b"AA", &ROLL_ACTIVE).is_err());
        assert!(matrix_lookup(b"AAAA", &ROLL_ACTIVE).is_err());
        assert!(matrix_lookup(b"AAN", &ROLL_ACTIVE).is_err());
    }

    #[test]
    fn test_matrix_lookup_error_display() {
        let error = MatrixLookupError {
            details: "Test error details".to_string(),
        };
        assert_eq!(format!("{}", error), "Error: Test error details");
    }
}
