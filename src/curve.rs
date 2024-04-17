//! This module contains functions for calculation of DNA curvature.
use std::collections::VecDeque;
use std::f64::consts::PI;
use std::fmt;
use std::iter::Iterator;

const NUC_MATRIX_DIMENSIONS: usize = 3;
pub type NucMatrix = [[[f64; 4]; 4]; 4];
#[allow(dead_code)]
const TWIST: NucMatrix = [[[0.598647428; 4]; 4]; 4];
#[allow(dead_code)]
const TILT: NucMatrix = [[[0.0; 4]; 4]; 4];
#[allow(dead_code)]
const ROLL_NUC: NucMatrix = [
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
#[allow(dead_code)]
const ROLL_DNASE: NucMatrix = [
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
struct MatrixLookupError {
    details: String,
}

impl fmt::Display for MatrixLookupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Error: {}", self.details)
    }
}

#[allow(dead_code)]
/// Given a triplet of nucleotides, look up the value in the given matrix.
fn matrix_lookup(triplet: &[u8], matrix: &NucMatrix) -> Result<f64, MatrixLookupError> {
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

#[allow(dead_code)]
#[derive(Debug, Clone)]
enum RollType {
    Simple,
    Activated,
}

#[allow(dead_code)]
#[derive(Clone)]
struct TripletData {
    twist: f64,
    roll: f64,
    tilt: f64,
    dx: f64,
    dy: f64,
    roll_type: RollType,
}

struct TripletWindowsIter<I: Iterator> {
    base_buffer: VecDeque<u8>,
    inner: I,
    twist_sum: f64,
    roll_type: RollType,
}

#[allow(dead_code)]
struct CoordsData {
    triplet_data: Option<TripletData>,
    x: f64,
    y: f64,
}

struct CoordsIter<I: Iterator> {
    inner: I,
    tail: bool,
    prev_x_coord: f64,
    prev_y_coord: f64,
    prev_dx: f64,
    prev_dy: f64,
}

impl<I> Iterator for TripletWindowsIter<I>
where
    I: Iterator<Item = u8>,
{
    type Item = TripletData;

    fn next(&mut self) -> Option<Self::Item> {
        while self.base_buffer.len() < NUC_MATRIX_DIMENSIONS {
            if let Some(item) = self.inner.next() {
                self.base_buffer.push_back(item);
            } else if self.base_buffer.is_empty() {
                return None;
            } else {
                break;
            }
        }

        if self.base_buffer.len() >= NUC_MATRIX_DIMENSIONS {
            let triplet: Vec<u8> = self.base_buffer.iter().cloned().take(3).collect();
            let twist = matrix_lookup(&triplet, &TWIST).unwrap();
            let roll = match self.roll_type {
                RollType::Simple => matrix_lookup(&triplet, &ROLL_DNASE).unwrap(),
                RollType::Activated => matrix_lookup(&triplet, &ROLL_NUC).unwrap(),
            };
            let tilt = matrix_lookup(&triplet, &TILT).unwrap();
            self.twist_sum += twist;

            let window = TripletData {
                twist,
                roll,
                tilt,
                dx: (roll * self.twist_sum.sin()) + (tilt * (self.twist_sum + PI / 2.0).sin()),
                dy: (roll * self.twist_sum.cos()) + (tilt * (self.twist_sum + PI / 2.0).cos()),
                roll_type: self.roll_type.clone(),
            };
            self.base_buffer.pop_front();
            Some(window)
        } else {
            None
        }
    }
}

impl<I> Iterator for CoordsIter<I>
where
    I: Iterator<Item = TripletData>,
{
    type Item = CoordsData;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(triplet_data) = self.inner.next() {
            self.prev_dx = triplet_data.dx;
            self.prev_dy = triplet_data.dy;
            Some(self.create_coords_data(Some(triplet_data.to_owned())))
        } else if !self.tail {
            self.tail = true;
            Some(self.create_coords_data(None))
        } else {
            None
        }
    }
}

impl<I> CoordsIter<I>
where
    I: Iterator<Item = TripletData>,
{
    fn create_coords_data(&mut self, triplet_data: Option<TripletData>) -> CoordsData {
        let x_coord = self.prev_x_coord + self.prev_dx;
        let y_coord = self.prev_y_coord + self.prev_dy;
        self.prev_x_coord = x_coord;
        self.prev_y_coord = y_coord;
        CoordsData {
            triplet_data,
            x: x_coord,
            y: y_coord,
        }
    }
}

trait TripletWindowsIterator: Iterator<Item = u8> + Sized {
    fn triplet_windows(self, roll_type: RollType) -> TripletWindowsIter<Self> {
        TripletWindowsIter {
            base_buffer: VecDeque::new(),
            inner: self,
            twist_sum: 0.0,
            roll_type,
        }
    }
}

trait CoordsIterator: Iterator<Item = TripletData> + Sized {
    fn coords_iter(self) -> CoordsIter<Self> {
        CoordsIter {
            inner: self,
            tail: false,
            prev_x_coord: 0.0,
            prev_y_coord: 0.0,
            prev_dx: 0.0,
            prev_dy: 0.0,
        }
    }
}

impl<I: Iterator<Item = u8>> TripletWindowsIterator for I {}
impl<I: Iterator<Item = TripletData>> CoordsIterator for I {}

#[cfg(test)]
mod tests {
    extern crate approx;
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_spot_check_indexing() {
        assert_relative_eq!(TWIST[0][0][0], 0.598647428);
        assert_relative_eq!(TWIST[1][1][1], 0.598647428);
        assert_relative_eq!(ROLL_DNASE[1][2][0], 10.0);
        assert_relative_eq!(matrix_lookup(b"AAA", &TWIST).unwrap(), 0.598647428);
        assert_relative_eq!(matrix_lookup(b"CCC", &TWIST).unwrap(), 0.598647428);
        assert!(matrix_lookup(b"AA", &ROLL_DNASE).is_err());
        assert!(matrix_lookup(b"AAAA", &ROLL_DNASE).is_err());
        assert!(matrix_lookup(b"AAN", &ROLL_DNASE).is_err());
    }

    // might do something like this with twist sum
    #[test]
    fn test_cumsum() {
        let nums = vec![1.0, 2.0, 3.0, 4.0];
        let cumsum: Vec<f64> = nums
            .iter()
            .scan(0.0, |acc, &x| {
                *acc += x;
                Some(*acc)
            })
            .collect();
        assert_eq!(cumsum, vec![1.0, 3.0, 6.0, 10.0]);
    }

    #[test]
    fn test_triplet_iter() {
        let dna = b"ACGTACGT";
        let windows: Vec<TripletData> = dna
            .iter()
            .cloned()
            .triplet_windows(RollType::Simple)
            .collect();
        assert_eq!(windows.len(), 6);
        let windows: Vec<TripletData> = dna
            .iter()
            .cloned()
            .triplet_windows(RollType::Activated)
            .collect();
        assert_eq!(windows.len(), 6);
    }

    #[test]
    fn test_coords_iter() {
        let dna = b"ACGTACGT";
        let windows: Vec<CoordsData> = dna
            .iter()
            .cloned()
            .triplet_windows(RollType::Simple)
            .coords_iter()
            .collect();
        assert_eq!(windows.len(), 7);
    }
}
