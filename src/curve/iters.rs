//! This module provides data structures and Iterator implementations for the calculation of DNA
//! curvature.
//!
//! It includes the necessary data structures for representing the DNA data and the traits and
//! implementations for iterating over this data. The iterators provided allow for efficient and
//! convenient traversal and manipulation of the DNA data for the purpose of curvature calculation.
use crate::curve::matrix;
use std::collections::VecDeque;
use std::f64::consts::PI;
use std::iter::Iterator;

/// Represents the data for a triplet of nucleotides.
///
/// # Layer 1
///
/// This struct contains the twist, roll, and tilt values for a triplet of nucleotides, as well as
/// the deltas `dx` and `dy` and the roll type. *`roll_type` may be removed from this struct in the
/// future to accommodate more-general matrix options.*
///
/// # Fields
///
/// * `twist`: The twist value for the triplet.
/// * `roll`: The roll value for the triplet.
/// * `tilt`: The tilt value for the triplet.
/// * `dx`: The delta x value, calculated based on the roll and tilt.
/// * `dy`: The delta y value, calculated based on the roll and tilt.
/// * `roll_type`: The type of roll (either simple or activated).
#[derive(Clone)]
struct TripletData {
    twist: f64,
    roll: f64,
    tilt: f64,
    dx: f64,
    dy: f64,
    roll_type: matrix::RollType,
}

/// An iterator-wrapping struct that yields TripletData from an inner `u8` iterator.
///
/// # Layer 1
///
/// `TripletWindowsIter` wraps around another iterator that yields `u8` (representing nucleotides),
/// and looks up the roll/tilt/twist values for each triplet of nucleotides in the inner iterator,
/// then populates a `TripletData` struct with these values.  While iterating, it also keeps track
/// of the sum of the twist values for the current window of triplets.
///
/// # Type Parameters
///
/// * `I`: The type of the inner iterator. Must be an iterator over `u8`.
///
/// # Fields
///
/// * `base_buffer`: A buffer that stores the current triplet of nucleotides.
/// * `inner`: The inner iterator that yields `u8`.
/// * `twist_sum`: The sum of the twist values for the current triplet.
/// * `roll_type`: The current roll type.
struct TripletWindowsIter<I: Iterator> {
    base_buffer: VecDeque<u8>,
    inner: I,
    twist_sum: f64,
    roll_type: matrix::RollType,
}

/// Represents the coordinates and associated data for a triplet of nucleotides.
///
/// # Layer 2
///
/// `CoordsData` contains the x and y coordinates calculated from the `TripletData`, as well as
/// the `TripletData` itself. The `TripletData` is optional, but is only None at the very end
/// of the associated iterator.
///
/// # Fields
///
/// * `triplet_data`: The `TripletData` associated with these coordinates. This is `None` if there
///   is no associated data.
/// * `x`: The x coordinate.
/// * `y`: The y coordinate.
struct CoordsData {
    triplet_data: Option<TripletData>,
    x: f64,
    y: f64,
}

/// An iterator-wrapping struct that yields `CoordsData` from another iterator.
///
/// # Layer 2
///
/// `CoordsIter` wraps around another iterator that yields `TripletData`, and yields `CoordsData`
/// calculated from the `TripletData` and the previous coordinates and deltas. It also keeps track
/// of whether it has yielded the tail coordinates yet.
///
/// # Type Parameters
///
/// * `I`: The type of the inner iterator. Must be an iterator over `TripletData`.
///
/// # Fields
///
/// * `inner`: The inner iterator that yields `TripletData`.
/// * `tail`: A boolean that indicates whether the end of the iterator has been reached,
///   at which point one more `CoordsData` is yielded with no associated `TripletData`.
/// * `prev_x_coord`: The x coordinate from the previous `CoordsData`.
/// * `prev_y_coord`: The y coordinate from the previous `CoordsData`.
/// * `prev_dx`: The delta x from the previous `TripletData`.
/// * `prev_dy`: The delta y from the previous `TripletData`.
struct CoordsIter<I: Iterator> {
    inner: I,
    tail: bool,
    prev_x_coord: f64,
    prev_y_coord: f64,
    prev_dx: f64,
    prev_dy: f64,
}

/// Implementation of the `Iterator` trait for `TripletWindowsIter` struct.
///
/// This iterator yields `TripletData` items, which are calculated based on the next three bases
/// as a sliding window from the inner iterator, as well as the current roll type.
///
/// # Type Parameters
///
/// * `I`: The type of the inner iterator. Must be an iterator over `u8`.
///
/// # Returns
///
/// The `next` method returns `Some(TripletData)` if there are enough items left in the inner
/// iterator, or `None` if there are not.
impl<I> Iterator for TripletWindowsIter<I>
where
    I: Iterator<Item = u8>,
{
    type Item = TripletData;

    fn next(&mut self) -> Option<Self::Item> {
        // Fill the buffer with the next three items from the inner iterator.
        while self.base_buffer.len() < matrix::TRIPLET_SIZE {
            if let Some(item) = self.inner.next() {
                self.base_buffer.push_back(item);
            } else {
                break;
            }
        }
        // When the buffer is full, calculate the twist, roll, and tilt values.
        if self.base_buffer.len() >= matrix::TRIPLET_SIZE {
            let triplet: Vec<u8> = self.base_buffer.iter().cloned().take(3).collect();
            let twist = matrix::matrix_lookup(&triplet, &matrix::TWIST).unwrap();
            let roll = match self.roll_type {
                matrix::RollType::Simple => {
                    matrix::matrix_lookup(&triplet, &matrix::ROLL_ACTIVE).unwrap()
                }
                matrix::RollType::Active => {
                    matrix::matrix_lookup(&triplet, &matrix::ROLL_SIMPLE).unwrap()
                }
            };
            let tilt = matrix::matrix_lookup(&triplet, &matrix::TILT).unwrap();
            self.twist_sum += twist;
            // Create a TripletData instance and return it.
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

    /// Implementation of `Iterator` trait for `CoordsIter` struct.
    ///
    /// # Layer 2
    ///
    /// This method first tries to get the next `TripletData` from the inner iterator. If there is a next item,
    /// it updates the previous deltas with the deltas from the current item and creates a new `CoordsData` with
    /// the current `TripletData`.
    ///
    /// If there are no more items in the inner iterator it yields one more new `CoordsData` without a
    /// `TripletData` but with `x` and `y` filled in.
    ///
    /// # Returns
    ///
    /// A `Some(CoordsData)` with the next coordinates and `TripletData`, or `None` if there are no more items.
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
    /// Creates a `CoordsData` instance from an optional `TripletData`.
    ///
    /// # Layer 2
    ///
    /// Helper to `CoordsIter::next()` that creates a `CoordsData` instance from the current
    /// `TripletData` and the previous coordinates.
    ///
    /// # Arguments
    ///
    /// * `triplet_data` - An optional `TripletData` that will be included in the created `CoordsData`.
    ///
    /// # Returns
    ///
    /// A `CoordsData` instance with the calculated coordinates and the given `TripletData`.
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

/// A trait for `u8` Iterators to yield `TripletData`.
///
/// # Layer 1
///
/// `TripletWindowsIterator` is a trait for iterators over `u8` that provides a method for
/// transforming the iterator into a `TripletWindowsIter`. This allows for convenient conversion
/// of any iterator over `u8` into an iterator that yields triplets of nucleotides.
///
/// # Type Parameters
///
/// * `Self`: The type implementing this trait. Must be an iterator over `u8`.
///
/// # Methods
///
/// * `triplet_windows_iter`: Takes a `RollType` and returns a `TripletWindowsIter` that yields
///   triplets of nucleotides from the original iterator.
trait TripletWindowsIterator: Iterator<Item = u8> + Sized {
    fn triplet_windows_iter(self, roll_type: matrix::RollType) -> TripletWindowsIter<Self> {
        TripletWindowsIter {
            base_buffer: VecDeque::new(),
            inner: self,
            twist_sum: 0.0,
            roll_type,
        }
    }
}

/// A trait for `TripletData` Iterators to yield `CoordsData`.
///
/// # Layer 2
///
/// `CoordsIterator` is a trait for iterators over `TripletData` that provides a method for
/// transforming the iterator into a `CoordsIter`. This allows for convenient conversion
/// of any iterator over `TripletData` into an iterator that yields `CoordsData`.
///
/// # Type Parameters
///
/// * `Self`: The type implementing this trait. Must be an iterator over `TripletData`.
///
/// # Methods
///
/// * `coords_iter`: Returns a `CoordsIter` that yields `CoordsData` calculated from the
///   `TripletData` yielded by the original iterator.
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

mod tests {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    fn test_triplet_iter() {
        let dna = b"ACGTACGT";
        let windows: Vec<TripletData> = dna
            .iter()
            .cloned()
            .triplet_windows_iter(matrix::RollType::Simple)
            .collect();
        assert_eq!(windows.len(), 6);
        let windows: Vec<TripletData> = dna
            .iter()
            .cloned()
            .triplet_windows_iter(matrix::RollType::Active)
            .collect();
        assert_eq!(windows.len(), 6);
    }

    #[test]
    fn test_triplet_iter_too_short() {
        let dna = b"AC";
        let windows: Vec<TripletData> = dna
            .iter()
            .cloned()
            .triplet_windows_iter(matrix::RollType::Simple)
            .collect();
        assert_eq!(windows.len(), 0);
    }

    #[test]
    fn test_coords_iter() {
        let dna = b"ACGTACGT";
        let windows: Vec<CoordsData> = dna
            .iter()
            .cloned()
            .triplet_windows_iter(matrix::RollType::Simple)
            .coords_iter()
            .collect();
        assert_eq!(windows.len(), 7);
    }
}
