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
#[derive(Clone, Debug)]
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
                    matrix::matrix_lookup(&triplet, &matrix::ROLL_SIMPLE).unwrap()
                }
                matrix::RollType::Active => {
                    matrix::matrix_lookup(&triplet, &matrix::ROLL_ACTIVE).unwrap()
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

/// A trait for `u8` Iterators to yield `TripletData`.
///
/// `TripletWindowsIterator` is a trait for iterators over `u8` that provides a method for
/// transforming the iterator into a `TripletWindowsIter`. This allows for convenient conversion
/// of any iterator over `u8` into an iterator that yields triplets of nucleotides. This is
/// **layer 1** of the iterator stack.
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

impl<I: Iterator<Item = u8>> TripletWindowsIterator for I {}

/// Represents the coordinates and associated data for a triplet of nucleotides.
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

impl CoordsData {
    /// Constructor for `CoordsData`.
    fn new(triplet_data: Option<TripletData>, x: f64, y: f64) -> Self {
        CoordsData { triplet_data, x, y }
    }
}

/// An iterator-wrapping struct that yields `CoordsData` from another iterator.
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
/// * `head`: A boolean that indicates whether the first `CoordsData` has been yielded yet.
/// * `tail`: A boolean that indicates whether the end of the iterator has been reached,
///   at which point one more `CoordsData` is yielded with no associated `TripletData`.
/// * `prev_x_coord`: The x coordinate from the previous `CoordsData`.
/// * `prev_y_coord`: The y coordinate from the previous `CoordsData`.
/// * `prev_dx`: The delta x from the previous `TripletData`.
/// * `prev_dy`: The delta y from the previous `TripletData`.
struct CoordsIter<I: Iterator> {
    inner: I,
    head: bool,
    tail: bool,
    prev_x_coord: f64,
    prev_y_coord: f64,
    prev_dx: f64,
    prev_dy: f64,
}

impl<I: Iterator<Item = TripletData>> CoordsIter<I> {
    /// Constructor for `CoordsIter`.
    fn new(inner: I) -> Self {
        CoordsIter {
            inner,
            head: false,
            tail: false,
            prev_x_coord: 0.0,
            prev_y_coord: 0.0,
            prev_dx: 0.0,
            prev_dy: 0.0,
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
            let result = Some(self.create_coords_data(Some(triplet_data.to_owned())));
            self.prev_dx = triplet_data.dx;
            self.prev_dy = triplet_data.dy;
            if !self.head {
                self.head = true;
                return self.next();
            }
            result
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

/// A trait for `TripletData` Iterators to yield `CoordsData`.
///
/// `CoordsIterator` is a trait for iterators over `TripletData` that provides a method for
/// transforming the iterator into a `CoordsIter`. This allows for convenient conversion
/// of any iterator over `TripletData` into an iterator that yields `CoordsData`. This
/// is **layer 2** of the iterator stack.
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
            head: false,
            tail: false,
            prev_x_coord: 0.0,
            prev_y_coord: 0.0,
            prev_dx: 0.0,
            prev_dy: 0.0,
        }
    }
}

impl<I: Iterator<Item = TripletData>> CoordsIterator for I {}

/// Represents the data for a rolling mean of the x and y coordinates.
///
/// # Fields
///
/// * `x_bar`: The weighted mean of the x coordinates.
/// * `y_bar`: The weighted mean of the y coordinates.
struct RollMeanData {
    x_bar: f64,
    y_bar: f64,
}

/// Represents the data for a rolling mean of the x and y coordinates.
///
/// The `RollMeanData` struct contains the weighted x and y means for a window of coordinates
/// that is 2 * `step_size` + 1 in length.
///
/// # Fields
///
/// * `inner`: The inner iterator that yields `CoordsData`.
/// * `buffer`: A buffer that stores the current window of coordinates.
/// * `step_size`: Half the size of the window minus one.  In other words,
///   2 * `step_size` + 1 is the size of the window.
/// * `x_roll_sum`: The sum of the x coordinates in the current window.
/// * `y_roll_sum`: The sum of the y coordinates in the current window.
struct RollMeanIter<I: Iterator> {
    inner: I,
    buffer: VecDeque<CoordsData>,
    step_size: usize,
    x_roll_sum: f64,
    y_roll_sum: f64,
}

/// Implementation of the `Iterator` trait for `RollMeanIter`.
///
/// This iterator wraps another iterator of items of type `CoordsData` and computes
/// a rolling mean of the `x` and `y` values of the items.
impl<I> Iterator for RollMeanIter<I>
where
    I: Iterator<Item = CoordsData>,
{
    type Item = RollMeanData;

    /// Computes the next item of the rolling mean iterator.
    ///
    /// This method computes the rolling mean of the `x` and `y` values of the next
    /// `window_size` items from the inner iterator, where `window_size` is `step_size * 2 + 1`.
    ///
    /// The method returns `Some(RollMeanData)` if there are enough items in the inner iterator,
    /// and `None` otherwise.
    fn next(&mut self) -> Option<Self::Item> {
        // Fill the buffer with the next three items from the inner iterator.
        let window_size = self.step_size * 2 + 1;
        while self.buffer.len() < window_size {
            if let Some(item) = self.inner.next() {
                self.x_roll_sum += item.x;
                self.y_roll_sum += item.y;
                self.buffer.push_back(item);
            } else {
                break;
            }
        }
        if self.buffer.len() >= window_size {
            // get the fron/back items without removing them and adjust the roll sum
            let adj_x_roll_sum = self.x_roll_sum
                - (0.5 * self.buffer.front().unwrap().x)
                - (0.5 * self.buffer.back().unwrap().x);
            let adj_y_roll_sum = self.y_roll_sum
                - (0.5 * self.buffer.front().unwrap().y)
                - (0.5 * self.buffer.back().unwrap().y);
            let x_bar = adj_x_roll_sum / (window_size as f64 - 1 as f64);
            let y_bar = adj_y_roll_sum / (window_size as f64 - 1 as f64);
            let result = Some(RollMeanData { x_bar, y_bar });
            let item = self.buffer.pop_front().unwrap();
            self.x_roll_sum -= item.x;
            self.y_roll_sum -= item.y;
            result
        } else {
            None
        }
    }
}

/// A trait for iterators that can compute a rolling mean of `CoordsData`.
///
/// This trait extends the `Iterator` trait, adding a `roll_mean_iter` method that
/// wraps the iterator in a `RollMeanIter`. The `RollMeanIter` computes a rolling mean
/// of the `x` and `y` values of the items from the original iterator.
trait RollMeanIterator: Iterator<Item = CoordsData> + Sized {
    /// Wraps the iterator in a `RollMeanIter`.
    ///
    /// This method takes ownership of the iterator and returns a `RollMeanIter` that
    /// computes a rolling mean of the `x` and `y` values of the items from the original iterator.
    ///
    /// # Parameters
    ///
    /// * `step_size`: half of the window size minus one. In other words, 2 * `step_size` + 1 is
    ///  the size of the window.
    ///
    /// # Returns
    ///
    /// A `RollMeanIter` that computes a rolling mean of the `x` and `y` values of the items.
    fn roll_mean_iter(self, step_size: usize) -> RollMeanIter<Self> {
        RollMeanIter {
            inner: self,
            buffer: VecDeque::new(),
            step_size,
            x_roll_sum: 0.0,
            y_roll_sum: 0.0,
        }
    }
}

impl<I: Iterator<Item = CoordsData>> RollMeanIterator for I {}

/// An iterator that computes the Euclidean distance between pairs of items from an inner iterator.
///
/// `EucDistIter` wraps another iterator that yields `RollMeanData`. It computes the Euclidean
/// distance between each pair of items from the inner iterator.
///
/// # Fields
///
/// * `inner`: The inner iterator that yields `RollMeanData`.
///
/// * `buffer`: A buffer that stores 2 * `curve_step_size` + 1 items from the inner iterator.
///
/// * `curve_step_size`: The distance from the midpoint base in the window.  
struct EucDistIter<I: Iterator> {
    inner: I,
    buffer: VecDeque<RollMeanData>,
    curve_step_size: usize,
}

impl<I> Iterator for EucDistIter<I>
where
    I: Iterator<Item = RollMeanData>,
{
    type Item = f64;

    /// Computes the next item of the Euclidean distance iterator.
    ///
    /// This method computes the Euclidean distance between each pair of consecutive items
    /// from the inner iterator. The Euclidean distance is computed as the square root of
    /// the sum of the squares of the differences of the `x_bar` and `y_bar` values of the items.
    ///
    /// The method returns `Some(f64)` if there are enough items in the inner iterator,
    /// and `None` otherwise.
    fn next(&mut self) -> Option<Self::Item> {
        // Fill the buffer with the next three items from the inner iterator.
        let window_size = self.curve_step_size * 2 + 1;
        while self.buffer.len() < window_size {
            if let Some(item) = self.inner.next() {
                self.buffer.push_back(item);
            } else {
                break;
            }
        }
        if self.buffer.len() >= window_size {
            let left = self.buffer.front().unwrap();
            let right = self.buffer.back().unwrap();
            let curve = ((right.y_bar - left.y_bar).powf(2.0)
                + (right.x_bar - left.x_bar).powf(2.0))
            .sqrt();
            self.buffer.pop_front();
            Some(curve)
        } else {
            None
        }
    }
}

trait EucDistIterator: Iterator<Item = RollMeanData> + Sized {
    fn euc_dist_iter(self, curve_step_size: usize) -> EucDistIter<Self> {
        EucDistIter {
            inner: self,
            buffer: VecDeque::new(),
            curve_step_size,
        }
    }
}

impl<I: Iterator<Item = RollMeanData>> EucDistIterator for I {}

/// An iterator that computes the curvature of a DNA sequence.
///
/// `CurveIter` wraps an iterator that yields `u8` and computes the curvature of the DNA sequence
/// represented by the nucleotides.
///
/// # Fields
///
/// * `inner`: The inner iterator that yields `u8`.
pub struct CurveIter<I: Iterator<Item = u8>> {
    inner: EucDistIter<RollMeanIter<CoordsIter<TripletWindowsIter<I>>>>,
}

impl<I: Iterator<Item = u8>> Iterator for CurveIter<I> {
    type Item = f64;

    /// Computes the next item of the curvature iterator.
    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

/// Construct a `CurveIter` from an iterator that yields `u8`.
///
/// This function constructs a `CurveIter` from an iterator that yields `u8`. The `CurveIter`
/// computes the curvature of the DNA sequence represented by the nucleotides.
///
/// # Parameters
///
/// * `seq_iter`: An iterator that yields `u8`.
/// * `roll_type`: The type of roll (either simple or activated).
/// * `step_b`: Half of the window size minus one. In other words, 2 * `step_size` + 1 is
///  the size of the window.
/// * `step_c`: The distance from the midpoint base to the sides in the curve window.
impl<I: Iterator<Item = u8>> CurveIter<I> {
    fn new(seq_iter: I, roll_type: matrix::RollType, step_b: usize, step_c: usize) -> Self {
        Self {
            inner: seq_iter
                .triplet_windows_iter(roll_type)
                .coords_iter()
                .roll_mean_iter(step_b)
                .euc_dist_iter(step_c),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Below is a table of some of the expected values for the coords iterator over the DNA
    ///
    /// | pos|nuc|trip | ixs |  twist |  roll_s |   tilt |twist_sum| dx_simp | dy_simp |
    /// | --:| -:| --: | --: | -----: | ------: | -----: | ------: | ------: | ------: |
    /// |  0 | C | CCA | 330 | 0.5986 |  0.7000 | 0.0000 |  0.5986 |  0.3945 |  0.5783 |
    /// |  1 | C | CAA | 300 | 0.5986 |  6.2000 | 0.0000 |  1.1973 |  5.7725 |  2.2622 |
    /// |  2 | A | AAC | 003 | 0.5986 |  1.6000 | 0.0000 |  1.7959 |  1.5596 | -0.3572 |
    /// |  3 | A | ACA | 030 | 0.5986 |  5.8000 | 0.0000 |  2.3946 |  3.9408 | -4.2556 |
    /// |  4 | C | CAT | 301 | 0.5986 |  8.7000 | 0.0000 |  2.9932 |  1.2860 | -8.6044 |
    /// |  5 | A | ATT | 011 | 0.5986 |  0.0000 | 0.0000 |  3.5919 |  0.0000 |  0.0000 |
    /// |  6 | T | TTT | 111 | 0.5986 |  0.1000 | 0.0000 |  4.1905 | -0.0867 | -0.0498 |
    /// |  7 | T | TTT | 111 | 0.5986 |  0.1000 | 0.0000 |  4.7892 | -0.0997 |  0.0077 |
    /// |  8 | T | TTG | 112 | 0.5986 |  6.2000 | 0.0000 |  5.3878 | -4.8387 |  3.8765 |
    /// |  9 | T | TGA | 120 | 0.5986 | 10.0000 | 0.0000 |  5.9865 | -2.9238 |  9.5630 |
    /// | 10 | G | GAC | 203 | 0.5986 |  5.6000 | 0.0000 |  6.5851 |  1.6653 |  5.3467 |
    /// | 11 | A | ACT | 031 | 0.5986 |  2.0000 | 0.0000 |  7.1838 |  1.5674 |  1.2423 |
    /// | 12 | C | CTT | 311 | 0.5986 |  4.2000 | 0.0000 |  7.7824 |  4.1892 |  0.3003 |
    /// | 13 | T | TTT | 111 | 0.5986 |  0.1000 | 0.0000 |  8.3811 |  0.0864 | -0.0503 |
    /// | 14 | T | TTT | 111 | 0.5986 |  0.1000 | 0.0000 |  8.9797 |  0.0431 | -0.0903 |
    /// | 15 | T | TTT | 111 | 0.5986 |  0.1000 | 0.0000 |  9.5784 | -0.0153 | -0.0988 |
    /// | 16 | T | TTG | 112 | 0.5986 |  6.2000 | 0.0000 | 10.1770 | -4.2363 | -4.5270 |
    /// | 17 | T | TGG | 122 | 0.5986 |  0.7000 | 0.0000 | 10.7757 | -0.6831 | -0.1527 |
    /// | 18 | G | GGG | 222 | 0.5986 |  5.7000 | 0.0000 | 11.3743 | -5.2961 |  2.1075 |
    /// | 19 | G | GGA | 220 | 0.5986 |  6.2000 | 0.0000 | 11.9729 | -3.4670 |  5.1400 |
    /// | 20 | G | GAG | 202 | 0.5986 |  6.6000 | 0.0000 | 12.5716 |  0.0345 |  6.5999 |
    /// | 21 | A | AGG | 022 | 0.5986 |  4.7000 | 0.0000 | 13.1702 |  2.6688 |  3.8688 |
    /// | 22 | G | GGG | 222 | 0.5986 |  5.7000 | 0.0000 | 13.7689 |  5.3178 |  2.0520 |
    /// | 23 | G | GGC | 223 | 0.5986 |  8.2000 | 0.0000 | 14.3675 |  7.9834 | -1.8724 |
    /// | 24 | G | GCA | 230 | 0.5986 |  7.5000 | 0.0000 | 14.9662 |  5.0670 | -5.5295 |
    /// | 25 | C | CAC | 303 | 0.5986 |  6.8000 | 0.0000 | 15.5648 |  0.9700 | -6.7305 |
    /// | 26 | A | ACT | 031 | 0.5986 |  2.0000 | 0.0000 | 16.1635 | -0.8799 | -1.7961 |
    /// | 27 | C | CTA | 310 | 0.5986 |  7.8000 | 0.0000 | 16.7621 | -6.7820 | -3.8528 |
    /// | 28 | T | TAG | 102 | 0.5986 |  7.8000 | 0.0000 | 17.3608 | -7.7738 |  0.6390 |
    /// | 29 | A | AGC | 023 | 0.5986 |  6.3000 | 0.0000 | 17.9594 | -4.8961 |  3.9646 |
    /// | 30 | G | GCA | 230 | 0.5986 |  7.5000 | 0.0000 | 18.5581 | -2.1553 |  7.1836 |
    /// | 31 | C | CAC | 303 | 0.5986 |  6.8000 | 0.0000 | 19.1567 |  2.0560 |  6.4817 |
    /// | 32 | A | ACC | 033 | 0.5986 |  5.2000 | 0.0000 | 19.7554 |  4.0920 |  3.2087 |
    /// | 33 | C | CCT | 331 | 0.5986 |  4.7000 | 0.0000 | 20.3540 |  4.6897 |  0.3116 |
    /// | 34 | C | CTA | 310 | 0.5986 |  7.8000 | 0.0000 | 20.9527 |  6.7208 | -3.9587 |
    /// | 35 | T | TAT | 101 | 0.5986 |  9.7000 | 0.0000 | 21.5513 |  4.1302 | -8.7767 |
    /// | 36 | A | ATC | 013 | 0.5986 |  3.6000 | 0.0000 | 22.1500 | -0.5693 | -3.5547 |
    /// | 37 | T | TCT | 131 | 0.5986 |  6.5000 | 0.0000 | 22.7486 | -4.4660 | -4.7228 |
    /// | 38 | C | CTA | 310 | 0.5986 |  7.8000 | 0.0000 | 23.3472 | -7.6209 | -1.6618 |
    /// | 39 | T | TAC | 103 | 0.5986 |  6.4000 | 0.0000 | 23.9459 | -5.9340 |  2.3974 |
    /// | 40 | A | ACC | 033 | 0.5986 |  5.2000 | 0.0000 | 24.5445 | -2.8853 |  4.3261 |
    /// | 41 | C | CCC | 333 | 0.5986 |  5.7000 | 0.0000 | 25.1432 |  0.0596 |  5.6997 |
    /// | 42 | C | CCT | 331 | 0.5986 |  4.7000 | 0.0000 | 25.7418 |  2.6890 |  3.8548 |
    /// | 43 | C | CTG | 312 | 0.5986 |  9.6000 | 0.0000 | 26.3405 |  8.9743 |  3.4092 |
    /// | 44 | T | TGA | 120 | 0.5986 | 10.0000 | 0.0000 | 26.9391 |  9.7238 | -2.3342 |
    /// | 45 | G | GAA | 200 | 0.5986 |  5.1000 | 0.0000 | 27.5378 |  3.4259 | -3.7780 |
    /// | 46 | A | AAT | 001 | 0.5986 |  0.0000 | 0.0000 | 28.1364 |  0.0000 |  0.0000 |
    /// | 47 | A | ATC | 013 | 0.5986 |  3.6000 | 0.0000 | 28.7351 | -1.6006 | -3.2246 |
    /// | 48 | T |     |     |         |        |        |         |         |         |
    /// | 49 | C |     |     |         |        |        |         |         |         |
    #[test]
    fn test_triplet_iter_long() {
        let dna = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let windows: Vec<TripletData> = dna
            .iter()
            .cloned()
            .triplet_windows_iter(matrix::RollType::Simple)
            .collect();
        assert_eq!(windows.len(), dna.len() - 2);
        // check first two
        assert_relative_eq!(windows[0].dx, 0.3945, epsilon = 1e-4);
        assert_relative_eq!(windows[0].dy, 0.5783, epsilon = 1e-4);
        assert_relative_eq!(windows[1].dx, 5.7725, epsilon = 1e-4);
        assert_relative_eq!(windows[1].dy, 2.2622, epsilon = 1e-4);
        // check last two
        assert_relative_eq!(windows[46].dx, 0.0000, epsilon = 1e-4);
        assert_relative_eq!(windows[46].dy, 0.0000, epsilon = 1e-4);
        assert_relative_eq!(windows[47].dx, -1.6006, epsilon = 1e-4);
        assert_relative_eq!(windows[47].dy, -3.2246, epsilon = 1e-4);
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

    /// Below is a table of some of the expected values for the coords iterator over the DNA
    ///
    /// | pos|nuc|trip | dx_simp | dy_simp |  x_coord |  y_coord |
    /// | --:| -:| --: | ------: | ------: | -------: | -------: |
    /// |  0 | C | CCA |  0.3945 |  0.5783 |          |          |
    /// |  1 | C | CAA |  5.7725 |  2.2622 |   0.3945 |   0.5783 |
    /// |  2 | A | AAC |  1.5596 | -0.3572 |   6.1670 |   2.8405 |
    /// |  3 | A | ACA |  3.9408 | -4.2556 |   7.7266 |   2.4833 |
    /// |  4 | C | CAT |  1.2860 | -8.6044 |  11.6674 |  -1.7723 |
    /// |  5 | A | ATT |  0.0000 |  0.0000 |  12.9534 | -10.3767 |
    /// |  6 | T | TTT | -0.0867 | -0.0498 |  12.9534 | -10.3767 |
    /// |  7 | T | TTT | -0.0997 |  0.0077 |  12.8667 | -10.4266 |
    /// |  8 | T | TTG | -4.8387 |  3.8765 |  12.7670 | -10.4189 |
    /// |  9 | T | TGA | -2.9238 |  9.5630 |   7.9283 |  -6.5424 |
    /// | 10 | G | GAC |  1.6653 |  5.3467 |   5.0045 |   3.0206 |
    /// | 11 | A | ACT |  1.5674 |  1.2423 |   6.6698 |   8.3673 |
    /// | 12 | C | CTT |  4.1892 |  0.3003 |   8.2372 |   9.6096 |
    /// | 13 | T | TTT |  0.0864 | -0.0503 |  12.4264 |   9.9099 |
    /// | 14 | T | TTT |  0.0431 | -0.0903 |  12.5128 |   9.8596 |
    /// | 15 | T | TTT | -0.0153 | -0.0988 |  12.5559 |   9.7693 |
    /// | 16 | T | TTG | -4.2363 | -4.5270 |  12.5406 |   9.6705 |
    /// | 17 | T | TGG | -0.6831 | -0.1527 |   8.3043 |   5.1435 |
    /// | 18 | G | GGG | -5.2961 |  2.1075 |   7.6212 |   4.9908 |
    /// | 19 | G | GGA | -3.4670 |  5.1400 |   2.3251 |   7.0983 |
    /// | 20 | G | GAG |  0.0345 |  6.5999 |  -1.1419 |  12.2383 |
    /// | 21 | A | AGG |  2.6688 |  3.8688 |  -1.1074 |  18.8382 |
    /// | 22 | G | GGG |  5.3178 |  2.0520 |   1.5614 |  22.7069 |
    /// | 23 | G | GGC |  7.9834 | -1.8724 |   6.8792 |  24.7590 |
    /// | 24 | G | GCA |  5.0670 | -5.5295 |  14.8626 |  22.8866 |
    /// | 25 | C | CAC |  0.9700 | -6.7305 |  19.9296 |  17.3571 |
    /// | 26 | A | ACT | -0.8799 | -1.7961 |  20.8995 |  10.6266 |
    /// | 27 | C | CTA | -6.7820 | -3.8528 |  20.0197 |   8.8305 |
    /// | 28 | T | TAG | -7.7738 |  0.6390 |  13.2377 |   4.9777 |
    /// | 29 | A | AGC | -4.8961 |  3.9646 |   5.4639 |   5.6167 |
    /// | 30 | G | GCA | -2.1553 |  7.1836 |   0.5678 |   9.5814 |
    /// | 31 | C | CAC |  2.0560 |  6.4817 |  -1.5875 |  16.7650 |
    /// | 32 | A | ACC |  4.0920 |  3.2087 |   0.4685 |  23.2467 |
    /// | 33 | C | CCT |  4.6897 |  0.3116 |   4.5605 |  26.4554 |
    /// | 34 | C | CTA |  6.7208 | -3.9587 |   9.2502 |  26.7669 |
    /// | 35 | T | TAT |  4.1302 | -8.7767 |  15.9709 |  22.8083 |
    /// | 36 | A | ATC | -0.5693 | -3.5547 |  20.1012 |  14.0315 |
    /// | 37 | T | TCT | -4.4660 | -4.7228 |  19.5319 |  10.4768 |
    /// | 38 | C | CTA | -7.6209 | -1.6618 |  15.0659 |   5.7540 |
    /// | 39 | T | TAC | -5.9340 |  2.3974 |   7.4450 |   4.0922 |
    /// | 40 | A | ACC | -2.8853 |  4.3261 |   1.5109 |   6.4896 |
    /// | 41 | C | CCC |  0.0596 |  5.6997 |  -1.3743 |  10.8157 |
    /// | 42 | C | CCT |  2.6890 |  3.8548 |  -1.3148 |  16.5154 |
    /// | 43 | C | CTG |  8.9743 |  3.4092 |   1.3742 |  20.3701 |
    /// | 44 | T | TGA |  9.7238 | -2.3342 |  10.3485 |  23.7794 |
    /// | 45 | G | GAA |  3.4259 | -3.7780 |  20.0722 |  21.4451 |
    /// | 46 | A | AAT |  0.0000 |  0.0000 |  23.4981 |  17.6671 |
    /// | 47 | A | ATC | -1.6006 | -3.2246 |  23.4981 |  17.6671 |
    /// | 48 | T |     |         |         |  21.8975 |  14.4425 |
    /// | 49 | C |     |         |         |          |          |
    #[test]
    fn test_coords_iter() {
        let dna = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let windows: Vec<CoordsData> = dna
            .iter()
            .cloned()
            .triplet_windows_iter(matrix::RollType::Simple)
            .coords_iter()
            .collect();
        assert_eq!(windows.len(), dna.len() - 2);
        assert_relative_eq!(windows[0].x, 0.3945, epsilon = 1e-4);
        assert_relative_eq!(windows[0].y, 0.5783, epsilon = 1e-4);
        assert_relative_eq!(windows[1].x, 6.1670, epsilon = 1e-4);
        assert_relative_eq!(windows[1].y, 2.8405, epsilon = 1e-4);
        assert_relative_eq!(windows[46].x, 23.4981, epsilon = 1e-4);
        assert_relative_eq!(windows[46].y, 17.6671, epsilon = 1e-4);
        assert_relative_eq!(windows[47].x, 21.8975, epsilon = 1e-4);
        assert_relative_eq!(windows[47].y, 14.4425, epsilon = 1e-4);
    }

    /// Helper for test_rollmean_iter()
    fn get_some_coords() -> Vec<CoordsData> {
        let x_values = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        let y_values = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
        ];

        x_values
            .into_iter()
            .zip(y_values.into_iter())
            .map(|(x, y)| CoordsData::new(None, x, y))
            .collect()
    }

    #[test]
    fn test_rollmean_iter() {
        let rolls: Vec<_> = get_some_coords().into_iter().roll_mean_iter(2).collect();
        assert_eq!(rolls.len(), 8);
        // x̄₃ = (½x₁ + x₂ + x₃ + x₄ + ½x₅)/4
        // x̄₃ = (0.5 + 2 + 3 + 4 + 2.5)/4 = 3
        assert_relative_eq!(rolls[0].x_bar, 3.0, epsilon = 1e-4);
        assert_relative_eq!(rolls[0].y_bar, 0.0, epsilon = 1e-4);
        // x̄₃ = (½x₂ + x₃ + x₄ + x₅ + ½x₆)/4
        // x̄₃ = (1 + 3 + 4 + 5 + 3)/4 = 16 / 4 = 4
        assert_relative_eq!(rolls[1].x_bar, 4.0, epsilon = 1e-4);
        assert_relative_eq!(rolls[2].x_bar, 5.0, epsilon = 1e-4);
        assert_relative_eq!(rolls[7].y_bar, 10.0, epsilon = 1e-4);
        let rolls: Vec<_> = get_some_coords().into_iter().roll_mean_iter(3).collect();
        // x̄₃ = (½x₁ + x₂ + x₃ + x₄ + x₅ + x₆+ ½x₇)/6
        // x̄₃ = (0.5 + 2 + 3 + 4 + 5 + 6 + 3.5)/6 = 24 / 6 = 4
        assert_relative_eq!(rolls[0].x_bar, 4.0, epsilon = 1e-4);
        assert_eq!(rolls.len(), 6);
    }

    /// Helper for test_eucdist_iter()
    fn get_some_means() -> Vec<RollMeanData> {
        let x_values = vec![3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.0, 5.0, 17.0];
        let y_values = vec![0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0];

        x_values
            .into_iter()
            .zip(y_values.into_iter())
            .map(|(x_bar, y_bar)| RollMeanData { x_bar, y_bar })
            .collect()
    }

    #[test]
    fn test_eucdist_iter() {
        let mean_rolls: Vec<_> = get_some_means();
        let vec_size = mean_rolls.len();
        let curve_step_size = 2;
        let euc_dists: Vec<_> = mean_rolls
            .into_iter()
            .euc_dist_iter(curve_step_size)
            .collect();
        // check curve_step_size number of items on both flanks are discarded
        assert_eq!(euc_dists.len(), vec_size - 2 * (curve_step_size));
        // √((7.0-3.0)² + (10.0-0.0)²) = √116 = 10.770329614269007
        assert_relative_eq!(euc_dists[0], 10.7703, epsilon = 1e-4);
        // √((8.0-4.0)² + (10.0-0.0)²) = √116 = 10.770329614269007
        assert_relative_eq!(euc_dists[1], 10.7703, epsilon = 1e-4);
        // √((8.0-5.0)² + (10.0-0.0)²) = √109 = 10.44031
        assert_relative_eq!(euc_dists[2], 10.44031, epsilon = 1e-4);
        // √((5.0-6.0)² + (10.0-0.0)²) = √101 = 10.04988
        assert_relative_eq!(euc_dists[3], 10.04988, epsilon = 1e-4);
        // √((17.0-7.0)² + (10.0-10.0)²) = √100 = 10.0
        assert_relative_eq!(euc_dists[4], 10.0, epsilon = 1e-4);
    }

    #[test]
    fn test_curve_iter() {
        let seq = b"CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let seq_len = seq.len();
        let curves: Vec<_> =
            CurveIter::new(seq.iter().cloned(), matrix::RollType::Simple, 5, 15).collect();
        assert_eq!(curves.len(), seq_len - (21 * 2));
    }
}
