//! # Curve module
//! This module contains functions for calculation of DNA curvature, and is divided into
//! several submodules.

// TODO: drop this allow once the unused TripletData/CoordsData fields are removed.
#[allow(dead_code)]
pub mod iters;
#[allow(dead_code)]
pub mod matrix;
pub mod scan;
