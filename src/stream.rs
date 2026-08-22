//! Streaming a FASTA file through the curvature calculation.
//!
//! Scores are handed to a callback as they are produced and then dropped, so the amount
//! held at once is a batch of chunks rather than a genome. What is not bounded here is
//! the sequence itself: without an index the FASTA reader yields a whole record at a
//! time, so the largest record is a floor on the footprint.

use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

use rayon::prelude::*;

use crate::curve::scan::{CurveParams, chunks, score_bases};
use crate::fasta::split_seq_by_gaps;

/// A chromosome name and its length, as needed for a bigWig header.
pub type ChromSize = (String, u32);

/// Read just the names and lengths of each record.
///
/// bigWig stores chromosome sizes in its header, so they must all be known before any
/// value is written. Sequences are read and dropped rather than retained, which costs a
/// second pass over the file but keeps only one record in memory at a time.
pub fn chrom_sizes(input: &Path) -> io::Result<Vec<ChromSize>> {
    let mut reader = open(input)?;
    let mut sizes = Vec::new();
    for result in reader.records() {
        let record = result?;
        sizes.push((
            String::from_utf8_lossy(record.name()).into_owned(),
            record.sequence().len() as u32,
        ));
    }
    Ok(sizes)
}

fn open(input: &Path) -> io::Result<noodles_fasta::io::Reader<BufReader<File>>> {
    Ok(noodles_fasta::io::Reader::new(BufReader::new(File::open(
        input,
    )?)))
}

/// How the scoring run went, for reporting.
#[derive(Debug, Default, Clone, Copy)]
pub struct RunStats {
    pub records: usize,
    pub pieces: usize,
    pub scores: usize,
    /// Bases in the largest record, which is the unavoidable part of the footprint.
    pub longest_record: usize,
}

/// Score every record in `input`, passing each score to `emit` as it is produced.
///
/// `emit` receives the record name, the 1-based position the score belongs to, and the
/// score. It is called in position order within a record and in file order across
/// records, which is the order a bigWig writer requires.
pub fn for_each_score<F>(
    input: &Path,
    params: &CurveParams,
    chunk_scores: usize,
    batch: usize,
    mut emit: F,
) -> io::Result<RunStats>
where
    F: FnMut(&str, usize, f64) -> io::Result<()>,
{
    let mut reader = open(input)?;
    let mut stats = RunStats::default();
    let batch = batch.max(1);

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let length = record.sequence().len();
        stats.records += 1;
        stats.longest_record = stats.longest_record.max(length);

        let pieces = split_seq_by_gaps(record);
        stats.pieces += pieces.len();

        for piece in &pieces {
            let bases = piece.bases();
            let piece_start = usize::from(piece.start);
            let plan = chunks(bases.len(), params, chunk_scores);

            // Score a batch of chunks at a time: enough to keep every thread busy, few
            // enough that only that many chunks of scores exist at once.
            for group in plan.chunks(batch) {
                let scored: Vec<Vec<f64>> = group
                    .par_iter()
                    .map(|chunk| score_bases(&bases[chunk.in_start..chunk.in_end], params))
                    .collect();
                for (chunk, values) in group.iter().zip(&scored) {
                    // out_start is an offset into the piece's scores; the first score
                    // sits lead_in bases into the piece.
                    let first = piece_start + params.lead_in() + chunk.out_start;
                    for (i, &value) in values.iter().enumerate() {
                        emit(&name, first + i, value)?;
                        stats.scores += 1;
                    }
                }
            }
        }
    }
    Ok(stats)
}
