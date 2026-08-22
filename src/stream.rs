//! Streaming a FASTA file through the curvature calculation.
//!
//! Scores are handed to a callback as they are produced and then dropped, so the amount
//! held at once is a batch of chunks rather than a genome. What is not bounded here is
//! the sequence itself: without an index the FASTA reader yields a whole record at a
//! time, so the largest record is a floor on the footprint.

use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use noodles_core::{Position, Region};
use noodles_fasta::fai;

use rayon::prelude::*;

use crate::curve::scan::{CurveParams, chunks, score_bases};
use crate::fasta::{RecordPiece, split_seq_by_gaps};

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

/// Score a run of already gap-split pieces, emitting each score with its absolute
/// position.
///
/// `base_offset` is how far the scored slice sits into the record: zero when the whole
/// record was read, and the window's start when reading through an index. `owned`
/// optionally restricts what is emitted to a 1-based inclusive absolute range, which is
/// how overlapping windows avoid emitting the same position twice.
#[allow(clippy::too_many_arguments)]
fn emit_pieces<F>(
    name: &str,
    pieces: &[RecordPiece],
    base_offset: usize,
    params: &CurveParams,
    chunk_scores: usize,
    batch: usize,
    owned: Option<(usize, usize)>,
    emit: &mut F,
) -> io::Result<usize>
where
    F: FnMut(&str, usize, f64) -> io::Result<()>,
{
    let mut emitted = 0usize;
    for piece in pieces {
        let bases = piece.bases();
        let piece_start = usize::from(piece.start);
        let plan = chunks(bases.len(), params, chunk_scores);

        // Score a batch of chunks at a time: enough to keep every thread busy, few
        // enough that only that many chunks of scores exist at once.
        for group in plan.chunks(batch.max(1)) {
            let scored: Vec<Vec<f64>> = group
                .par_iter()
                .map(|chunk| score_bases(&bases[chunk.in_start..chunk.in_end], params))
                .collect();
            for (chunk, values) in group.iter().zip(&scored) {
                // out_start is an offset into the piece's scores; the first score sits
                // lead_in bases into the piece.
                let first = base_offset + piece_start + params.lead_in() + chunk.out_start;
                for (i, &value) in values.iter().enumerate() {
                    let position = first + i;
                    if let Some((lo, hi)) = owned
                        && (position < lo || position > hi)
                    {
                        continue;
                    }
                    emit(name, position, value)?;
                    emitted += 1;
                }
            }
        }
    }
    Ok(emitted)
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

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let length = record.sequence().len();
        stats.records += 1;
        stats.longest_record = stats.longest_record.max(length);

        let pieces = split_seq_by_gaps(record);
        stats.pieces += pieces.len();
        stats.scores += emit_pieces(
            &name,
            &pieces,
            0,
            params,
            chunk_scores,
            batch,
            None,
            &mut emit,
        )?;
    }
    Ok(stats)
}

/// Score every record, reading each one a window at a time through a FASTA index.
///
/// Without an index the reader hands over a whole record, so the largest record is a
/// floor on memory. Querying windows lowers that floor to the window size.
///
/// Windows overlap by `lead_in` bases and each emits only the range it owns. The overlap
/// is what makes a windowed run agree with a whole-record one: a piece straddling a
/// boundary is truncated within the window and so loses `lead_in` scores at that
/// artificial edge, but those positions belong to the neighbouring window, which sees
/// them with full context.
pub fn for_each_score_indexed<F>(
    input: &Path,
    index: &fai::Index,
    params: &CurveParams,
    chunk_scores: usize,
    batch: usize,
    window_bases: usize,
    mut emit: F,
) -> io::Result<RunStats>
where
    F: FnMut(&str, usize, f64) -> io::Result<()>,
{
    let inner = BufReader::new(File::open(input)?);
    let mut reader = noodles_fasta::io::IndexedReader::new(inner, index.clone());
    let lead = params.lead_in();
    let window = window_bases.max(2 * lead + 1);
    let mut stats = RunStats::default();

    for record in index.as_ref() {
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let length = record.length() as usize;
        stats.records += 1;
        stats.longest_record = stats.longest_record.max(length);
        if length == 0 {
            continue;
        }

        let mut own_start = 1usize; // 1-based, inclusive
        while own_start <= length {
            let own_end = (own_start + window - 1).min(length);
            // Read the owned range plus context on each side, clamped to the record.
            let read_start = own_start.saturating_sub(lead).max(1);
            let read_end = (own_end + lead).min(length);

            let region = Region::new(
                name.as_str(),
                Position::try_from(read_start).map_err(io::Error::other)?
                    ..=Position::try_from(read_end).map_err(io::Error::other)?,
            );
            let slice = reader.query(&region)?;
            let pieces = split_seq_by_gaps(slice);
            if own_start == 1 {
                stats.pieces += pieces.len();
            }
            stats.scores += emit_pieces(
                &name,
                &pieces,
                read_start - 1,
                params,
                chunk_scores,
                batch,
                Some((own_start, own_end)),
                &mut emit,
            )?;

            own_start = own_end + 1;
        }
    }
    Ok(stats)
}

/// Read `<input>.fai` if it is there.
///
/// A missing index is not an error: it only means falling back to whole-record reads.
/// A malformed one is an error, since silently ignoring it would quietly give up the
/// lower memory the user was expecting.
pub fn load_index(input: &Path) -> io::Result<Option<fai::Index>> {
    let mut path = input.as_os_str().to_owned();
    path.push(".fai");
    let path = PathBuf::from(path);
    if !path.exists() {
        return Ok(None);
    }
    fai::fs::read(&path)
        .map(Some)
        .map_err(|e| io::Error::other(format!("cannot read {}: {e}", path.display())))
}

/// Chromosome names and lengths straight from an index, with no sequence read at all.
pub fn chrom_sizes_from_index(index: &fai::Index) -> Vec<ChromSize> {
    index
        .as_ref()
        .iter()
        .map(|r| {
            (
                String::from_utf8_lossy(r.name()).into_owned(),
                r.length() as u32,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::matrix::RollType;
    use noodles_fasta::io::Indexer;

    fn params() -> CurveParams {
        CurveParams {
            roll_type: RollType::Simple,
            step_b: 5,
            step_c: 15,
            curve_scale: 0.33335,
        }
    }

    fn tmp(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("symcurve-stream-{}-{}", std::process::id(), name));
        p
    }

    /// Write a FASTA with several records, gaps and soft-masked runs, plus its index.
    fn write_fasta_and_index(tag: &str) -> PathBuf {
        let unit = "CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
        let mut text = String::new();
        // Deliberately uneven: gaps and masked runs at awkward offsets so pieces
        // straddle window boundaries rather than lining up with them.
        for (name, spec) in [
            (
                "chrA",
                vec![
                    (unit.repeat(31), "NNNNN"),
                    (unit.repeat(17), "RR"),
                    (unit.repeat(23), ""),
                ],
            ),
            (
                "chrB",
                vec![(unit.repeat(9).to_lowercase(), "NN"), (unit.repeat(41), "")],
            ),
            ("chrC", vec![(unit.repeat(5), "")]),
        ] {
            text.push_str(&format!(">{name}\n"));
            let mut seq = String::new();
            for (body, gap) in spec {
                seq.push_str(&body);
                seq.push_str(gap);
            }
            for line in seq.as_bytes().chunks(60) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
        }
        let path = tmp(&format!("{tag}.fa"));
        std::fs::write(&path, text).unwrap();

        // Build the .fai alongside it.
        let mut indexer = Indexer::new(BufReader::new(File::open(&path).unwrap()));
        let mut records = Vec::new();
        while let Some(record) = indexer.index_record().unwrap() {
            records.push(record);
        }
        let index = fai::Index::from(records);
        let mut fai_path = path.as_os_str().to_owned();
        fai_path.push(".fai");
        fai::fs::write(PathBuf::from(fai_path), &index).unwrap();
        path
    }

    fn collect(
        run: impl FnOnce(&mut dyn FnMut(&str, usize, f64) -> io::Result<()>) -> io::Result<RunStats>,
    ) -> (Vec<(String, usize, f64)>, RunStats) {
        let mut out = Vec::new();
        let stats = run(&mut |name, position, score| {
            out.push((name.to_string(), position, score));
            Ok(())
        })
        .unwrap();
        (out, stats)
    }

    #[test]
    fn test_indexed_and_whole_record_reads_agree() {
        // The windowed path is only worth having if it computes the same thing. Windows
        // are chosen small and awkward so that pieces and gaps straddle their edges.
        let path = write_fasta_and_index("agree");
        let index = load_index(&path).unwrap().expect("index should be found");
        let p = params();

        let (plain, plain_stats) = collect(|emit| for_each_score(&path, &p, 1000, 4, emit));
        assert!(plain.len() > 5000, "test input too small to be meaningful");

        for window in [2 * p.lead_in() + 1, 97, 1000, 7919, 1 << 20] {
            let (windowed, windowed_stats) =
                collect(|emit| for_each_score_indexed(&path, &index, &p, 1000, 4, window, emit));
            assert_eq!(
                windowed.len(),
                plain.len(),
                "score count differs at window {window}"
            );
            assert_eq!(windowed_stats.records, plain_stats.records);
            for (a, b) in plain.iter().zip(&windowed) {
                assert_eq!(a.0, b.0, "chromosome differs at window {window}");
                assert_eq!(a.1, b.1, "position differs at window {window}");
                // Not bit-identical: a window accumulates twist over a different run.
                let rel = (a.2 - b.2).abs() / a.2.abs().max(1e-12);
                assert!(rel < 1e-6, "value {} vs {} at window {window}", a.2, b.2);
            }
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_positions_are_emitted_in_order_and_without_duplicates() {
        let path = write_fasta_and_index("order");
        let index = load_index(&path).unwrap().unwrap();
        let (out, _) =
            collect(|emit| for_each_score_indexed(&path, &index, &params(), 500, 4, 331, emit));
        let mut by_chrom: Vec<(&str, usize)> =
            out.iter().map(|(c, p, _)| (c.as_str(), *p)).collect();
        let before = by_chrom.len();
        by_chrom.dedup();
        assert_eq!(by_chrom.len(), before, "a position was emitted twice");
        // Within each chromosome positions must strictly increase.
        for pair in out.windows(2) {
            if pair[0].0 == pair[1].0 {
                assert!(
                    pair[0].1 < pair[1].1,
                    "out of order: {:?} then {:?}",
                    pair[0],
                    pair[1]
                );
            }
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_chrom_sizes_from_index_match_reading_the_file() {
        let path = write_fasta_and_index("sizes");
        let index = load_index(&path).unwrap().unwrap();
        assert_eq!(chrom_sizes_from_index(&index), chrom_sizes(&path).unwrap());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_missing_index_is_not_an_error() {
        let path = tmp("no-index.fa");
        std::fs::write(&path, ">c\nACGT\n").unwrap();
        assert!(load_index(&path).unwrap().is_none());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_malformed_index_is_an_error() {
        // Silently ignoring a broken index would quietly give up the lower memory the
        // user asked for.
        let path = tmp("bad-index.fa");
        std::fs::write(&path, ">c\nACGT\n").unwrap();
        let mut fai = path.as_os_str().to_owned();
        fai.push(".fai");
        let fai = PathBuf::from(fai);
        std::fs::write(&fai, "this is not an index\n").unwrap();
        assert!(load_index(&path).is_err());
        std::fs::remove_file(&path).ok();
        std::fs::remove_file(&fai).ok();
    }
}
