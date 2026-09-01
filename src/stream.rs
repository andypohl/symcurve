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

use crate::curve::calls::{CallParams, NucleosomeCall, call_nucleosomes, greedy_non_overlapping};
use crate::curve::scan::{
    CurveParams, Stage, SymParams, chunks, score_bases, score_symmetry_bases, sym_chunks,
};
use crate::fasta::{RecordPiece, split_seq_by_gaps};

/// Everything a scoring run needs beyond the input file itself.
///
/// Grouped rather than passed one at a time: these travel together through every entry
/// point here, and a function taking eight positional values is easy to call wrongly.
#[derive(Debug, Clone, Copy)]
pub struct ScanConfig {
    pub params: CurveParams,
    pub sym: SymParams,
    pub stage: Stage,
    /// Scores produced per chunk, which is how the memory budget is applied.
    pub chunk_scores: usize,
    /// Chunks scored at once, normally the worker count.
    pub batch: usize,
    /// Bases held at once when reading through an index.
    pub window_bases: usize,
}

impl ScanConfig {
    /// How many bases at each end of a piece receive no score, for the selected stage.
    pub fn lead_in(&self) -> usize {
        match self.stage {
            Stage::Curvature => self.params.lead_in(),
            Stage::Symmetry => self.sym.lead_in(&self.params),
        }
    }

    /// The gap between consecutive output positions.
    pub fn stride(&self) -> usize {
        match self.stage {
            Stage::Curvature => 1,
            Stage::Symmetry => self.sym.step.max(1),
        }
    }
}

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
    config: &ScanConfig,
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
        let (lead, stride) = (config.lead_in(), config.stride());
        let plan = match config.stage {
            Stage::Curvature => chunks(bases.len(), &config.params, config.chunk_scores),
            Stage::Symmetry => sym_chunks(
                bases.len(),
                &config.params,
                &config.sym,
                config.chunk_scores,
            ),
        };

        // Score a batch of chunks at a time: enough to keep every thread busy, few
        // enough that only that many chunks of scores exist at once.
        for group in plan.chunks(config.batch.max(1)) {
            let scored: Vec<Vec<f64>> = group
                .par_iter()
                .map(|chunk| {
                    let slice = &bases[chunk.in_start..chunk.in_end];
                    match config.stage {
                        Stage::Curvature => score_bases(slice, &config.params),
                        Stage::Symmetry => score_symmetry_bases(slice, &config.params, &config.sym),
                    }
                })
                .collect();
            for (chunk, values) in group.iter().zip(&scored) {
                // out_start is an offset into the piece's scores; the first score sits
                // lead bases into the piece, and successive scores are `stride` apart.
                let first = base_offset + piece_start + lead + chunk.out_start * stride;
                for (i, &value) in values.iter().enumerate() {
                    let position = first + i * stride;
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
pub fn for_each_score<F>(input: &Path, config: &ScanConfig, mut emit: F) -> io::Result<RunStats>
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
        stats.scores += emit_pieces(&name, &pieces, 0, config, None, &mut emit)?;
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
    config: &ScanConfig,
    mut emit: F,
) -> io::Result<RunStats>
where
    F: FnMut(&str, usize, f64) -> io::Result<()>,
{
    let inner = BufReader::new(File::open(input)?);
    let mut reader = noodles_fasta::io::IndexedReader::new(inner, index.clone());
    let lead = config.lead_in();
    let window = config.window_bases.max(2 * lead + 1);
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
                config,
                Some((own_start, own_end)),
                &mut emit,
            )?;

            own_start = own_end + 1;
        }
    }
    Ok(stats)
}

/// Score every record and hand its nucleosome calls to `on_record`, with the sequence.
///
/// Calls need two things the per-score path does not provide: the record's bases, for the
/// sequence the reference puts in the GFF attribute column, and every call for a record at
/// once, because the greedy selection ranks them against each other. Both are per record,
/// so this reads whole records rather than windows even when an index is available.
///
/// Only dyads that actually score are retained, which is a small fraction of positions, so
/// what is held is the calls for one record rather than its scores.
pub fn for_each_record_calls<F>(
    input: &Path,
    config: &ScanConfig,
    call_params: &CallParams,
    greedy: bool,
    mut on_record: F,
) -> io::Result<RunStats>
where
    F: FnMut(&str, &[u8], &[NucleosomeCall], usize) -> io::Result<()>,
{
    let mut reader = open(input)?;
    let mut stats = RunStats::default();

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let length = record.sequence().len();
        stats.records += 1;
        stats.longest_record = stats.longest_record.max(length);

        // Keep the sequence alive for the attribute column while the pieces borrow it.
        let sequence: Vec<u8> = record.sequence().as_ref().to_vec();
        let pieces = split_seq_by_gaps(record);
        stats.pieces += pieces.len();

        let mut scored: Vec<(usize, f64)> = Vec::new();
        emit_pieces(
            &name,
            &pieces,
            0,
            config,
            None,
            &mut |_, position, score| {
                if score > 0.0 {
                    // The reference indexes its arrays from zero; positions here are 1-based.
                    scored.push((position - 1, score));
                }
                Ok(())
            },
        )?;

        let calls = call_nucleosomes(&scored, length, call_params);
        let calls = if greedy {
            greedy_non_overlapping(&calls, call_params)
        } else {
            calls
        };
        stats.scores += calls.len();
        on_record(&name, &sequence, &calls, call_params.half_width)?;
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

    /// A small symmetry window, so test inputs need not be thousands of bases long.
    fn sym() -> SymParams {
        SymParams { win: 20, step: 1 }
    }

    fn config(stage: Stage, chunk_scores: usize, window_bases: usize) -> ScanConfig {
        ScanConfig {
            params: params(),
            sym: sym(),
            stage,
            chunk_scores,
            batch: 4,
            window_bases,
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

        let (plain, plain_stats) =
            collect(|emit| for_each_score(&path, &config(Stage::Curvature, 1000, 0), emit));
        assert!(plain.len() > 5000, "test input too small to be meaningful");

        for window in [2 * params().lead_in() + 1, 97, 1000, 7919, 1 << 20] {
            let (windowed, windowed_stats) = collect(|emit| {
                for_each_score_indexed(&path, &index, &config(Stage::Curvature, 1000, window), emit)
            });
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
        let (out, _) = collect(|emit| {
            for_each_score_indexed(&path, &index, &config(Stage::Curvature, 500, 331), emit)
        });
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

    #[test]
    fn test_zero_length_index_records_are_skipped() {
        let path = write_fasta_and_index("zerolen");
        let index = load_index(&path).unwrap().unwrap();
        let config = config(Stage::Curvature, 1000, 700);
        let count = |index: &fai::Index| {
            let n = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
            let seen = std::sync::Arc::clone(&n);
            let stats = for_each_score_indexed(&path, index, &config, move |_, _, _| {
                seen.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                Ok(())
            })
            .unwrap();
            (stats, n.load(std::sync::atomic::Ordering::Relaxed))
        };
        let (plain, plain_n) = count(&index);

        // An empty record must be counted but never queried, since there is no range.
        let mut records: Vec<fai::Record> = index.as_ref().to_vec();
        let one = std::num::NonZero::new(1).unwrap();
        records.push(fai::Record::new("empty", 0, 0, one, one));
        let (padded, n) = count(&fai::Index::from(records));
        assert_eq!(padded.records, plain.records + 1);
        assert_eq!(padded.scores, plain.scores);
        assert_eq!(n, plain_n);

        std::fs::remove_file(&path).ok();
        let mut fai_path = path.as_os_str().to_owned();
        fai_path.push(".fai");
        std::fs::remove_file(PathBuf::from(fai_path)).ok();
    }
}
