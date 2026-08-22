//! Writing curvature scores out, in the format implied by the output path.

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use bigtools::BigWigWrite;
use bigtools::Value;
use bigtools::beddata::BedParserStreamingIterator;

use crate::curve::scan::PieceCurves;

/// The output formats this tool can write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// Indexed binary bigWig, via bigtools.
    BigWig,
    /// Plain-text bedGraph: `chrom<TAB>start<TAB>end<TAB>value`, 0-based half-open.
    BedGraph,
}

impl OutputFormat {
    /// Pick a format from the output path's extension.
    ///
    /// Chosen from the extension rather than a flag so that the output path alone says
    /// what the file is, and an unrecognized extension is refused rather than guessed at:
    /// silently writing text into a file named `.bw` would produce something no genome
    /// browser can read.
    pub fn from_path(path: &Path) -> Result<Self, UnknownFormat> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or_default()
            .to_ascii_lowercase();
        match ext.as_str() {
            "bw" | "bigwig" => Ok(OutputFormat::BigWig),
            "bedgraph" | "bg" => Ok(OutputFormat::BedGraph),
            _ => Err(UnknownFormat {
                path: path.display().to_string(),
            }),
        }
    }
}

/// The output path's extension did not name a format this tool can write.
#[derive(Debug)]
pub struct UnknownFormat {
    path: String,
}

impl fmt::Display for UnknownFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "cannot tell the output format from {:?}: expected one of .bw, .bigWig, .bedGraph, .bg",
            self.path
        )
    }
}

impl std::error::Error for UnknownFormat {}

/// One scored record: its name, its full length, and the scores for each of its pieces.
pub struct ScoredRecord {
    pub name: String,
    pub length: usize,
    pub pieces: Vec<PieceCurves>,
}

/// Flatten scored records into per-base intervals, in the order they must be written.
///
/// bigWig requires values sorted by chromosome and then by position; records arrive in
/// file order and pieces within a record in position order, so emitting them as they come
/// preserves that. Positions are converted from the 1-based inclusive convention used
/// throughout the curve code to the 0-based half-open convention both formats use.
fn intervals(records: &[ScoredRecord]) -> impl Iterator<Item = (String, Value)> + '_ {
    records.iter().flat_map(|record| {
        record.pieces.iter().flat_map(move |piece| {
            piece.curves.iter().enumerate().map(move |(i, &score)| {
                let start = (piece.curve_start - 1 + i) as u32;
                (
                    record.name.clone(),
                    Value {
                        start,
                        end: start + 1,
                        value: score as f32,
                    },
                )
            })
        })
    })
}

/// Write scores to `path`, choosing the format from its extension.
pub fn write(
    path: &Path,
    records: &[ScoredRecord],
) -> Result<OutputFormat, Box<dyn std::error::Error>> {
    let format = OutputFormat::from_path(path)?;
    match format {
        OutputFormat::BedGraph => write_bedgraph(path, records)?,
        OutputFormat::BigWig => write_bigwig(path, records)?,
    }
    Ok(format)
}

/// Write a plain-text bedGraph.
pub fn write_bedgraph(path: &Path, records: &[ScoredRecord]) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    for (chrom, value) in intervals(records) {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            chrom, value.start, value.end, value.value
        )?;
    }
    out.flush()
}

/// Write an indexed bigWig.
pub fn write_bigwig(
    path: &Path,
    records: &[ScoredRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    // bigWig carries chromosome sizes in its header, so they must all be known before
    // any value is written. Full record lengths are used, not piece lengths: the
    // coordinates in the file refer to the record.
    let chrom_sizes: HashMap<String, u32> = records
        .iter()
        .map(|r| (r.name.clone(), r.length as u32))
        .collect();

    let writer = BigWigWrite::create_file(path, chrom_sizes)?;
    let data = BedParserStreamingIterator::wrap_infallible_iter(intervals(records), false);
    let runtime = tokio::runtime::Builder::new_multi_thread().build()?;
    writer.write(data, runtime)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::scan::PieceCurves;
    use std::path::PathBuf;

    fn scored() -> Vec<ScoredRecord> {
        vec![
            ScoredRecord {
                name: "chr1".to_string(),
                length: 100,
                pieces: vec![PieceCurves {
                    start: 1,
                    curve_start: 22,
                    curves: vec![1.5, 2.5, 3.5],
                }],
            },
            ScoredRecord {
                name: "chr2".to_string(),
                length: 60,
                pieces: vec![PieceCurves {
                    start: 10,
                    curve_start: 31,
                    curves: vec![4.0],
                }],
            },
        ]
    }

    fn tmp(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("symcurve-test-{}-{}", std::process::id(), name));
        p
    }

    #[test]
    fn test_format_from_extension() {
        let cases = [
            ("out.bw", Some(OutputFormat::BigWig)),
            ("out.bigWig", Some(OutputFormat::BigWig)),
            ("OUT.BW", Some(OutputFormat::BigWig)),
            ("out.bedGraph", Some(OutputFormat::BedGraph)),
            ("out.bg", Some(OutputFormat::BedGraph)),
            ("out.txt", None),
            ("out", None),
        ];
        for (name, expected) in cases {
            let got = OutputFormat::from_path(Path::new(name)).ok();
            assert_eq!(got, expected, "for {name}");
        }
    }

    #[test]
    fn test_unknown_extension_message_lists_the_options() {
        let err = OutputFormat::from_path(Path::new("out.txt")).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("out.txt"), "{msg}");
        assert!(msg.contains(".bw") && msg.contains(".bedGraph"), "{msg}");
    }

    #[test]
    fn test_intervals_are_zero_based_half_open_and_in_order() {
        let records = scored();
        let got: Vec<_> = intervals(&records).collect();
        assert_eq!(got.len(), 4);
        // curve_start 22 is 1-based, so the first interval starts at 21 zero-based.
        assert_eq!(got[0].0, "chr1");
        assert_eq!((got[0].1.start, got[0].1.end), (21, 22));
        assert_eq!((got[1].1.start, got[1].1.end), (22, 23));
        assert_eq!((got[2].1.start, got[2].1.end), (23, 24));
        assert_eq!(got[3].0, "chr2");
        assert_eq!((got[3].1.start, got[3].1.end), (30, 31));
    }

    #[test]
    fn test_bedgraph_contents() {
        let path = tmp("out.bedGraph");
        let format = write(&path, &scored()).unwrap();
        assert_eq!(format, OutputFormat::BedGraph);
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "chr1\t21\t22\t1.5");
        assert_eq!(lines[3], "chr2\t30\t31\t4");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_bigwig_round_trips() {
        // Write a bigWig and read it back, rather than checking the magic bytes: a file
        // with the right header but wrong contents would still be useless.
        let path = tmp("out.bw");
        let format = write(&path, &scored()).unwrap();
        assert_eq!(format, OutputFormat::BigWig);

        let mut read = bigtools::BigWigRead::open_file(&path).unwrap();
        let mut chroms: Vec<_> = read
            .chroms()
            .iter()
            .map(|c| (c.name.clone(), c.length))
            .collect();
        chroms.sort();
        assert_eq!(
            chroms,
            vec![("chr1".to_string(), 100u32), ("chr2".to_string(), 60u32)]
        );

        let values: Vec<_> = read.values("chr1", 21, 24).unwrap();
        assert_eq!(values.len(), 3);
        assert!((values[0] - 1.5).abs() < 1e-6, "{values:?}");
        assert!((values[2] - 3.5).abs() < 1e-6, "{values:?}");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_write_refuses_an_unknown_extension_without_creating_a_file() {
        let path = tmp("out.txt");
        assert!(write(&path, &scored()).is_err());
        assert!(
            !path.exists(),
            "no file should be created for a rejected format"
        );
    }
}
