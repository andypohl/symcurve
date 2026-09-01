//! Writing curvature scores out, in the format implied by the output path.

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;

use bigtools::BigWigWrite;
use bigtools::Value;

use crate::curve::calls::NucleosomeCall;
use bigtools::beddata::BedParserStreamingIterator;

/// The output formats this tool can write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// Indexed binary bigWig, via bigtools.
    BigWig,
    /// Plain-text bedGraph: `chrom<TAB>start<TAB>end<TAB>value`, 0-based half-open.
    BedGraph,
    /// Plain-text GFF of nucleosome calls, one feature per call.
    Gff,
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
            "gff" | "gff2" | "gff3" => Ok(OutputFormat::Gff),
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

/// A score at a position, as the writers want it.
///
/// Positions arrive 1-based and inclusive from the curve code and are converted here to
/// the 0-based half-open convention both output formats use.
fn interval(name: &str, position: usize, score: f64) -> (String, Value) {
    let start = (position - 1) as u32;
    (
        name.to_string(),
        Value {
            start,
            end: start + 1,
            value: score as f32,
        },
    )
}

/// Write a bedGraph, pulling scores from `produce` as they are needed.
///
/// Nothing is retained: each score is formatted and written as it arrives.
pub fn write_bedgraph_streaming<P>(
    path: &Path,
    produce: P,
) -> Result<(), Box<dyn std::error::Error>>
where
    P: FnOnce(&mut dyn FnMut(&str, usize, f64) -> io::Result<()>) -> io::Result<()>,
{
    let mut out = BufWriter::new(File::create(path)?);
    produce(&mut |name, position, score| {
        let (chrom, value) = interval(name, position, score);
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            chrom, value.start, value.end, value.value
        )
    })?;
    out.flush()?;
    Ok(())
}

/// How many values travel through the channel at once.
///
/// Sending a hundred million values one at a time costs more in synchronization than in
/// work, so they move in batches. Memory stays bounded because both the batch size and
/// the number of batches in flight are fixed.
const SEND_BATCH: usize = 8192;

/// Write a bigWig, pulling scores from `produce` as they are needed.
///
/// bigtools takes an iterator, but the scoring loop pushes, so the producer runs on its
/// own thread and feeds a bounded channel. The bound is what makes this streaming rather
/// than merely deferred: when the writer falls behind, the producer blocks instead of
/// queueing the genome up in memory.
pub fn write_bigwig_streaming<P>(
    path: &Path,
    chrom_sizes: Vec<(String, u32)>,
    queue_depth: usize,
    produce: P,
) -> Result<(), Box<dyn std::error::Error>>
where
    P: FnOnce(&mut dyn FnMut(&str, usize, f64) -> io::Result<()>) -> io::Result<()>
        + Send
        + 'static,
{
    let sizes: HashMap<String, u32> = chrom_sizes.into_iter().collect();
    let writer = BigWigWrite::create_file(path, sizes)?;

    let batches_in_flight = queue_depth.div_ceil(SEND_BATCH).max(1);
    let (tx, rx) = mpsc::sync_channel::<Vec<(String, Value)>>(batches_in_flight);
    // The producer's error is reported after the write finishes: a send failure only
    // says the consumer went away, not why, and a read error must not look like a
    // successful but truncated file.
    let failure = Arc::new(Mutex::new(None::<String>));
    let producer_failure = Arc::clone(&failure);

    let producer = thread::spawn(move || {
        let mut batch: Vec<(String, Value)> = Vec::with_capacity(SEND_BATCH);
        let result = (|| {
            produce(&mut |name, position, score| {
                batch.push(interval(name, position, score));
                if batch.len() == SEND_BATCH {
                    let full = std::mem::replace(&mut batch, Vec::with_capacity(SEND_BATCH));
                    // A closed channel means the writer stopped; report it so a partial
                    // file is never mistaken for a complete one.
                    tx.send(full)
                        .map_err(|_| io::Error::other("bigWig writer stopped accepting values"))?;
                }
                Ok(())
            })?;
            if !batch.is_empty() {
                tx.send(std::mem::take(&mut batch))
                    .map_err(|_| io::Error::other("bigWig writer stopped accepting values"))?;
            }
            Ok::<(), io::Error>(())
        })();
        if let Err(err) = result {
            *producer_failure.lock().unwrap() = Some(err.to_string());
        }
    });

    let data = BedParserStreamingIterator::wrap_infallible_iter(rx.into_iter().flatten(), false);
    let runtime = tokio::runtime::Builder::new_multi_thread().build()?;
    let write_result = writer.write(data, runtime);

    producer.join().map_err(|_| "scoring thread panicked")?;
    if let Some(err) = failure.lock().unwrap().take() {
        return Err(err.into());
    }
    write_result?;
    Ok(())
}

/// Write nucleosome calls as GFF, pulling records from `produce` as they are needed.
///
/// One feature per call. The attribute column carries the called sequence, as the
/// reference implementation does.
///
/// Coordinates follow the reference rather than the GFF specification. It prints
/// `dyad - half_width` directly, which is a zero-based index into the record, where GFF
/// expects one-based inclusive coordinates. Everything it emits is therefore one base to
/// the left of where a browser will read it. Reproduced so that positions can be compared
/// against the reference; see the Algorithm Issues page.
pub fn write_gff_streaming<P>(
    path: &Path,
    feature: &str,
    produce: P,
) -> Result<(), Box<dyn std::error::Error>>
where
    P: FnOnce(
        &mut dyn FnMut(&str, &[u8], &[NucleosomeCall], usize) -> io::Result<()>,
    ) -> io::Result<()>,
{
    let mut out = BufWriter::new(File::create(path)?);
    produce(&mut |name, bases, calls, half_width| {
        for call in calls {
            let start = call.dyad - half_width;
            let end = call.dyad + half_width;
            let sequence = std::str::from_utf8(&bases[start..=end]).unwrap_or("");
            writeln!(
                out,
                "{name}\tevidence\t{feature}\t{start}\t{end}\t{}\t+\t.\t{sequence}",
                call.reported_score()
            )?;
        }
        Ok(())
    })?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// A producer emitting a fixed set of scores, standing in for a scoring run.
    fn produce(emit: &mut dyn FnMut(&str, usize, f64) -> io::Result<()>) -> io::Result<()> {
        emit("chr1", 22, 1.5)?;
        emit("chr1", 23, 2.5)?;
        emit("chr1", 24, 3.5)?;
        emit("chr2", 31, 4.0)?;
        Ok(())
    }

    fn sizes() -> Vec<(String, u32)> {
        vec![("chr1".to_string(), 100), ("chr2".to_string(), 60)]
    }

    fn tmp(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("symcurve-out-{}-{}", std::process::id(), name));
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
            assert_eq!(
                OutputFormat::from_path(Path::new(name)).ok(),
                expected,
                "for {name}"
            );
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
    fn test_interval_is_zero_based_half_open() {
        // Position 22 is 1-based, so the interval covers [21, 22).
        let (chrom, value) = interval("chr1", 22, 1.5);
        assert_eq!(chrom, "chr1");
        assert_eq!((value.start, value.end), (21, 22));
        assert!((value.value - 1.5).abs() < 1e-6);
    }

    #[test]
    fn test_bedgraph_contents() {
        let path = tmp("out.bedGraph");
        write_bedgraph_streaming(&path, produce).unwrap();
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
        write_bigwig_streaming(&path, sizes(), 8, produce).unwrap();

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
    fn test_bigwig_survives_a_queue_depth_of_one() {
        // A depth of 1 makes the producer block on nearly every value, which is the
        // backpressure path; it must produce the same file, not deadlock.
        let path = tmp("depth1.bw");
        write_bigwig_streaming(&path, sizes(), 1, produce).unwrap();
        let mut read = bigtools::BigWigRead::open_file(&path).unwrap();
        let values: Vec<_> = read.values("chr1", 21, 24).unwrap();
        assert_eq!(values.len(), 3);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_producer_error_is_reported_not_silently_truncated() {
        // A read failure partway through must not look like a short but successful file.
        let path = tmp("fail.bw");
        let failing = |emit: &mut dyn FnMut(&str, usize, f64) -> io::Result<()>| {
            emit("chr1", 22, 1.5)?;
            Err(io::Error::other("synthetic read failure"))
        };
        let err = write_bigwig_streaming(&path, sizes(), 8, failing).unwrap_err();
        assert!(err.to_string().contains("synthetic read failure"), "{err}");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_bigwig_streaming_flushes_full_batches() {
        // More values than one send batch, so the mid-run flush runs and not only the
        // final drain.
        let n = 2 * SEND_BATCH + 7;
        let path = tmp("batches.bw");
        let many = move |emit: &mut dyn FnMut(&str, usize, f64) -> io::Result<()>| {
            for i in 0..n {
                emit("chr1", i + 1, i as f64)?;
            }
            Ok(())
        };
        write_bigwig_streaming(&path, vec![("chr1".to_string(), n as u32)], 8, many).unwrap();

        let mut read = bigtools::BigWigRead::open_file(&path).unwrap();
        let values: Vec<_> = read.values("chr1", 0, n as u32).unwrap();
        assert_eq!(values.len(), n);
        assert_eq!(values[0], 0.0);
        assert_eq!(values[n - 1], (n - 1) as f32);
        std::fs::remove_file(&path).ok();
    }
}
