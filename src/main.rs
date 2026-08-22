//! The symcurve command line tool.

use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use clap::Parser;

use symcurve::cli::Cli;
use symcurve::curve::matrix::RollType;
use symcurve::curve::scan::{CurveParams, score_pieces};
use symcurve::fasta::split_seq_by_gaps;
use symcurve::output::{self, OutputFormat, ScoredRecord};

fn main() {
    if let Err(err) = run() {
        // Print Display rather than Debug: returning the error from main would render it
        // with Debug, which shows struct internals instead of the message.
        eprintln!("error: {err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();

    // Resolve the output format before doing any work, so an unusable output path fails
    // immediately rather than after scoring an entire genome.
    let format = OutputFormat::from_path(&args.output)?;

    let params = curve_params(&args);
    warn_about_unused_arguments(&args, &params);

    let file = File::open(&args.input)
        .map_err(|e| format!("cannot open {}: {e}", args.input.display()))?;
    let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));

    let mut scored = Vec::new();
    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let length = record.sequence().len();
        let pieces = split_seq_by_gaps(record);
        let curves = score_pieces(&pieces, &params);
        if args.verbose {
            let scores: usize = curves.iter().map(|p| p.curves.len()).sum();
            eprintln!(
                "{name}: {length} bases, {} pieces, {scores} scores",
                pieces.len()
            );
        }
        scored.push(ScoredRecord {
            name,
            length,
            pieces: curves,
        });
    }

    output::write(&args.output, &scored)?;

    if args.verbose {
        let total: usize = scored
            .iter()
            .map(|r| r.pieces.iter().map(|p| p.curves.len()).sum::<usize>())
            .sum();
        eprintln!(
            "wrote {total} scores for {} records to {} as {}",
            scored.len(),
            args.output.display(),
            match format {
                OutputFormat::BigWig => "bigWig",
                OutputFormat::BedGraph => "bedGraph",
            }
        );
    }
    Ok(())
}

/// Translate the CLI arguments into the parameters the curve code takes.
///
/// The reference implementation has two rolling-mean parameters, `stepone` and `steptwo`,
/// but its weighting only works out when `stepone == steptwo + 2`, and the Rust rolling
/// mean derives both ends from a single `step_b`. So `curve_step_one` is the source of
/// truth, giving `step_b = curve_step_one - 1`, and `curve_step_two` is checked for
/// consistency rather than used.
fn curve_params(args: &Cli) -> CurveParams {
    CurveParams {
        roll_type: RollType::Simple,
        step_b: usize::from(args.curve_step_one) - 1,
        step_c: usize::from(args.curve_step),
        curve_scale: f64::from(args.curve_scale),
    }
}

/// Tell the user about arguments that will not affect the result.
///
/// Several flags describe stages that are not implemented yet. Accepting them silently
/// would let someone believe they had changed the output when they had not.
fn warn_about_unused_arguments(args: &Cli, params: &CurveParams) {
    let implied_step_two = params.step_b.saturating_sub(1);
    if usize::from(args.curve_step_two) != implied_step_two {
        eprintln!(
            "warning: --curve-step-two {} is inconsistent with --curve-step-one {} \
             (which implies {}) and is being ignored",
            args.curve_step_two, args.curve_step_one, implied_step_two
        );
    }
    if args.matrices.is_some() {
        eprintln!("warning: --matrices is not implemented yet and is being ignored");
    }
    for (name, is_default) in [
        ("--symcurve-win", args.symcurve_win == 101),
        ("--symcurve-step", args.symcurve_step == 1),
        ("--min-linker-size", args.min_linker_size == 30),
    ] {
        if !is_default {
            eprintln!("warning: {name} applies to the SymCurv stage, which is not implemented yet");
        }
    }
}
