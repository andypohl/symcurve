//! The symcurve command line tool.

use std::error::Error;

use clap::Parser;

use symcurve::cli::Cli;
use symcurve::curve::matrix::RollType;
use symcurve::curve::scan::{CurveParams, DEFAULT_CHUNK_SCORES};
use symcurve::output::{self, OutputFormat};
use symcurve::stream;

/// Upper bound on how much sequence to hold when reading through an index. Past this a
/// larger window costs memory without saving meaningful work.
const MAX_WINDOW_BASES: usize = 64 << 20;

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

    let threads = rayon::current_num_threads();
    let chunk_scores = args.max_memory.chunk_scores(threads, DEFAULT_CHUNK_SCORES);
    if args.verbose {
        eprintln!(
            "budget {} over {threads} threads: {chunk_scores} scores per chunk",
            args.max_memory
        );
    }

    // An index lets each record be read a window at a time, which lowers the floor set
    // by the largest record, and supplies chromosome sizes without a pass over the file.
    let index = stream::load_index(&args.input)?;
    let window_bases = args.max_memory.window_bases(MAX_WINDOW_BASES);
    if args.verbose {
        match &index {
            Some(_) => eprintln!(
                "using {}.fai: reading in windows of {window_bases} bases",
                args.input.display()
            ),
            None => eprintln!(
                "no {}.fai found: reading whole records, so the largest record sets the floor",
                args.input.display()
            ),
        }
    }

    let input = args.input.clone();
    let verbose = args.verbose;
    let index_for_run = index.clone();
    let produce = move |emit: &mut dyn FnMut(&str, usize, f64) -> std::io::Result<()>| {
        let stats = match &index_for_run {
            Some(index) => stream::for_each_score_indexed(
                &input,
                index,
                &params,
                chunk_scores,
                threads,
                window_bases,
                emit,
            )?,
            None => stream::for_each_score(&input, &params, chunk_scores, threads, emit)?,
        };
        if verbose {
            eprintln!(
                "{} records, {} scores; largest record {} bases",
                stats.records, stats.scores, stats.longest_record
            );
        }
        Ok(())
    };

    match format {
        OutputFormat::BedGraph => output::write_bedgraph_streaming(&args.output, produce)?,
        OutputFormat::BigWig => {
            // The header needs every chromosome size before any value, so names and
            // lengths are read in a first pass that keeps no sequence.
            let sizes = match &index {
                // The index already carries every name and length.
                Some(index) => stream::chrom_sizes_from_index(index),
                None => stream::chrom_sizes(&args.input)
                    .map_err(|e| format!("cannot read {}: {e}", args.input.display()))?,
            };
            // Keep roughly a chunk's worth of values queued: enough to keep the writer
            // fed, bounded so a slow writer cannot let the queue grow without limit.
            let queue_depth = chunk_scores.min(1 << 16);
            output::write_bigwig_streaming(&args.output, sizes, queue_depth, produce)?
        }
    }

    if args.verbose {
        eprintln!(
            "wrote {} as {}",
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
