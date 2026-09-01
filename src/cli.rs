//! Command line interface symcurve tool.
//!
//! The main arguments to the symcurve CLI are the input and output file paths.
//! These should be provided as positional arguments. The other arguments are optional
//! but have constraints and default values.
//!
//! ```text
//! Symmetry of DNA curvature.
//!
//! Usage: symcurve [OPTIONS] <INPUT> <OUTPUT>
//!
//! Arguments:
//!   <INPUT>   FASTA input file path
//!   <OUTPUT>  output file path; the extension picks the format: .bw/.bigWig, .bedGraph/.bg, or .gff
//!
//! Options:
//!   -v, --verbose                            verbose setting
//!   -m, --matrices <MATRICES>                optional matrices YAML file
//!       --curve-step <CURVE_STEP>            curve step [default: 15]
//!       --curve-scale <CURVE_SCALE>          curve scale [default: 0.33335]
//!       --curve-step-one <CURVE_STEP_ONE>    curve step one [default: 6]
//!       --curve-step-two <CURVE_STEP_TWO>    curve step two [default: 4]
//!       --symcurve-win <SYMCURVE_WIN>        symcurve window [default: 101]
//!       --symcurve-step <SYMCURVE_STEP>      symcurve step [default: 1]
//!       --min-linker-size <MIN_LINKER_SIZE>  minimum linker size [default: 30]
//!       --max-memory <MAX_MEMORY>            score buffer budget [default: 8G]
//!       --stage <STAGE>                      curvature, symmetry, calls, final-calls [default: curvature]
//!       --roll <ROLL>                        simple or active [default: simple]
//!   -h, --help                               Print help
//!   -V, --version                            Print version
//! ```

use clap::{Parser, ValueEnum};
use std::path::PathBuf;

/// Which stage of the calculation to write out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Stage {
    /// Curvature values.
    Curvature,
    /// Symmetry of curvature around each dyad.
    Symmetry,
    /// Every nucleosome call, which may overlap one another.
    Calls,
    /// Non-overlapping nucleosome calls, chosen greedily by score.
    FinalCalls,
}

impl Stage {
    /// Whether this stage produces nucleosome calls rather than per-base scores.
    pub fn is_calls(self) -> bool {
        matches!(self, Stage::Calls | Stage::FinalCalls)
    }
}

impl From<Stage> for crate::curve::scan::Stage {
    fn from(stage: Stage) -> Self {
        match stage {
            Stage::Curvature => Self::Curvature,
            // Calls are derived from symmetry, so the scoring stage is the same.
            Stage::Symmetry | Stage::Calls | Stage::FinalCalls => Self::Symmetry,
        }
    }
}

/// Which roll matrix to score with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Roll {
    /// The simple matrix, used for the DNase state in the reference implementation.
    Simple,
    /// The activated matrix, used for the nucleosome state.
    Active,
}

impl From<Roll> for crate::curve::matrix::RollType {
    fn from(roll: Roll) -> Self {
        match roll {
            Roll::Simple => Self::Simple,
            Roll::Active => Self::Active,
        }
    }
}

#[derive(Parser, Debug)]
#[command(version = env!("CARGO_PKG_VERSION"), about = "Symmetry of DNA curvature.", long_about = None)]
pub struct Cli {
    /// FASTA input file path
    pub input: PathBuf,

    /// output file path; the extension picks the format: .bw/.bigWig, .bedGraph/.bg, or .gff
    pub output: PathBuf,

    /// verbose setting
    #[arg(short, long)]
    pub verbose: bool,

    /// optional matrices YAML file
    #[arg(short, long)]
    pub matrices: Option<PathBuf>,

    /// curve step
    #[arg(long, default_value = "15", value_parser = clap::value_parser!(u16).range(1..))]
    pub curve_step: u16,

    /// curve scale
    #[arg(long, default_value = "0.33335", value_parser = parse_float_in_range)]
    pub curve_scale: f64,

    /// curve step one
    #[arg(long, default_value = "6", value_parser = clap::value_parser!(u16).range(1..))]
    pub curve_step_one: u16,

    /// curve step two
    #[arg(long, default_value = "4", value_parser = clap::value_parser!(u16).range(1..))]
    pub curve_step_two: u16,

    /// symcurve window
    #[arg(long, default_value = "101", value_parser = clap::value_parser!(u16).range(1..))]
    pub symcurve_win: u16,

    /// symcurve step
    #[arg(long, default_value = "1", value_parser = clap::value_parser!(u16).range(1..))]
    pub symcurve_step: u16,

    /// minimum linker size
    #[arg(long, default_value = "30", value_parser = clap::value_parser!(u16).range(1..))]
    pub min_linker_size: u16,

    /// upper bound on memory used to buffer scores, e.g. 8G, 512M, 64K
    #[arg(long, default_value = "8G")]
    pub max_memory: crate::memory::MemoryBudget,

    /// which stage to write out
    #[arg(long, value_enum, default_value_t = Stage::Curvature)]
    pub stage: Stage,

    /// which roll matrix to score with
    #[arg(long, value_enum, default_value_t = Roll::Simple)]
    pub roll: Roll,
}

fn parse_float_in_range(s: &str) -> Result<f64, String> {
    let value = s
        .parse::<f64>()
        .map_err(|_| "Value must be a floating-point number")?;
    if (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err("The value must be between 0 and 1".to_owned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::error::*;

    #[test]
    fn test_cli_args() {
        // Your test code will go here
        let args = Cli::parse_from([
            "symcurve",
            "input.fasta",
            "output.bw",
            "--verbose",
            "--matrices",
            "matrices.yaml",
            "--curve-step",
            "20",
        ]);
        assert_eq!(args.input.to_str().unwrap(), "input.fasta");
        assert_eq!(args.output.to_str().unwrap(), "output.bw");
        assert!(args.verbose);
        assert_eq!(args.matrices.unwrap().to_str().unwrap(), "matrices.yaml");
        assert_eq!(args.curve_step, 20);
    }

    #[test]
    fn test_missing_matrix_file() {
        let args_result =
            Cli::try_parse_from(["symcurve", "input.fasta", "output.bw", "--matrices"]);
        // construct the Error object manually, this probably
        // overkill and the error .to_string() is enough but it's
        // just some practice
        let cmd = clap::Command::new("symcurve");
        let mut err = clap::Error::new(ErrorKind::InvalidValue).with_cmd(&cmd);
        err.insert(
            ContextKind::InvalidArg,
            ContextValue::String("--matrices <MATRICES>".to_owned()),
        );
        err.insert(
            ContextKind::InvalidValue,
            ContextValue::String("".to_owned()),
        );
        assert!(args_result.is_err());
        assert_eq!(args_result.unwrap_err().to_string(), err.to_string());
    }

    // test when the curve_step argument is not > 0
    #[test]
    fn test_zero_curve_step() {
        let args_result =
            Cli::try_parse_from(["symcurve", "input.fasta", "output.bw", "--curve-step", "0"]);
        assert!(args_result.is_err());
        assert!(
            args_result
                .unwrap_err()
                .to_string()
                .starts_with("error: invalid value '0' for '--curve-step")
        );
    }

    // helper to test_curve_scale()
    fn get_different_curve_scale_parsings(curve_scale_s: &str) -> Result<Cli, clap::error::Error> {
        Cli::try_parse_from([
            "symcurve",
            "input.fasta",
            "output.bw",
            "--curve-scale",
            curve_scale_s,
        ])
    }

    #[test]
    fn test_curve_scale_keeps_full_precision() {
        // Parsed as f32 this became 0.33334999..., a relative error of 8e-9 that every
        // curvature value carried and that symmetry squared into 1.7e-8. That was enough
        // to move agreement with the reference implementation by three orders.
        let args = Cli::parse_from(["symcurve", "in.fa", "out.bw"]);
        assert_eq!(args.curve_scale, 0.33335_f64);

        let args = Cli::parse_from(["symcurve", "in.fa", "out.bw", "--curve-scale", "0.1"]);
        assert_eq!(args.curve_scale, 0.1_f64);
        assert_ne!(args.curve_scale, f64::from(0.1_f32));
    }

    #[test]
    fn test_max_memory_parses_and_defaults() {
        let args = Cli::parse_from(["symcurve", "in.fa", "out.bw"]);
        assert_eq!(args.max_memory.to_string(), "8G");

        let args = Cli::parse_from(["symcurve", "in.fa", "out.bw", "--max-memory", "512M"]);
        assert_eq!(args.max_memory.bytes(), 512 * 1024 * 1024);

        let bad = Cli::try_parse_from(["symcurve", "in.fa", "out.bw", "--max-memory", "lots"]);
        assert!(bad.is_err());
        assert!(
            bad.unwrap_err().to_string().contains("8G"),
            "the error should show the expected form"
        );
    }

    #[test]
    fn test_curve_scale() {
        // test different passed in curve scales
        assert!(get_different_curve_scale_parsings("0").is_ok());
        assert!(get_different_curve_scale_parsings("0.33").is_ok());
        assert!(get_different_curve_scale_parsings("1").is_ok());
        assert!(get_different_curve_scale_parsings("1.1").is_err());
        assert!(get_different_curve_scale_parsings("-1").is_err());
        assert!(get_different_curve_scale_parsings("abc").is_err());
    }
}
