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
//!   <OUTPUT>  bigWig output file path
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
//!   -h, --help                               Print help
//!   -V, --version                            Print version
//! ```

use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version = env!("CARGO_PKG_VERSION"), about = "Symmetry of DNA curvature.", long_about = None)]
pub struct Cli {
    /// FASTA input file path
    pub input: PathBuf,

    /// bigWig output file path
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
    pub curve_scale: f32,

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
}

fn parse_float_in_range(s: &str) -> Result<f32, String> {
    let value = s
        .parse::<f32>()
        .map_err(|_| "Value must be a floating-point number")?;
    if value >= 0.0 && value <= 1.0 {
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
        let args = Cli::parse_from(&[
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
        assert_eq!(args.verbose, true);
        assert_eq!(args.matrices.unwrap().to_str().unwrap(), "matrices.yaml");
        assert_eq!(args.curve_step, 20);
    }

    #[test]
    fn test_missing_matrix_file() {
        let args_result =
            Cli::try_parse_from(&["symcurve", "input.fasta", "output.bw", "--matrices"]);
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
        assert_eq!(args_result.is_err(), true);
        assert_eq!(args_result.unwrap_err().to_string(), err.to_string());
    }

    // test when the curve_step argument is not > 0
    #[test]
    fn test_zero_curve_step() {
        let args_result =
            Cli::try_parse_from(&["symcurve", "input.fasta", "output.bw", "--curve-step", "0"]);
        assert_eq!(args_result.is_err(), true);
        assert!(args_result
            .unwrap_err()
            .to_string()
            .starts_with("error: invalid value '0' for '--curve-step"));
    }

    // helper to test_curve_scale()
    fn get_different_curve_scale_parsings(curve_scale_s: &str) -> Result<Cli, clap::error::Error> {
        return Cli::try_parse_from(&[
            "symcurve",
            "input.fasta",
            "output.bw",
            "--curve-scale",
            curve_scale_s,
        ]);
    }

    #[test]
    fn test_curve_scale() {
        // test different passed in curve scales
        assert_eq!(get_different_curve_scale_parsings("0").is_ok(), true);
        assert_eq!(get_different_curve_scale_parsings("0.33").is_ok(), true);
        assert_eq!(get_different_curve_scale_parsings("1").is_ok(), true);
        assert_eq!(get_different_curve_scale_parsings("1.1").is_err(), true);
        assert_eq!(get_different_curve_scale_parsings("-1").is_err(), true);
        assert_eq!(get_different_curve_scale_parsings("abc").is_err(), true);
    }
}
