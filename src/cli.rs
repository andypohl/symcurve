//! Command line interface symcurve tool.
//!
//! The main arguments to the symcurve CLI are the input and output file paths.
//! These should be provided as positional arguments. The other arguments are optional
//! but have constraints and default values.

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
    #[arg(short, long, name = "YAML")]
    pub matrices: Option<PathBuf>,

    /// curve step
    #[arg(long, default_value = "15", value_parser = clap::value_parser!(u16).range(1..), name = "INT")]
    pub curve_step: u16,
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
            ContextValue::String("--matrices <YAML>".to_owned()),
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
            .starts_with("error: invalid value '0' for '--curve-step <INT>'"));
    }
}
