mod cli;
use cli::Cli;
use clap::Parser;

fn main() {
    let cli = Cli::parse();
    let input = cli.input.to_str().unwrap();
    if !input.is_empty() {
        println!("Value for input: {input}");
    }
}
