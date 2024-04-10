mod cli;
use cli::Cli;
use clap::Parser;

// still basically a hello-world
fn main() {
    let cli = Cli::parse();
    let input = cli.input.to_str().unwrap();
    if !input.is_empty() {
        println!("Value for input: {input}");
    }
}
