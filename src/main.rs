mod cli;
use clap::Parser;
use cli::Cli;

// still basically a hello-world
fn main() {
    Cli::parse();
}
