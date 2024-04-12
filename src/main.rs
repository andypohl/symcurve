mod cli;
use cli::Cli;
use clap::Parser;

// still basically a hello-world
fn main() {
    Cli::parse();
}
