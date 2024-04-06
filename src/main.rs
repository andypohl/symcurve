mod cli;
use cli::Cli;
use clap::Parser;

fn main() {
    let cli = Cli::parse();

    // You can check the value provided by positional arguments, or option arguments
    let input = cli.input.to_str().unwrap();
    if !input.is_empty() {
        println!("Value for input: {input}");
    }

    // Continued program logic goes here...
}
