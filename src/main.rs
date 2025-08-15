use clap::Parser;
use std::{fs, io::Error, io::ErrorKind, path::PathBuf, result::Result};

#[derive(Parser)]
struct Cli {
    #[arg(
        short = 'i',
        long = "input-file",
        help = "The input FASTA file to analyze."
    )]
    input_file: PathBuf,

    #[arg(
        short = 'o',
        long = "output-dir",
        default_value = "./output",
        help = "The output directory for the BED and summary files."
    )]
    output_dir: PathBuf,
}

impl Cli {
    fn validate(self: &Cli) -> Result<(), Error> {
        if !self.input_file.is_file() {
            Err(Error::new(
                ErrorKind::InvalidInput,
                format!("The input file '{:?}' is not a file.", self.input_file),
            ))
        } else if self.output_dir.is_file() {
            Err(Error::new(
                ErrorKind::InvalidInput,
                format!("The output directory '{:?}' is a file.", self.output_dir),
            ))
        } else if !self.output_dir.exists() {
            println!(
                "Creating output directory '{:?}', as it does not yet exist.",
                self.output_dir
            );
            fs::create_dir_all(&self.output_dir)
        } else {
            Ok(())
        }
    }
}

fn main() {
    let args = Cli::parse();
    let validation_result = args.validate();
    if validation_result.is_err() {
        eprintln!(
            "Stopping, as input parameters are invalid: '{:?}'.",
            validation_result.err()
        );
    } else {
        println!(
            "Input file: {:?}, Output directory: {:?}",
            args.input_file, args.output_dir
        );
        println!("Done.");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_validation() {
        let cli = Cli {
            input_file: PathBuf::from("test.fasta"),
            output_dir: PathBuf::from("output"),
        };
        // Test invalid input file
        assert!(cli.validate().is_err());
    }
}

