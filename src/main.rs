use clap::Parser;
use noodles_fasta as fasta;
use std::error::Error;
use std::fs::File;
use std::{fs, io::BufReader, io::ErrorKind, path::PathBuf, result::Result};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    fasta_file: PathBuf,

    #[arg(
        short = 'o',
        long = "output-dir",
        default_value = ".",
        help = "The output directory for the BED and summary files."
    )]
    output_dir: PathBuf,
}

impl Cli {
    fn validate(self: &Cli) -> Result<(), std::io::Error> {
        if !self.fasta_file.is_file() {
            Err(std::io::Error::new(
                ErrorKind::InvalidInput,
                format!("The input file '{:?}' is not a file.", self.fasta_file),
            ))
        } else if self.output_dir.is_file() {
            Err(std::io::Error::new(
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

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    args.validate().expect("Failed to validate CLI arguments");
    println!(
        "Input file: {:?}, Output directory: {:?}",
        args.fasta_file, args.output_dir
    );

    let mut reader = File::open(args.fasta_file)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)?;

    for result in reader.records() {
        let record = result?;
        let record_name = record.definition().name();
        println!(
            "Processing record: {:?} with length {:?}",
            record.definition().name(),
            record.sequence().len()
        );

        let sequence: &[u8] = record.sequence().as_ref();
        let mut hard_mask_counter= 0;
        for base in sequence {
            if *base == b'N' {
                hard_mask_counter += 1;
            }
        }
        println!("Found {} 'N' bases in sequence: {}", hard_mask_counter, record_name);
    }
    println!("Done.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_validation() {
        let cli = Cli {
            fasta_file: PathBuf::from("test.fasta"),
            output_dir: PathBuf::from("output"),
        };
        // Test invalid input file
        assert!(cli.validate().is_err());
    }
}
