use clap::Parser;
use fastats::*;
use noodles_fasta as fasta;
use noodles_fasta::Record as FastaRecord;
use rayon::prelude::*;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::io::ErrorKind;
use std::path::PathBuf;
use std::result::Result;

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

    #[arg(
        short = 'q',
        long = "quiet",
        default_value = "false",
        help = "Do not print results on stdout."
    )]
    quiet: bool,

    #[arg(
        long = "ignore-iupac",
        default_value = "false",
        help = "Enable this to avoid failing when encountering a sequence character that is not in ('A', 'C', 'T', 'G', 'N', 'a', 'c', 't', 'g', 'n')."
    )]
    ignore_iupac: bool,

    #[arg(
        long = "no-bed-output",
        default_value = "false",
        help = "Do not store masking regions into BED files."
    )]
    no_bed_output: bool,

    #[arg(
        long = "match-regex",
        default_value = ".*",
        help = "Regular expression to focus the analysis on sequences matching a specific regular expression."
    )]
    sequence_match_regex: String,

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
            fs::create_dir_all(&self.output_dir)
        } else {
            Ok(())
        }
    }

    fn bed_output_dir(&self) -> Option<&PathBuf> {
        if self.no_bed_output {
            None
        } else {
            Some(&self.output_dir)
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    args.validate().expect("Failed to validate CLI arguments");

    let mut reader = File::open(&args.fasta_file)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)?;
    let records: Vec<FastaRecord> = reader.records().collect::<Result<_, _>>()?;

    let mut sequence_statistics: Vec<SequenceStatistics> = records
        .par_iter()
        .flat_map(process_fasta(args.bed_output_dir().map(|pb| pb.as_path()), args.sequence_match_regex.as_str(), args.ignore_iupac))
        .collect();
    sequence_statistics.sort_by_key(|s| s.sequence_name.clone());

    let json_output = serde_json::to_string_pretty(&sequence_statistics).unwrap();
    if !args.quiet {
        println!("{}", json_output.clone());
    }
    fs::write(args.output_dir.join("summary.json"), json_output)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_validation() {
        let cli = Cli {
            fasta_file: PathBuf::from("does-not-exist.fasta"),
            output_dir: PathBuf::from("output"),
            quiet: false,
            ignore_iupac: false,
            no_bed_output: false,
            sequence_match_regex: ".*".to_string(),
        };
        // Test invalid input file
        assert!(cli.validate().is_err());
    }
}
