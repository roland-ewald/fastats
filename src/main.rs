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

    // 1) quiet mode, 2) contig regex match, 3) json output file
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
    //TODO: print only if JSON output file given and not in quiet mode
/*     println!(
        "Input file: '{:?}'. Output directory: '{:?}'.\nNow loading the FASTA records...",
        args.fasta_file, args.output_dir
    ); */

    let mut reader = File::open(args.fasta_file).map(BufReader::new).map(fasta::io::Reader::new)?;
    let records: Vec<FastaRecord> = reader.records().collect::<Result<_, _>>()?;
    
    let record_processing = process_fasta(&args.output_dir);
    let mut sequence_statistics: Vec<SequenceStatistics> = records.par_iter().flat_map(record_processing).collect();
    sequence_statistics.sort_by_key(|s| s.sequence_name.clone());

    let json_output = serde_json::to_string_pretty(&sequence_statistics).unwrap();

    // TODO: store results a JSON summary file OR print to stdout (so that eg jq can be used downstream)
    println!("{}", json_output);
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
