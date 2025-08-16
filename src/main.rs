use bstr::ByteSlice;
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
        "Input file: '{:?}'. Output directory: '{:?}'.\nNow loading the FASTA records...",
        args.fasta_file, args.output_dir
    );

    let mut reader = File::open(args.fasta_file).map(BufReader::new).map(fasta::io::Reader::new)?;
    let records: Vec<FastaRecord> = reader.records().collect::<Result<_, _>>()?;
    
    let _results: Vec<()> = records.par_iter().map(|record| {

        let record_name: &str = record.definition().name().to_str()
            .expect(format!("Failed to convert record name to string: '{}'", record.definition().name()).as_str());
        
        if record.sequence().len() == 0 {
            println!("Skipping empty sequence: {}", record.definition().name());
            return;
        }

        let mut non_masked_bed_writer = create_bed_writer(&format!("{}-non-masked.bed", record_name));
        let mut soft_masked_bed_writer = create_bed_writer(&format!("{}-soft-masked.bed", record_name));
        let mut hard_masked_bed_writer = create_bed_writer(&format!("{}-hard-masked.bed", record_name));

        let sequence: &[u8] = record.sequence().as_ref();
        
        let mut index1: usize = 0;
        let mut gc_counter:usize = 0;

        let mut non_mask_counter: usize = 0;
        let mut soft_mask_counter: usize = 0;
        let mut hard_mask_counter: usize = 0;

        let mut non_mask_region_start1: Option<usize> = None;
        let mut soft_masked_region_start1: Option<usize> = None;
        let mut hard_masked_region_start1: Option<usize> = None;

        for base in sequence {
            index1 += 1;
            let mut non_masking: bool = false;
            let mut soft_masking: bool = false;
            let mut hard_masking: bool = false;
            match *base {
                b'C' | b'G' => {                    
                    gc_counter += 1;
                    non_mask_counter += 1;
                    non_masking = true;
                },
                b'c' | b'g' => {
                    gc_counter += 1;
                    soft_mask_counter += 1;
                    soft_masking = true;
                }
                b'A' | b'T' => {
                    non_mask_counter += 1;
                    non_masking = true;
                }
                b'a' | b't' => {
                    soft_mask_counter += 1;
                    soft_masking = true;
                },
                b'N' => {
                    hard_mask_counter += 1;
                    hard_masking = true;
                },
                b'n' => {
                    soft_mask_counter += 1;
                    hard_mask_counter += 1;
                    soft_masking = true;
                    hard_masking = true;
                },
                _ => panic!("Unexpected base: '{}'", *base as char),
            }

            update_mask_region(&mut non_mask_region_start1, non_masking, &mut non_masked_bed_writer, record_name, index1);
            update_mask_region(&mut soft_masked_region_start1, soft_masking, &mut soft_masked_bed_writer, record_name, index1);
            update_mask_region(&mut hard_masked_region_start1, hard_masking, &mut hard_masked_bed_writer, record_name, index1);
        }
        assert!(non_mask_counter + soft_mask_counter + hard_mask_counter == sequence.len(),
            "The sum of masked bases does not match the sequence length ({}) for '{}'. This seems to be a bug", sequence.len(), record_name);

        println!(
            "Found {} hard-masked bases, {} soft-masked bases, {} non-masked bases in sequence '{}' (GC ratio: {:.2}, overall length: {})",
            hard_mask_counter,
            soft_mask_counter,
            non_mask_counter,
            record_name,
            gc_counter as f64 / sequence.len() as f64,
            sequence.len()
        );
    }).collect();
    // TODO: store results in BED files AND a JSON summary file (so that eg jq can be used downstream)
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
