use bed::feature::RecordBuf;
use bed::io::writer::Writer as BedWriter;
use clap::Parser;
use noodles_bed as bed;
use noodles_core::Position;
use noodles_fasta as fasta;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
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

    {
        let mut writer: BedWriter<3, BufWriter<File>> = bed::io::writer::Builder::default()
            .build_from_path(args.output_dir.join("output.bed"))?;
        write_bed_record(&mut writer, "chr2", 123, 124)?;
    }
    let mut reader = File::open(args.fasta_file)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)?;

    for result in reader.records() {
        let record = result?;
        let record_name = record.definition().name();
        if record.sequence().len() == 0 {
            println!("Skipping empty sequence: {}", record.definition().name());
            continue;
        }

        // TODO: run this in parallel
        // TODO: store results in BED files AND a JSON summary file (so that eg jq can be used downstream)

        let sequence: &[u8] = record.sequence().as_ref();
        let mut hard_mask_counter = 0;
        let mut soft_mask_counter = 0;
        let mut gc_counter = 0;

        for base in sequence {
            match *base {
                b'C' | b'G' => {
                    gc_counter += 1;
                },
                b'c' | b'g' => {
                    gc_counter += 1;
                    soft_mask_counter += 1;
                }
                b'A' | b'T' | b'a' | b't' => {                    
                }
                b'N' => {
                    hard_mask_counter += 1;
                },
                b'n' => {
                    soft_mask_counter += 1;
                    hard_mask_counter += 1;
                },
                _ => panic!("Unexpected base: {}", *base as char),
            }
        }

        println!(
            "Found {} 'N' bases, {} soft-masked bases in sequence: {} (GC ratio: {:.2}, overall length: {})",
            hard_mask_counter,
            soft_mask_counter,
            record_name,
            gc_counter as f64 / sequence.len() as f64,
            sequence.len()
        );
    }
    println!("Done.");
    Ok(())
}

fn write_bed_record<X: std::io::Write>(
    writer: &mut BedWriter<3, X>,
    sequence_name: &str,
    start1: usize,
    end1: usize,
) -> Result<(), Box<dyn Error>> {
    let record = RecordBuf::<3>::builder()
        .set_reference_sequence_name(sequence_name)
        .set_feature_start(Position::try_from(start1)?)
        .set_feature_end(Position::try_from(end1)?)
        .build();
    writer.write_feature_record(&record)?;
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
