use bed::feature::RecordBuf;
use bed::io::writer::Writer as BedWriter;
use clap::Parser;
use noodles_bed as bed;
use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_fasta::Record as FastaRecord;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::{fs, io::BufReader, io::ErrorKind, path::PathBuf, result::Result};
use rayon::prelude::*;

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
        let record_name = record.definition().name();
        if record.sequence().len() == 0 {
            println!("Skipping empty sequence: {}", record.definition().name());
            return;
        }

        let soft_masked_bed_path = format!("{}-soft-masked.bed", record_name);
        let mut soft_masked_bed_writer: BedWriter<3, BufWriter<File>> = bed::io::writer::Builder::default()
        .build_from_path(soft_masked_bed_path.clone())
        .expect(format!("Could not write to output BED file '{}'.", soft_masked_bed_path).as_str());

        let hard_masked_bed_path = format!("{}-hard-masked.bed", record_name);
        let mut hard_masked_bed_writer: BedWriter<3, BufWriter<File>> = bed::io::writer::Builder::default()
        .build_from_path(hard_masked_bed_path.clone())
        .expect(format!("Could not write to output BED file '{}'.", hard_masked_bed_path).as_str());

        let sequence: &[u8] = record.sequence().as_ref();
        let mut hard_mask_counter:usize = 0;
        let mut soft_mask_counter:usize = 0;
        let mut gc_counter:usize = 0;
        let mut index1: usize = 0;
        let mut soft_masked_region_start1:Option<usize> = None;
        let mut hard_masked_region_start1:Option<usize> = None;

        for base in sequence {
            index1 += 1;
            let mut soft_masking = false;
            let mut hard_masking = false;
            match *base {
                b'C' | b'G' => {
                    gc_counter += 1;
                },
                b'c' | b'g' => {
                    gc_counter += 1;
                    soft_mask_counter += 1;
                    soft_masking = true;
                }
                b'A' | b'T' => {
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

            if soft_masking {
                if soft_masked_region_start1.is_none() {
                    soft_masked_region_start1 = Some(index1);
                }
            } else if let Some(start1) = soft_masked_region_start1 {
                let _result = write_bed_record( &mut soft_masked_bed_writer,
                    record_name.to_string().as_str(),
                    start1,
                    index1,
                );
                soft_masked_region_start1 = None;
            }

            if hard_masking {
                if hard_masked_region_start1.is_none() {
                    hard_masked_region_start1 = Some(index1);
                }
            } else if let Some(start1) = hard_masked_region_start1 {
                let _result = write_bed_record( &mut hard_masked_bed_writer,
                    record_name.to_string().as_str(),
                    start1,
                    index1,
                );
                hard_masked_region_start1 = None;
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
    }).collect();
    // TODO: store results in BED files AND a JSON summary file (so that eg jq can be used downstream)
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

    #[test]
    fn write_bed_record_ok() -> Result<(), Box<dyn Error>> {
        {
        let mut writer: BedWriter<3, BufWriter<File>> = bed::io::writer::Builder::default().build_from_path("test.bed")?;
        write_bed_record(&mut writer, "chr2", 123, 124)?;
        }
    Ok(())
    }
}
