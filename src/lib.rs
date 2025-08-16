use bed::feature::RecordBuf;
use bed::io::writer::Writer as BedWriter;
use bstr::ByteSlice;
use clap::Parser;
use noodles_bed as bed;
use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_fasta::Record as FastaRecord;
use rayon::prelude::*;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::ErrorKind; 
use std::path::PathBuf;
use std::result::Result;


pub fn create_bed_writer(path: &str) -> BedWriter<3, BufWriter<File>> {
    bed::io::writer::Builder::default()
        .build_from_path(path)
        .expect(&format!("Could not write to output BED file '{}'.", path))
}

pub fn write_bed_record<X: std::io::Write>(
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

pub fn update_mask_region<X: std::io::Write>(
    region_start: &mut Option<usize>,
    masking: bool,
    writer: &mut BedWriter<3, X>,
    record_name: &str,
    index1: usize
) {
    if masking {
        if region_start.is_none() {
            *region_start = Some(index1);
        }
    } else if let Some(start1) = *region_start {
        let _ = write_bed_record(writer, record_name, start1, index1);
        *region_start = None;
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn write_bed_record_ok() -> Result<(), Box<dyn Error>> {
        {
        let mut writer: BedWriter<3, BufWriter<File>> = bed::io::writer::Builder::default().build_from_path("test.bed")?;
        write_bed_record(&mut writer, "chr2", 123, 124)?;
        }
    Ok(())
    }
}