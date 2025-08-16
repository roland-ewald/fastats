use bed::feature::RecordBuf;
use bed::io::writer::Writer as BedWriter;
use bstr::ByteSlice;
use noodles_bed as bed;
use noodles_core::Position;
use noodles_fasta::Record as FastaRecord;
use serde::Serialize;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::result::Result;

#[derive(Debug, Clone, Serialize)]
pub struct SequenceStatistics {
    pub sequence_name: String,
    pub non_masked_bases: usize,
    pub soft_masked_bases: usize,
    pub hard_masked_bases: usize,
    pub gc_content: f64,
    pub sequence_length: usize,
}

pub fn process_fasta(output_dir: &PathBuf) -> impl Fn(&FastaRecord) -> Option<SequenceStatistics> {
    move |record| process_fasta_record(record, output_dir)
}

fn process_fasta_record(record: &FastaRecord, output_dir: &PathBuf) -> Option<SequenceStatistics> {
    
    let record_name: &str = record.definition().name().to_str().expect(
        format!(
            "Failed to convert record name to string: '{}'",
            record.definition().name()
        )
        .as_str(),
    );

    if record.sequence().len() == 0 {
        return Some(SequenceStatistics {
            sequence_name: record_name.to_string(),
            non_masked_bases: 0,
            soft_masked_bases: 0,
            hard_masked_bases: 0,
            gc_content: 0.0,
            sequence_length: 0,
        });
    }

    let mut non_masked_bed_writer = create_bed_writer(
        output_dir
            .join(&format!("{}-non-masked.bed", record_name))
            .to_str()
            .unwrap(),
    );
    let mut soft_masked_bed_writer = create_bed_writer(
        output_dir
            .join(&format!("{}-soft-masked.bed", record_name))
            .to_str()
            .unwrap(),
    );
    let mut hard_masked_bed_writer = create_bed_writer(
        output_dir
            .join(&format!("{}-hard-masked.bed", record_name))
            .to_str()
            .unwrap(),
    );

    let sequence: &[u8] = record.sequence().as_ref();

    let mut index1: usize = 0;
    let mut gc_counter: usize = 0;

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
            }
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
            }
            b'N' => {
                hard_mask_counter += 1;
                hard_masking = true;
            }
            b'n' => {
                soft_mask_counter += 1;
                hard_mask_counter += 1;
                soft_masking = true;
                hard_masking = true;
            }
            _ => panic!("Unexpected base: '{}'", *base as char),
        }

        update_mask_region(
            &mut non_mask_region_start1,
            non_masking,
            &mut non_masked_bed_writer,
            record_name,
            index1,
        );
        update_mask_region(
            &mut soft_masked_region_start1,
            soft_masking,
            &mut soft_masked_bed_writer,
            record_name,
            index1,
        );
        update_mask_region(
            &mut hard_masked_region_start1,
            hard_masking,
            &mut hard_masked_bed_writer,
            record_name,
            index1,
        );
    }
    assert!(
        non_mask_counter + soft_mask_counter + hard_mask_counter >= sequence.len(), // '>=' because 'n' counts as both soft and hard masked
        "The sum of masked bases does not match the sequence length ({}) for '{}'. This seems to be a bug.",
        sequence.len(),
        record_name
    );
    return Some(SequenceStatistics {
        sequence_name: record_name.to_string(),
        non_masked_bases: non_mask_counter,
        soft_masked_bases: soft_mask_counter,
        hard_masked_bases: hard_mask_counter,
        gc_content: gc_counter as f64 / sequence.len() as f64,
        sequence_length: sequence.len(),
    });
}

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
    index1: usize,
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
            let mut writer: BedWriter<3, BufWriter<File>> =
                bed::io::writer::Builder::default().build_from_path("test.bed")?;
            write_bed_record(&mut writer, "chr2", 123, 124)?;
        }
        Ok(())
    }
}
