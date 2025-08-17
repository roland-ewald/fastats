use bed::feature::RecordBuf;
use bed::io::writer::Writer as BedWriter;
use bstr::ByteSlice;
use noodles_bed as bed;
use noodles_core::Position;
use noodles_fasta::Record as FastaRecord;
use regex::Regex;
use serde::Serialize;
use sha2::Digest;
use sha2::Sha256;
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::BufWriter;
use std::result::Result;

#[derive(Debug, Clone, Serialize)]
pub struct SequenceStatistics {
    pub sequence_name: String,
    pub non_masked_bases: usize,
    pub soft_masked_bases: usize,
    pub hard_masked_bases: usize,
    pub non_masked_ratio: f64,
    pub soft_masked_ratio: f64,
    pub hard_masked_ratio: f64,
    pub gc_content: f64,
    pub sequence_length: usize,
    pub checksum_sha256: String,
}

pub fn process_fasta(
    output_dir: Option<&Path>,
    sequence_match_regex: &str,
) -> impl Fn(&FastaRecord) -> Option<SequenceStatistics> {
    move |record| process_fasta_record(record, output_dir, sequence_match_regex)
}

fn process_fasta_record(
    record: &FastaRecord,
    output_dir: Option<&Path>,
    sequence_match_regex: &str,
) -> Option<SequenceStatistics> {
    let record_name: &str = record.definition().name().to_str().expect(
        format!(
            "Failed to convert record name to string: '{}'",
            record.definition().name()
        )
        .as_str(),
    );

    // Ignore records that do not match the regex
    let regex_matcher = Regex::new(ensure_full_match_regex(sequence_match_regex).as_str());
    if regex_matcher.is_err() {
        panic!("Invalid regular expression: '{}'", sequence_match_regex);
    } else if !regex_matcher.unwrap().is_match(record_name) {
        return None;
    }

    // Report empty sequences with all statistics set to zero.
    if record.sequence().len() == 0 {
        return Some(SequenceStatistics {
            sequence_name: record_name.to_string(),
            non_masked_bases: 0,
            soft_masked_bases: 0,
            hard_masked_bases: 0,
            non_masked_ratio: 0.0,
            soft_masked_ratio: 0.0,
            hard_masked_ratio: 0.0,
            gc_content: 0.0,
            sequence_length: 0,
            checksum_sha256: "".to_string(),
        });
    }

    let mut sha256_hasher = Sha256::new();

    let mut non_masked_bed_writer = create_bed_writer(output_dir, "non-masked", record_name);
    let mut soft_masked_bed_writer = create_bed_writer(output_dir, "soft-masked", record_name);
    let mut hard_masked_bed_writer = create_bed_writer(output_dir, "hard-masked", record_name);

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

        sha256_hasher.update([*base]);
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
            b'n' | b'N' => {
                hard_mask_counter += 1;
                hard_masking = true;
            }
            _ => panic!("Unexpected base: '{}' in sequence '{}'.", *base as char, record_name),
        }

        update_mask_region(&mut non_mask_region_start1, non_masking, non_masked_bed_writer.as_mut(), record_name, index1);
        update_mask_region(&mut soft_masked_region_start1, soft_masking, soft_masked_bed_writer.as_mut(), record_name, index1);
        update_mask_region(&mut hard_masked_region_start1, hard_masking, hard_masked_bed_writer.as_mut(), record_name, index1);
    }

    // Write the last regions if they were not closed yet (with index of base after last base).
    update_mask_region(&mut non_mask_region_start1, false, non_masked_bed_writer.as_mut(), record_name, index1 + 1);
    update_mask_region(&mut soft_masked_region_start1, false, soft_masked_bed_writer.as_mut(), record_name, index1 + 1);
    update_mask_region(&mut hard_masked_region_start1, false, hard_masked_bed_writer.as_mut(), record_name, index1 + 1);

    assert!(
        non_mask_counter + soft_mask_counter + hard_mask_counter == sequence.len(),
        "The sum of masked bases does not match the sequence length ({}) for '{}'. This seems to be a bug.",
        sequence.len(),
        record_name
    );
    return Some(SequenceStatistics {
        sequence_name: record_name.to_string(),
        non_masked_bases: non_mask_counter,
        soft_masked_bases: soft_mask_counter,
        hard_masked_bases: hard_mask_counter,
        non_masked_ratio: non_mask_counter as f64 / sequence.len() as f64,
        soft_masked_ratio: soft_mask_counter as f64 / sequence.len() as f64,
        hard_masked_ratio: hard_mask_counter as f64 / sequence.len() as f64,
        gc_content: gc_counter as f64 / sequence.len() as f64,
        sequence_length: sequence.len(),
        checksum_sha256: format!("{:x}", sha256_hasher.finalize()),
    });
}

fn create_bed_writer(
    output_dir: Option<&Path>,
    bed_ending: &str,
    record_name: &str,
) -> Option<BedWriter<3, BufWriter<File>>> {
    output_dir.and_then(|output_dir| {
        let output_path = output_dir.join(&format!("{}.{}.bed", record_name, bed_ending));
        Some(
            bed::io::writer::Builder::default()
                .build_from_path(output_path.clone())
                .expect(&format!(
                    "Could not write to output BED file '{}'.",
                    output_path.to_str().unwrap()
                )),
        )
    })
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
        .set_feature_end(Position::try_from(end1 - 1)?) // Method is called *after* the end of region is reached, so we need to subtract 1.
        .build();
    writer.write_feature_record(&record)?;
    Ok(())
}

fn update_mask_region<X: std::io::Write>(
    region_start: &mut Option<usize>,
    masking: bool,
    writer_opt: Option<&mut BedWriter<3, X>>,
    record_name: &str,
    index1: usize,
) {
    if let Some(writer) = writer_opt {
        if masking {
            if region_start.is_none() {
                *region_start = Some(index1);
            }
        } else if let Some(start1) = *region_start {
            let _ = write_bed_record(writer, record_name, start1, index1);
            *region_start = None;
        }
    }
}

fn ensure_full_match_regex(regex: &str) -> String {
    let start_ok = regex.starts_with('^');
    let end_ok = regex.ends_with('$');
    if start_ok && end_ok {
        regex.to_string()
    } else if start_ok && !end_ok {
        format!("{}$", regex)
    } else if !start_ok && end_ok {
        format!("^{}", regex)
    } else {
        format!("^{}$", regex)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn ensure_full_match_regex_ok() {
        assert_eq!(ensure_full_match_regex(".*"), "^.*$");
        assert_eq!(ensure_full_match_regex("^.*"), "^.*$");
        assert_eq!(ensure_full_match_regex(".*$"), "^.*$");
        assert_eq!(ensure_full_match_regex("^.*$"), "^.*$");
        assert_eq!(ensure_full_match_regex("abc"), "^abc$");
        assert_eq!(ensure_full_match_regex("^abc"), "^abc$");
        assert_eq!(ensure_full_match_regex("abc$"), "^abc$");
        assert_eq!(ensure_full_match_regex("^abc$"), "^abc$");
    }

    #[test]
    fn process_fasta_record_ok() -> Result<(), Box<dyn Error>> {
        let record = FastaRecord::new(
            noodles_fasta::record::Definition::new("test_sequence", None),
            noodles_fasta::record::Sequence::from(
                b"CTGTGCTGGCATAGTGGTCTCACCTCCGGCAGtatcaccaccactgggcacaagcttctccagcacagcaNNnNnactgtgtcttatttctccttgtactcccagtgttcacaccatgctgcactcacagaagactcttcgttgatattt".to_vec()),
        );

        let tmpdir = tempfile::tempdir()?;
        let stats = process_fasta_record(&record, Some(&tmpdir.path()), ".*");
        assert!(stats.is_some());
        let stats = stats.unwrap();

        // Check the output BED files
        let non_masked_bed_path = tmpdir.path().join("test_sequence.non-masked.bed");
        let non_masked_bed = std::fs::read_to_string(non_masked_bed_path)?;
        assert_eq!(non_masked_bed, "test_sequence\t0\t32\n");
        
        let soft_masked_bed_path = tmpdir.path().join("test_sequence.soft-masked.bed");
        let soft_masked_bed = std::fs::read_to_string(soft_masked_bed_path)?;
        assert_eq!(soft_masked_bed, "test_sequence\t32\t70\ntest_sequence\t75\t150\n");

        let hard_masked_bed_path = tmpdir.path().join("test_sequence.hard-masked.bed");
        let hard_masked_bed = std::fs::read_to_string(hard_masked_bed_path)?;
        assert_eq!(hard_masked_bed, "test_sequence\t70\t75\n");

        assert_eq!(stats.sequence_name, "test_sequence");
        assert_eq!(stats.non_masked_bases, 32);
        assert_eq!(stats.non_masked_ratio, 32.0 / 150.0);

        assert_eq!(stats.soft_masked_bases, 113);
        assert_eq!(stats.soft_masked_ratio, 113.0 / 150.0);

        assert_eq!(stats.hard_masked_bases, 5);
        assert_eq!(stats.hard_masked_ratio, 5.0 / 150.0);
        
        assert_eq!(stats.gc_content, 73.0 / 150.0);
        assert_eq!(stats.sequence_length, 150);
        assert_eq!(stats.checksum_sha256, "39d1aba0a51cb46ce7cef81ac808bf46e3594b0ccabcd77fabf19c2a395651fa");
        Ok(())
    }
    
}
