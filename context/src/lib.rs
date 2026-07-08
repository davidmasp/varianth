use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles::core;
use noodles::fasta;
use noodles::fasta::io::indexed_reader::Builder;

use noodles::vcf;
use noodles::vcf::header::record::value::map::info::Number;
use noodles::vcf::header::record::value::map::Info;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::variant::io::Write;
use noodles::vcf::variant::record_buf::info::field::Value;

// regions and positions are 1-based (!!!)
// for how to write custom fields in header and in the record see
// https://github.com/zaeleus/noodles/issues/160#issuecomment-1509508247

fn write_nnn_string(k: usize) -> String {
    let size = (2 * k) + 1;
    std::iter::repeat("N").take(size).collect()
}

fn get_ntp_from_record(
    vcf_record: &vcf::variant::RecordBuf,
    fasta_index_reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,
    k: usize,
) -> Result<String> {
    let pos1 = vcf_record
        .variant_start()
        .context("VCF record is missing a variant start position")?;
    let end = pos1
        .checked_add(k)
        .with_context(|| format!("variant position {pos1} with context size {k} overflows"))?;

    // If the requested left flank falls before the start of the reference, keep the
    // existing fallback behavior and emit an all-N context string.
    let Some(start) = usize::from(pos1)
        .checked_sub(k)
        .and_then(|position| core::Position::try_from(position).ok())
    else {
        return Ok(write_nnn_string(k));
    };

    let chrom = vcf_record.reference_sequence_name().to_string();
    let tntp_region = core::Region::new(chrom, start..=end);

    let tntp = match fasta_index_reader.query(&tntp_region) {
        Ok(v) => v,
        Err(_e) => return Ok(write_nnn_string(k)),
    };

    std::str::from_utf8(tntp.sequence().as_ref())
        .context("FASTA query returned non-UTF-8 sequence data")
        .map(str::to_owned)
}

pub fn addms(
    genome: PathBuf,
    variants_in: PathBuf,
    variants_out: PathBuf,
    kval: usize,
    key_name: String,
    key_description: String,
) -> Result<()> {
    let reference_path: PathBuf = genome;
    let vcf_path: PathBuf = variants_in;
    let vcf_path_out: PathBuf = variants_out;

    let mut reference_reader = Builder::default()
        .build_from_path(&reference_path)
        .with_context(|| {
            format!(
                "failed to open indexed FASTA reference {}",
                reference_path.display()
            )
        })?;

    /* here we need to decide if stdin is used, not sure how to do that yet */

    let mut variants_reader = vcf::io::reader::Builder::default()
        .build_from_path(&vcf_path)
        .with_context(|| format!("failed to open input VCF {}", vcf_path.display()))?;

    let header = variants_reader
        .read_header()
        .with_context(|| format!("failed to read VCF header from {}", vcf_path.display()))?;

    let mut writer = vcf::io::writer::Builder::default()
        .build_from_path(&vcf_path_out)
        .with_context(|| format!("failed to create output VCF {}", vcf_path_out.display()))?;

    let mut header_out = header.clone();
    // Parse non-standard keys using `info::Key::from_str`.
    let ms_key = key_name;
    // Create structured header records using `Map<I>`.
    let ms_value = Map::<Info>::new(
        Number::Count(1),
        noodles::vcf::header::record::value::map::info::Type::String,
        key_description,
    );
    header_out.infos_mut().insert(ms_key.clone(), ms_value);
    writer
        .write_header(&header_out)
        .with_context(|| format!("failed to write VCF header to {}", vcf_path_out.display()))?;

    for (record_index, result) in variants_reader.record_bufs(&header).enumerate() {
        let mut record_out = result.with_context(|| {
            format!(
                "failed to read VCF record {} from {}",
                record_index + 1,
                vcf_path.display()
            )
        })?;
        let tntp_results = get_ntp_from_record(&record_out, &mut reference_reader, kval)
            .with_context(|| {
                format!(
                    "failed to calculate {ms_key} for VCF record {}",
                    record_index + 1
                )
            })?;
        record_out
            .info_mut()
            .insert(ms_key.clone(), Some(Value::String(tntp_results)));

        writer
            .write_variant_record(&header_out, &record_out)
            .with_context(|| {
                format!(
                    "failed to write VCF record {} to {}",
                    record_index + 1,
                    vcf_path_out.display()
                )
            })?;
    }

    Ok(())
}
