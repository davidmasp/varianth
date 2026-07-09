use anyhow::{Context, Result, bail};
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;

use noodles::vcf::variant::record::info::field::Value::Array as VcfArrayValue;
use noodles::vcf::variant::record::info::field::value::Array as VcfArray;

use noodles::vcf::Record as VcfRecord;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn process_record(record: VcfRecord, header: &vcf::Header) -> Result<Vec<String>> {
    let chr_name = record.reference_sequence_name();
    let pos = record
        .variant_start()
        .with_context(|| format!("VCF record on {chr_name} is missing a variant start"))?
        .with_context(|| format!("failed to parse variant start for record on {chr_name}"))?;
    let reference = record.reference_bases();
    let alt_str = record.alternate_bases();

    if alt_str.len() != 1 {
        bail!(
            "expected exactly one alternate base for record at {}:{}, found {}",
            chr_name,
            pos,
            alt_str.len()
        );
    }

    let alt = alt_str
        .iter()
        .next()
        .with_context(|| format!("no alternate base found for record at {chr_name}:{pos}"))?
        .with_context(|| {
            format!("failed to parse alternate base for record at {chr_name}:{pos}")
        })?;

    let position_str = format!("{}|{}|{}|{}", chr_name, pos, reference, alt);

    //println!("{:?}", record);
    let info_record = record.info();

    let csq_field_option = info_record.get(&header, "CSQ");

    let field_value_raw = csq_field_option
        .with_context(|| format!("CSQ INFO field not found for record {position_str}"))?;
    let value = field_value_raw
        .with_context(|| format!("failed to parse CSQ INFO field for record {position_str}"))?
        .with_context(|| format!("CSQ INFO field has no value for record {position_str}"))?;

    match value {
        VcfArrayValue(VcfArray::String(values)) => {
            let mut results = Vec::new();
            for item in values.iter() {
                if let Some(s) = item.with_context(|| {
                    format!("failed to parse CSQ entry for record {position_str}")
                })? {
                    let result_string = format!("{}|{}", position_str, s);
                    results.push(result_string);
                }
            }
            Ok(results)
        }
        VcfArrayValue(_) => bail!("CSQ array is not of String type for record {position_str}"),
        _ => bail!("CSQ is not array type for record {position_str}"),
    }
}

pub fn vep2table(vcf_path: impl AsRef<Path>, output_path: impl AsRef<Path>) -> Result<()> {
    // let vcf_path = "../path_to_mutation.annotated.vcf.gz";
    let vcf_path = vcf_path.as_ref();
    let output_path = output_path.as_ref();

    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(vcf_path)
        .with_context(|| format!("failed to open input VCF {}", vcf_path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read VCF header from {}", vcf_path.display()))?;
    let infos = header.infos();
    let csq_info = infos
        .get("CSQ")
        .context("VCF header is missing CSQ INFO definition")?;
    let csq_description = csq_info.description();

    // Create output file writer
    let file = File::create(output_path)
        .with_context(|| format!("failed to create output file {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // Extract format fields after "Format: "
    if let Some(format_pos) = csq_description.find("Format: ") {
        let format_str = &csq_description[format_pos + 8..];
        let header_str = format!("chrom|pos|ref|alt|{}", format_str);
        writeln!(writer, "{}", header_str).with_context(|| {
            format!("failed to write output header to {}", output_path.display())
        })?;
    } else {
        bail!("Format field not found in CSQ description");
    }

    for (record_index, result) in reader.records().enumerate() {
        let record = result.with_context(|| {
            format!(
                "failed to read VCF record {} from {}",
                record_index + 1,
                vcf_path.display()
            )
        })?;
        let output_lines = process_record(record, &header)
            .with_context(|| format!("failed to process VCF record {}", record_index + 1))?;
        for line in output_lines {
            writeln!(writer, "{}", line).with_context(|| {
                format!(
                    "failed to write output line for VCF record {} to {}",
                    record_index + 1,
                    output_path.display()
                )
            })?;
        }
    }

    writer
        .flush()
        .with_context(|| format!("failed to flush output file {}", output_path.display()))?;

    Ok(())
}
