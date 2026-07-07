use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;

use noodles::vcf::variant::record::info::field::Value::Array as VcfArrayValue;
use noodles::vcf::variant::record::info::field::value::Array as VcfArray;

use noodles::vcf::Record as VcfRecord;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn process_record(record: VcfRecord, header: &vcf::Header) -> Vec<String> {
    let chr_name = record.reference_sequence_name();
    let pos = record.variant_start().unwrap().unwrap();
    let reference = record.reference_bases();
    let alt_str = record.alternate_bases();

    if alt_str.len() != 1 {
        eprintln!(
            "More than one alternate base found for record at {}:{}",
            chr_name, pos
        );
        std::process::exit(1);
    }

    let alt = alt_str
        .iter()
        .next()
        .unwrap_or_else(|| {
            eprintln!("No alternate base found");
            std::process::exit(1);
        })
        .unwrap_or_else(|e| {
            eprintln!("Error parsing alternate base: {}", e);
            std::process::exit(1);
        });

    let position_str = format!("{}|{}|{}|{}", chr_name, pos, reference, alt);

    //println!("{:?}", record);
    let info_record = record.info();

    let csq_field_option = info_record.get(&header, "CSQ");

    match csq_field_option {
        Some(field_value_raw) => {
            let value = field_value_raw.unwrap().unwrap();
            //println!("CSQ: {:?}", value);
            if let VcfArrayValue(arr) = value {
                match arr {
                    VcfArray::String(values) => {
                        let mut results = Vec::new();
                        for item in values.iter() {
                            if let Ok(Some(s)) = item {
                                let result_string = format!("{}|{}", position_str, s);
                                results.push(result_string);
                            }
                        }
                        return results;
                    }
                    _ => {
                        eprintln!("CSQ array is not of String type");
                        std::process::exit(1);
                    }
                }
            } else {
                eprintln!("CSQ is not array type");
                std::process::exit(1);
            }
        }
        None => {
            eprintln!("CSQ not found");
            std::process::exit(1);
        }
    }
}

pub fn vep2table(vcf_path: impl AsRef<Path>, output_path: impl AsRef<Path>) {
    // let vcf_path = "../path_to_mutation.annotated.vcf.gz";
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(vcf_path)
        .unwrap();
    let header = reader.read_header().unwrap();
    let infos = header.infos();
    let csq_info = infos.get("CSQ").unwrap();
    let csq_description = csq_info.description();

    // Create output file writer
    let file = File::create(output_path).unwrap();
    let mut writer = BufWriter::new(file);

    // Extract format fields after "Format: "
    if let Some(format_pos) = csq_description.find("Format: ") {
        let format_str = &csq_description[format_pos + 8..];
        let header_str = format!("chrom|pos|ref|alt|{}", format_str);
        writeln!(writer, "{}", header_str).unwrap();
    } else {
        eprintln!("Format field not found in CSQ description");
        std::process::exit(1);
    }

    for result in reader.records() {
        let record = result.unwrap();
        let output_lines = process_record(record, &header);
        for line in output_lines {
            writeln!(writer, "{}", line).unwrap();
        }
    }
}
