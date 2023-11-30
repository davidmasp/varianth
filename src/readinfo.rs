
use crate::core::fromvcfrecord2region;

use std::path::PathBuf;

use noodles::vcf;
use noodles::vcf::header::record::value::map::Info;
use noodles::vcf::header::record::value::Map;

use noodles::sam;
use noodles::bam;

use noodles::core;

// the idea is to go from a vcf file to a vcf file with a new info field
// that contains read information

// another thing that I want to do related to this might be to get the
// vaf from a random position in many bams. Similar but unrelated.

pub fn readinfo(reads: PathBuf, variants_in: PathBuf, variants_out: PathBuf, key_name: String, key_description: String) {

    // block to define the variant reader
    let vcf_path: PathBuf = variants_in;
    let mut variants_reader = vcf::reader::Builder::default()
    .build_from_path(vcf_path)
    .unwrap();
    let variants_header: vcf::Header = variants_reader.read_header().unwrap();

    // a block to define the bam reader
    let bam_path = reads;
    let mut bam_reader: bam::IndexedReader<noodles::bgzf::Reader<std::fs::File>> = bam::indexed_reader::Builder::default()
        .build_from_path(bam_path).unwrap();
    let bam_header: sam::Header = bam_reader.read_header().unwrap();

    // a block for the variants writer
    let vcf_path_out: PathBuf = variants_out;
    let mut writer = vcf::writer::Builder
        .build_from_path(vcf_path_out)
        .unwrap();

    let mut variants_header_out: vcf::Header = variants_header.clone();
    // Parse non-standard keys using `info::Key::from_str`.
    let rdinfo_key: vcf::record::info::field::Key = key_name.parse().unwrap();
    // Create structured header records using `Map<I>`.
    let rdinfo_value = Map::<Info>::new(
        vcf::header::Number::Count(1),
        noodles::vcf::header::record::value::map::info::Type::String,
        key_description,
    );

    variants_header_out.infos_mut().insert(rdinfo_key.clone(), rdinfo_value);
    writer.write_header(&variants_header_out).unwrap();

    for result in variants_reader.records(&variants_header){
        let variant = result.unwrap();
            let readinfo_result = get_readinfo_from_record(
                &variant,
                &mut bam_reader,
                &bam_header
            );
        println!("{}", readinfo_result);

    }


        /*
        // write block
        let mut record_out = record.clone();
        record_out.info_mut().insert(
            rdinfo_key.clone(),
            Some(vcf::record::info::field::Value::String(
                readinfo_result.clone(),
            )),
        );

        writer.write_record(&header, &record_out).unwrap();
         */

}


fn get_readinfo_from_record(
    variant: &vcf::Record,
    bam_reader: &mut bam::IndexedReader<noodles::bgzf::Reader<std::fs::File>>,
    bam_header: &sam::Header,
) -> String {
    let variant_position = variant.position();
    let variant_position_usize = usize::from(variant_position);
    let region = fromvcfrecord2region(variant);
    let query = bam_reader.query(bam_header, &region).unwrap();

    let mut out_str = String::new();

    for result in query {
        let record = result.unwrap();
        let alig_start = record.alignment_start().unwrap();
        // let alig_end  = record.alignment_end().unwrap();
        let pos_start = usize::from(alig_start);
        let size_to_start = variant_position_usize - pos_start;
        let size_to_start_str = size_to_start.to_string();
        out_str.push_str(&size_to_start_str);
        out_str.push_str(":");
    }

    out_str
}

