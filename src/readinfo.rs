
use crate::core::fromvcfrecord2region;

use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use serde_json;

use std::collections::HashMap;

use noodles::vcf;
use noodles::sam;
use noodles::bam;

pub fn readinfo(reads: PathBuf, variants_in: PathBuf, json_out: PathBuf) {

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

    let mut hash_readpos: HashMap<usize, u32> = HashMap::new(); 

    for result in variants_reader.records(&variants_header){
        let variant = result.unwrap();
        get_readinfo_from_record(
            &variant,
            &mut bam_reader,
            &bam_header,
            &mut hash_readpos,
        );
    }

    let json = serde_json::to_string(&hash_readpos).unwrap();
    let mut file = File::create(json_out).unwrap();
    file.write_all(json.as_bytes()).unwrap();

}


fn get_readinfo_from_record(
    variant: &vcf::Record,
    bam_reader: &mut bam::IndexedReader<noodles::bgzf::Reader<std::fs::File>>,
    bam_header: &sam::Header,
    hash_ref: &mut HashMap<usize, u32>,
) -> () {
    let variant_position = variant.position();
    let variant_position_usize = usize::from(variant_position);
    let region = fromvcfrecord2region(variant);
    let query = bam_reader.query(bam_header, &region).unwrap();

    for result in query {
        let record = result.unwrap();
        let alig_start = record.alignment_start().unwrap();
        // let alig_end  = record.alignment_end().unwrap();
        let pos_start = usize::from(alig_start);
        let size_to_start = variant_position_usize - pos_start;

        let stat = hash_ref.entry(size_to_start).or_insert(0);
        *stat += 1;
    }
}

