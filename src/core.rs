

use noodles::core;
use noodles::vcf;
use noodles::bed;

pub fn fromvcfrecord2region(vcf_record: &vcf::Record) -> core::Region {
    let pos1 = core::Position::try_from(usize::from(
        vcf_record.position(),
    ))
    .unwrap();

    let chrom = vcf_record.chromosome().to_string();
    let region = core::Region::new(chrom, pos1..=pos1);

    region
}

pub fn bed_record_to_region(record: bed::Record<3>) -> core::Region {
    let chr_name = record.reference_sequence_name();
    let pos1 = record.start_position();
    let pos2 = record.end_position();
    let region = core::Region::new(chr_name, pos1..=pos2);
    region
}

/* 
pub fn build_bed_record3(chr_name: &str, pos1: core::Position, pos2: core::Position) -> bed::Record<3> {

    let record_result = bed::Record::<3>::builder()
        .set_reference_sequence_name(chr_name)
        .set_start_position(pos1)
        .set_end_position(pos2)
        .build();
        
    let record = match record_result {
        Ok(record) => record,
        Err(error) => panic!("Problem building bed record: {:?}", error),
    };
    record
}


*/
