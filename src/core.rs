

use noodles::core;
use noodles::vcf;

pub fn fromvcfrecord2region(vcf_record: &vcf::Record) -> core::Region {
    let pos1 = core::Position::try_from(usize::from(
        vcf_record.position(),
    ))
    .unwrap();

    let chrom = vcf_record.chromosome().to_string();
    let region = core::Region::new(chrom, pos1..=pos1);

    region
}


