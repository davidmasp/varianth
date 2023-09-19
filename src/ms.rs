

use std::path::PathBuf;

use noodles::core;
use noodles::fasta;
use noodles::fasta::indexed_reader::Builder;

use noodles::vcf;
use noodles::vcf::header::record::value::map::Info;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::header::Number;

// regions and positions are 1-based (!!!)
// for how to write custom fields in header and in the record see
// https://github.com/zaeleus/noodles/issues/160#issuecomment-1509508247

fn get_ntp_from_record(
    vcf_record: vcf::Record,
    fasta_index_reader: &mut fasta::IndexedReader<
        Box<dyn noodles::fasta::io::BufReadSeek>,
    >,
    k: usize,
) -> String {
    let pos1 = core::Position::try_from(usize::from(
        vcf_record.position(),
    ))
    .unwrap();
    let end = pos1.checked_add(k).unwrap();
    // for some reason there is no substr method in the noodles?
    let start = core::Position::try_from(
        usize::from(pos1).checked_sub(k).unwrap(),
    )
    .unwrap();

    let chrom = vcf_record.chromosome().to_string();
    let tntp_region = core::Region::new(chrom, start..=end);

    let tntp =
        fasta_index_reader.query(&tntp_region).unwrap();

    let out_str = String::try_from(
        std::str::from_utf8(tntp.sequence().as_ref())
            .unwrap(),
    )
    .unwrap();
    out_str
}

pub fn addms(genome: PathBuf, variants_in: PathBuf, variants_out: PathBuf, kval: usize, key_name: String, key_description: String, use_stdin: bool, use_stdout: bool) {
    /*
    let reference_path: PathBuf = Into::into("reference.fa");
    let vcf_path: PathBuf = Into::into("sample.vcf.gz");
    let vcf_path_out: PathBuf = Into::into("out.vcf.gz");
     */
    let reference_path: PathBuf = genome;
    let vcf_path: PathBuf = variants_in;
    let vcf_path_out: PathBuf = variants_out;

    let mut reference_reader = Builder::default()
        .build_from_path(reference_path)
        .unwrap();

    /* here we need to decide if stdin is used, not sure how to do that yet */

    let mut variants_reader = vcf::reader::Builder::default()
            .build_from_path(vcf_path)
            .unwrap();

    let header = variants_reader.read_header().unwrap();

    let mut writer = vcf::writer::Builder
            .build_from_path(vcf_path_out)
            .unwrap();

    let mut header_out = header.clone();
    // Parse non-standard keys using `info::Key::from_str`.
    let ms_key: vcf::record::info::field::Key = key_name.parse().unwrap();
    // Create structured header records using `Map<I>`.
    let ms_value = Map::<Info>::new(
        Number::Count(1),
        noodles::vcf::header::record::value::map::info::Type::String,
        key_description,
    );
    header_out.infos_mut().insert(ms_key.clone(), ms_value);
    writer.write_header(&header_out).unwrap();

    // i think we should map this
    for result in variants_reader.records(&header) {
        let record = result.unwrap();
        let mut record_out = record.clone();
        let tntp_results = get_ntp_from_record(
            record,
            &mut reference_reader,
            kval,
        );
        record_out.info_mut().insert(
            ms_key.clone(),
            Some(vcf::record::info::field::Value::String(
                tntp_results.clone(),
            )),
        );

        writer.write_record(&header, &record_out).unwrap();

    }
}

