

mod rmpio;
use rmpio::bamread::reader_from_path;

mod column;
use column::iterators::MpileupIterator;

use noodles::sam::alignment::record::MappingQuality;

use std::io::{self, BufWriter};
use std::io::Write;
use std::str::from_utf8;

use noodles::bam::io::Reader as BamReader;
use noodles::bgzf::Reader as GzReader;
use std::fs::File;
use noodles::sam::Header;

use noodles::fasta::indexed_reader::Builder;
use noodles::fasta::indexed_reader::IndexedReader;
use noodles::fasta::io::BufReader;
use noodles::core::{Position, Region};
use noodles::fasta::record::Sequence;

fn main() {

    let mut fasta_reader: IndexedReader<BufReader<File>> = Builder::default().build_from_path("./hg38.fa").unwrap();

    let stdout = io::stdout();
    let mut handle = BufWriter::new(stdout.lock());

    // PARAMS
    let exclude_flags: u16 = 3584;

    // I am not entirely sure if this is the 
    // default in samtools, because in the docs says NULL.
    let include_flags_option: Option<u16> = None;

    // So FOR MAPPING QUALITY we should create a value outside the loop and
    // then compare it inside
    let min_mapq: MappingQuality = MappingQuality::new(24).unwrap();

    // let path = "minsample.bam";
    let path = "minsample.bam";

    let (reader, bam_header): (BamReader<GzReader<File>>, Header) = reader_from_path(path).expect("msg");

    let iterator = match MpileupIterator::new(reader, min_mapq, exclude_flags, include_flags_option) {
        Ok(iterator) => iterator,
        Err(e) => panic!("Failed to create iterator: {}", e),
    };

    for mut i in iterator {
        let _ = i.add_reference_base(&bam_header, &mut fasta_reader);
        let row = i.format_column(&bam_header);
        writeln!(handle, "{}", row).unwrap();
    }

}
