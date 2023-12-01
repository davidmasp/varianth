

use crate::core::bed_record_to_region;

use std::ops::Bound;
use std::path::PathBuf;
use std::fs::File;
use std::io::BufReader;

use std::collections::HashMap;

use noodles::core;
use noodles::sam;
use noodles::bam;
use noodles::bed;

use std::io::Write;
use std::fs::OpenOptions;


pub fn readfreq(reads: PathBuf, sites_in: PathBuf, tsv_out: PathBuf) -> () {

    // a block to define the output file
    let mut out_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(tsv_out)
        .unwrap();

    let sites_in_file_result = File::open(sites_in);
    let sites_in_file = match sites_in_file_result {
        Ok(file) => BufReader::new(file),
        Err(error) => panic!("Problem opening the file: {:?}", error),
    };

    let mut variants_reader = bed::Reader::new(sites_in_file);

    // a block to define the bam reader
    let bam_path = reads;
    let mut bam_reader: bam::IndexedReader<noodles::bgzf::Reader<std::fs::File>> = bam::indexed_reader::Builder::default()
        .build_from_path(bam_path).unwrap();
    let bam_header: sam::Header = bam_reader.read_header().unwrap();

    for bed_record_result in variants_reader.records::<3>() {
        let bed_record = match bed_record_result {
            Ok(bed_record) => bed_record,
            Err(error) => panic!("Problem reading the file: {:?}", error),
        };

        let region = bed_record_to_region(bed_record);

        let hash_freq = get_readfrequency(
            &region,
            &mut bam_reader,
            &bam_header,
        );

        let region_seqname = region.name();
        let region_start_bpos = region.start();
        let region_start_pos = get_position(region_start_bpos).unwrap();
        let region_start_usize = usize::from(region_start_pos);
        let region_start_usize_0base = region_start_usize.checked_sub(1).unwrap();
        let region_start_string = region_start_usize_0base.to_string();
        
        let region_end_bpos = region.end();
        let region_end_pos = get_position(region_end_bpos).unwrap();
        let region_end_string = usize::from(region_end_pos).to_string();

        for (key, value) in &hash_freq {
            let output = format!("{}\t{}\t{}\t{}\t{}",
                region_seqname,
                region_start_string,
                region_end_string,
                key,
                value);
            //println!("{}", output);
            writeln!(out_file, "{}", output).unwrap();
        }

        //println!("{:?}", hash_freq);

    } 

    /*
    let mut records = variants_reader.records::<3>()
        .map(|bed_record_result| {
            

        });
    */


}

fn get_position(bound: Bound<core::Position>) -> Option<core::Position>  {
    let position_result: Option<core::Position> = match bound {
        Bound::Included(position) => Some(position),
        Bound::Excluded(position) => Some(position),
        Bound::Unbounded => None,
    };
    position_result
}

fn get_readfrequency(
    region_in: &core::Region,
    bam_reader: &mut bam::IndexedReader<noodles::bgzf::Reader<std::fs::File>>,
    bam_header: &sam::Header,
) -> HashMap<String, usize> {

    let region_start_usize = usize::from(get_position(region_in.start()).unwrap());

    let region_end_usize = usize::from(get_position(region_in.end()).unwrap());
    let region_end_usize_plus1_result = region_end_usize.checked_add(1); // this is for 0-based range
    let region_end_usize_plus1 = match region_end_usize_plus1_result {
        Some(region_end_usize_plus1) => region_end_usize_plus1,
        None => panic!("Problem with usize subtraction, negative position?"),
    };


    let query = bam_reader.query(bam_header, &region_in).unwrap();

    let mut hash_read_counts: HashMap<String, usize> = HashMap::new();

    for result in query {
        //println!("{} {}", alig_start_usize, region_start_usize);

        let record = result.unwrap();

        let record_flags = record.flags();

        if record_flags.is_qc_fail() || record_flags.is_duplicate() || record_flags.is_secondary() || record_flags.is_supplementary() {
            continue;
        }

        let alig_start = record.alignment_start().unwrap();
        let alig_start_usize = usize::from(alig_start);
        let alig_start_usize_minus1_result = alig_start_usize.checked_sub(1); // this is for 1-based substraction later
        
        let alig_start_usize_minus1 = match alig_start_usize_minus1_result {
            Some(alig_start_usize_minus1) => alig_start_usize_minus1,
            None => panic!("Problem with usize subtraction, negative position?"),
        };

        if alig_start_usize > region_start_usize {
            // this means that the aligment starts after the region, thus we can not get the full sequence and
            // we skip this read
            continue;
        }

        let alig_sequence = record.sequence();
        
        let mut tmp_seq_str = "".to_string();
        for i in region_start_usize..region_end_usize_plus1 {

            let sequence_idx_result = i.checked_sub(alig_start_usize_minus1);

            let sequence_idx = match sequence_idx_result {
                Some(sequence_idx) => sequence_idx,
                None => panic!("Problem with usize subtraction, negative position?"),
            };
            

            let sequence_idx_position = core::Position::try_from(sequence_idx).unwrap();

            let seq_val_result = alig_sequence.get(sequence_idx_position);
            let seq_val = match seq_val_result {
                Some(seq_val) => seq_val,
                None => panic!("Problem with sequence index"),
            };

            let seq_val_char = match seq_val {
                sam::record::sequence::Base::A => 'A',
                sam::record::sequence::Base::C => 'C',
                sam::record::sequence::Base::G => 'G',
                sam::record::sequence::Base::T => 'T',
                sam::record::sequence::Base::N => 'N',
                _ => 'O'
            };

            tmp_seq_str.push(seq_val_char);

        }
        //println!("{}", tmp_seq_str);
        let stat = hash_read_counts.entry(tmp_seq_str).or_insert(0);
        *stat += 1;
    }

    hash_read_counts


}

