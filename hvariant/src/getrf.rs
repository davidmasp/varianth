

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
use noodles::sam::record::cigar::Cigar;
use noodles::sam::record::cigar::op::kind::Kind;

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
    for alig_result in query {
        let alig = alig_result.unwrap();
        let record_flags = alig.flags();

        if record_flags.is_qc_fail() || record_flags.is_duplicate() || record_flags.is_secondary() || record_flags.is_supplementary() {
            continue;
        }

        let alig_start = alig.alignment_start().unwrap();
        let alig_start_usize = usize::from(alig_start);
        let alig_end =  alig.alignment_end().unwrap();
        let alig_end_usize = usize::from(alig_end);

        let alig_sequence = alig.sequence().to_string();
        let alig_cigar = alig.cigar();
        println!("{:?}", alig_cigar);

        let sel_seq_result = get_seq(
            &alig_sequence,
            &alig_start_usize,
            &alig_end_usize,
            &region_start_usize,
            &region_end_usize_plus1,
            alig_cigar,
        );

        let sel_seq = match sel_seq_result {
            Some(sel_seq) => sel_seq,
            None => continue,
        };
        println!("seq: {}, pos:{} ", sel_seq, region_end_usize);
        let stat = hash_read_counts.entry(sel_seq).or_insert(0);
        *stat += 1;
    }
    hash_read_counts
}


/*

pub enum Kind {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SequenceMatch,
    SequenceMismatch,
}

From what I can gather,

- M or Match/Mismatch in the read should do nothing.
- SoftClip in the read should do nothing. ?? VERIFY!
- SequenceMatch & SequenceMismatch (these are the new M, = & X)

- Deletion (in the read) should decrease the positions of the region.
- Skip regions (in the read) should decrease the positions of the region.
        note that these aren't going to be common in WGS data but yes in RNA
- HardClip (in the read) should decrease the positions of the region. ?? VERIFY!

- Insertion (in the read) should increase the positions of the region.
*/


fn get_seq(sequence: &str, astart: &usize, aend: &usize, rstart: &usize, rend: &usize, alig_cigar: &Cigar) -> Option<String> {
    if rstart < astart {
        return None;
    }

    if rend > aend {
        return None;
    }

    let mut rel_start = rstart - astart + 1;
    let mut rel_end = rend - astart + 1;

    /*
    let mut mod_rel_start = rel_start;
    let mut mod_rel_end = rel_end;
     */
    let mut del_pos_modifier: usize = 0;
    let mut softclip_pos_modifier: usize = 0;
    let mut cursor : usize = 1;
    let mut new_sequence = "".to_string();
    for i in 0..alig_cigar.len() {
        let cigar_op = alig_cigar[i]; // operation
        let chunk_size: usize = cigar_op.len(); // block size
        let chunk_type = cigar_op.kind();
        match chunk_type {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let idx_seq = cursor - 1 - del_pos_modifier + softclip_pos_modifier;
                let idx_seq_end = cursor + chunk_size - 1 - del_pos_modifier + softclip_pos_modifier;
                let extra_seq = sequence[idx_seq..idx_seq_end].to_string();
                new_sequence.push_str(&extra_seq);
                cursor += chunk_size;
            },
            Kind::SoftClip => {
                softclip_pos_modifier += chunk_size;
            },
            Kind::HardClip | Kind::Pad  => {
                return None;
            },
            Kind::Insertion => {
                let idx_seq = cursor - 1 - del_pos_modifier + softclip_pos_modifier;
                let idx_seq_end = cursor + chunk_size - 1 - del_pos_modifier + softclip_pos_modifier;
                let extra_seq = sequence[idx_seq..idx_seq_end].to_string();
                new_sequence.push_str(&extra_seq);
                // cursor does not move, but relative positions do
                if rel_start > cursor {
                    rel_start += chunk_size;
                }
                if rel_end > cursor {
                    rel_end += chunk_size;
                }
            },
            Kind::Deletion => {
                // here adds "-"
                let deletion_seq = "-".repeat(chunk_size);
                new_sequence.push_str(&deletion_seq);
                // cursor advances
                cursor += chunk_size;
                del_pos_modifier += chunk_size;
            },
            Kind::Skip => {
                // here adds "N"
                let skip_seq = "N".repeat(chunk_size);
                new_sequence.push_str(&skip_seq);
                // cursor advances
                cursor += chunk_size;
                del_pos_modifier += chunk_size;
            },
        };
    }
    // this is only to convert 1-based input to 0-based substring
    let rel_start_idx: usize = rel_start - 1;
    let rel_end_idx: usize = rel_end - 1;
    println!("{}",new_sequence);
    let result_seq = &new_sequence[rel_start_idx..rel_end_idx];
    Some(result_seq.to_string())
}

