use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;

use crate::column::column::build_mpileupcolumn;
use crate::column::column::MpileupColumn;
use noodles::bam::record::Record;
use noodles::core::Position;
use std::collections::BTreeMap;


fn parse_cigar(
    cigar_operation: Op,
    record_object: &Record,
    btmap: &mut BTreeMap<Position, MpileupColumn>, // this is global
    ref_pos: &mut Position,
    read_pos: &mut usize,
) -> () {
    let kind_op = cigar_operation.kind();

    match kind_op {
        Kind::Match
        | Kind::SequenceMatch
        | Kind::SequenceMismatch
        | Kind::Skip
        | Kind::Deletion => {
            // consumes query and reference
            let mut op_len = cigar_operation.len();
            while op_len > 0 {
                if btmap.contains_key(ref_pos) {
                    let entry = btmap.get_mut(ref_pos).unwrap();
                    entry.add_record_to_column(&record_object, &kind_op, read_pos, None);
                } else {
                    let def_mpc = build_mpileupcolumn(
                        &record_object.reference_sequence_id().unwrap().unwrap(),
                        &ref_pos,
                        &record_object,
                        &kind_op,
                        read_pos,
                        None,
                    );
                    btmap.insert(*ref_pos, def_mpc);
                }
                *ref_pos = ref_pos.checked_add(1).unwrap();
                if (kind_op != Kind::Skip) & (kind_op != Kind::Deletion) {
                    // skip does not consume the "read"
                    // skips are introns, this might be a problem
                    *read_pos = read_pos.checked_add(1).unwrap();
                }
                op_len -= 1;
            }
        }
        Kind::HardClip | Kind::Pad => {
            // consumes nothing, also not included in the Pileup
        }
        Kind::SoftClip => {
            // consumes query but not reference, is not included in the Pileup
            *read_pos = read_pos.checked_add(cigar_operation.len()).unwrap();
        }
        Kind::Insertion => {
            // consumes query but not reference
            // the problem is that because the consumption has already happen
            // we need to go back to the previous position
            // and it also needs to get attached to the right reads
            // to find the "same read" we can use the query name
            // BWA might include insertion sites in the first CIGAR operation
            // see https://github.com/annalam/mutato/issues/1
            // now i have some sort of control here: ^jdhfksdj
            let ref_pos_insertion_raw: usize = Position::try_into(ref_pos.clone()).unwrap(); 
            let ref_pos_insertion_usize = ref_pos_insertion_raw.checked_sub(1).unwrap();
            let ref_pos_insertion = Position::try_from(ref_pos_insertion_usize).unwrap();
            let op_len = cigar_operation.len(); // this is equivalent to the insertion length
            if btmap.contains_key(&ref_pos_insertion) {
                let entry = btmap.get_mut(&ref_pos_insertion).unwrap();
                entry.add_record_to_column(&record_object, &kind_op, read_pos, Some(&op_len));
            } else {
                panic!("Insrtion problem, I think this should not happen");
                /*
                let def_mpc = build_mpileupcolumn(
                    &record_object.reference_sequence_id().unwrap().unwrap(),
                    &ref_pos,
                    None,
                    &record_object,
                    &kind_op,
                    read_pos,
                    Some(&op_len),
                );
                btmap.insert(*ref_pos, def_mpc);
                */
            }
            *read_pos = read_pos.checked_add(op_len).unwrap();
        }
    }
    ()
}

pub fn process_read(
    btree_map: &mut BTreeMap<Position, MpileupColumn>,
    record: &Record,
    exclude_flags: u16,
    include_flags_option: Option<u16>,
) -> Result<(), std::io::Error> {

    let flag_bits = record.flags().bits();
    // does the order matter? I don't think so
    let exclude_flags_test = flag_bits & exclude_flags;
    if exclude_flags_test > 0 {
        // if bit-wise AND is not zero, it means that from all the flags in the
        // read, at least one coincides with the exclude test and thus should be removed
        return Ok(());
    }

    match include_flags_option {
        Some(include_flags) => {
            let include_flags_test = flag_bits & include_flags;
            if include_flags_test == 0 {
                // if bit-wise AND is zero, it means that from all the flags in the
                // read, none coincides with the include test and thus should be removed
                // from samtools: rejecting reads that have none of the mask bits set
                return Ok(());
            }
        },
        None => {
            // if there is no include flags, then we should include all reads
        }
    }

    let cig = record.cigar();
    let cig_iter = cig.iter();
    let mut read_pos: usize = 0;
    let ref_pos_opt = record.alignment_start();
    let mut ref_pos: Position = match ref_pos_opt {
        Some(x) => x?,
        None => panic!("No alignment start found"), // fix later
    };

    let mut first_cigar_op = true; 

    for x in cig_iter {
        let cig_op = x?;

        if first_cigar_op & (cig_op.kind() == Kind::Insertion) {
            // ^jdhfksdj
            // this never happens in my test/real data and when I try to fake it
            // and test how samtools behaves it does seem to be ignoring the
            // operation, so should be similar to what I am doing here.
            read_pos = read_pos.checked_add(cig_op.len()).unwrap();
            first_cigar_op = false;
            eprintln!("First CIGAR operation is an insertion, skipping it"); // needs to be a warning!
            continue;
        } else if first_cigar_op {
            first_cigar_op = false;
        } else {
            // do nothing
        }

        parse_cigar(cig_op, 
            record, 
            btree_map, 
            &mut ref_pos, 
            &mut read_pos,);
    }

    Ok(())
}



