
use std::num::{NonZeroUsize, NonZeroU8};
use std::path::PathBuf;
use crate::position::Contig;
use crate::position::Position;
use crate::position::VariantPosition;
use noodles::bed::io::reader::Builder;
use noodles::bed;

pub fn read_positions_from_bed3(filename: &PathBuf) -> Vec<Position> {
    let reader_result = Builder::<3>::default()
                .build_from_path(filename);
    let mut reader = reader_result.unwrap();
    let mut record = bed::Record::default();
    let mut positions = Vec::new();

    while reader.read_record(&mut record).unwrap() != 0 {
        let contig_raw_name = record.reference_sequence_name();
        let contig = Contig::new_from_noodles(contig_raw_name, None);
        // this is already 1-based
        let pos_start = record.feature_start().unwrap().get();
        let pstart: NonZeroUsize = pos_start.try_into().unwrap();
        let pos_end_option = record.feature_end();
        let current_position = match pos_end_option {
            Some(pos_end_result) => {
                let pend = pos_end_result.unwrap().get();
                let width_usize = pend - pos_start + 1;
                let width_u8: u8 = width_usize.try_into().unwrap();
                println!("position start: {}", pos_start);
                println!("position end: {}", pend);
                let width: NonZeroU8 = width_u8.try_into().unwrap();
                Position::new(contig, pstart, width)
            }
            None => {
                // not sure why there is an option in the noodles here...
                Position::new(contig,
                              pstart,
                              NonZeroU8::new(1).unwrap())
            },
        };
        positions.push(current_position);
    }
    positions
}

pub fn read_variant_positions_from_bed4(filename: &PathBuf) -> Vec<VariantPosition> {
    let reader_result = Builder::<4>::default()
                .build_from_path(filename);
    let mut reader = reader_result.unwrap();
    let mut record = bed::Record::default();
    let mut positions = Vec::new();

    while reader.read_record(&mut record).unwrap() != 0 {
        let contig_raw_name = record.reference_sequence_name();
        let mutstring = record.name().unwrap();
        let contig = Contig::new_from_noodles(contig_raw_name, None);
        // this is already 1-based
        let pos_start = record.feature_start().unwrap().get();
        let pstart: NonZeroUsize = pos_start.try_into().unwrap();
        let pos_end_option = record.feature_end();
        if pos_end_option.is_none() {
            panic!("Variant positions should have a width of at least 1");
        }
        let pend = pos_end_option.unwrap().unwrap().get();
        let width_usize = pend - pos_start + 1;
        let width_u8: u8 = width_usize.try_into().unwrap();
        let width: NonZeroU8 = width_u8.try_into().unwrap();
        let position = Position::new(contig, pstart, width);
        let current_position = VariantPosition::new_from_mut(position, mutstring);
        positions.push(current_position);
    }
    positions
}


pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}

