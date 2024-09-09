
use noodles::bam;
use bam::io::Reader as BamReader;
use bam::record::Record;
use noodles::bam::io::reader::Query;
use noodles::bgzf::Reader as GzReader;
use noodles::sam::alignment::record::MappingQuality;
use noodles::core::Position;
use noodles::sam::alignment::record::cigar::op::Kind;

use std::fs::File;
use crate::column::operations::process_read;
use std::collections::BTreeMap;

use crate::column::column::MpileupColumn;


pub struct MpileupColumnIterator<'a> {
    mpileup_column: &'a MpileupColumn,
    current_index: usize,
}

impl<'a> MpileupColumnIterator<'a> {
    pub fn new(mpileup_column: &'a MpileupColumn) -> Self {
        MpileupColumnIterator {
            mpileup_column,
            current_index: 0,
        }
    }
}

impl<'a> Iterator for MpileupColumnIterator<'a> {
    type Item = (&'a Record,&'a Kind,&'a usize, Option<&'a usize>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_index < self.mpileup_column.records.len() {
            let record = &self.mpileup_column.records[self.current_index];
            let read_position = &self.mpileup_column.read_position[self.current_index];
            let kind_ref = &self.mpileup_column.cigar_kind[self.current_index];
            let isize_opt = match &self.mpileup_column.insertion_size[self.current_index] {
                Some(x) => Some(x),
                None => None,
            };
            self.current_index += 1;
            Some((record, kind_ref, read_position, isize_opt))
        } else {
            None // End of iteration.
        }
    }
}

impl<'a> Iterator for MpileupIterator {
    type Item = MpileupColumn;

    fn next(&mut self) -> Option<Self::Item> {
        // here we check if the buffer is not empty, and thus should
        // return the last element in the buffer.
        if let Some(item) = self.buffer.pop() {
            return Some(item);
        }

        let mut record_obj = bam::Record::default();
        // note here about the usize
        //  block size of 0 is returned, the stream reached EOF.
        // https://docs.rs/noodles-bam/latest/noodles_bam/io/reader/struct.Reader.html#method.read_record
        while let Ok(ret_usize) = self.reader.read_record(&mut record_obj) {
            if ret_usize == 0 {
                break;
            } else {
                if record_obj.mapping_quality().expect("Mapping qual not set?") < self.min_mapq {
                    continue;
                }
                let _ = process_read(&mut self.bt_pileup,
                     &record_obj,
                      self.exclude_flags,
                       self.include_flags_option);
                if self.bt_pileup.len() > self.cache_size {
                    // push to buffer
                    push_buffer(&mut self.bt_pileup,
             Some(&record_obj.alignment_start().unwrap().unwrap()),
                         &mut self.buffer);
                    // return the last element
                    if let Some(item) = self.buffer.pop() {
                        return Some(item);
                    }
                }
            }
        }
        // this marks the first 
        push_buffer(&mut self.bt_pileup,
            None,
                        &mut self.buffer);
        if let Some(item) = self.buffer.pop() {
            return Some(item);
        } else {
            None
        }
    }
}

// this is the important one:
pub struct MpileupIterator {
    reader: BamReader<GzReader<File>>,
    min_mapq: MappingQuality,
    exclude_flags: u16,
    include_flags_option: Option<u16>,
    bt_pileup: BTreeMap<Position, MpileupColumn>,
    cache_size: usize,
    buffer: Vec<MpileupColumn>,
}

impl MpileupIterator {
    pub fn new(reader: BamReader<GzReader<File>>,
               min_mapq: MappingQuality,
               exclude_flags: u16, 
               include_flags_option: Option<u16>,
            ) -> Result<Self, std::io::Error> {
        let bt_pileup: BTreeMap<Position, MpileupColumn> = BTreeMap::new();
        Ok(Self {
            reader,
            min_mapq,
            exclude_flags,
            include_flags_option,
            bt_pileup: bt_pileup,
            cache_size: 500,
            buffer: Vec::new(),
        })
    }
}

/// Queried version, is there a way to merge them??
impl<'a> Iterator for MpileupQueriedIterator<'a> {
    type Item = MpileupColumn;

    fn next(&mut self) -> Option<Self::Item> {
        // here we check if the buffer is not empty, and thus should
        // return the last element in the buffer.
        if let Some(item) = self.buffer.pop() {
            return Some(item);
        }

        while let Some(record_obj_res) = self.reader.next() {
            let record_obj = record_obj_res.expect("Error reading record");
            if record_obj.mapping_quality().expect("Mapping qual not set?") < self.min_mapq {
                continue;
            }
            let _ = process_read(&mut self.bt_pileup,
                    &record_obj,
                    self.exclude_flags,
                    self.include_flags_option);
            if self.bt_pileup.len() > self.cache_size {
                // push to buffer
                push_buffer(&mut self.bt_pileup,
            Some(&record_obj.alignment_start().unwrap().unwrap()),
                        &mut self.buffer);
                // return the last element
                if let Some(item) = self.buffer.pop() {
                    return Some(item);
                }
            }
            
        }
        // this marks the first 
        push_buffer(&mut self.bt_pileup,
            None,
                        &mut self.buffer);
        if let Some(item) = self.buffer.pop() {
            return Some(item);
        } else {
            None
        }
    }
}


// this is the important one:
pub struct MpileupQueriedIterator<'a> {
    reader: Query<'a , GzReader<File>>,
    min_mapq: MappingQuality,
    exclude_flags: u16,
    include_flags_option: Option<u16>,
    bt_pileup: BTreeMap<Position, MpileupColumn>,
    cache_size: usize,
    buffer: Vec<MpileupColumn>,
}

impl<'a> MpileupQueriedIterator<'a> {
    pub fn new(reader: Query<'a , GzReader<File>>,
               min_mapq: MappingQuality,
               exclude_flags: u16, 
               include_flags_option: Option<u16>,
            ) -> Result<Self, std::io::Error> {
        let bt_pileup: BTreeMap<Position, MpileupColumn> = BTreeMap::new();
        Ok(Self {
            reader,
            min_mapq,
            exclude_flags,
            include_flags_option,
            bt_pileup: bt_pileup,
            cache_size: 500,
            buffer: Vec::new(),
        })
    }
}

pub fn push_buffer(
    btmap: &mut BTreeMap<Position, MpileupColumn>,
    current_position_option: Option<&Position>,
    buffer: &mut Vec<MpileupColumn>
) -> () {
    let mut keys_to_delete = Vec::new();

    let _ = match current_position_option {
        Some(current_position) => {
            for (keypos, _column) in btmap.iter() {
                if keypos < current_position {
                    keys_to_delete.push(*keypos);
                } else {
                    continue;
                }
            }
        }
        None => {
            for (keypos, _column) in btmap.iter() {
                keys_to_delete.push(*keypos);
            }
        }
    };

    // reorder keys_to_delete
    //keys_to_delete.sort(); so I am pretty sure i don't need to sort here.
    keys_to_delete.reverse();

    // delete the keys
    for key in keys_to_delete {
        if let Some(column) = btmap.remove(&key) {
            buffer.push(column);
        }
    }

    // should deal with errors at some point!
    ()
}

