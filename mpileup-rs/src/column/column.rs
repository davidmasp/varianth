// the ref_base can be N that will be the result when the
// option is None
//
// The depth will be the lenght of the records
// these are re-imported from the io file. can we isolate them there?
// maybe with a type?



// use crate::column::errors::NotValidBam;

use noodles::bam;
use bam::record::Record;
use noodles::core::Position;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::Header;

use crate::column::iterators::MpileupColumnIterator;

use std::fs::File;
use noodles::fasta::indexed_reader::IndexedReader;
use noodles::fasta::io::BufReader;
use noodles::core::Region;
use noodles::fasta::record::Sequence;


pub struct MpileupColumn {
    pub contig: usize,
    pub position: Position, // 1-based
    pub ref_base: Option<u8>,
    pub records: Vec<Record>, // so I would like to use ref here so I can avoid cloning, will try later
    pub cigar_kind: Vec<Kind>,
    pub read_position: Vec<usize>, // might regret this
    pub insertion_size: Vec<Option<usize>>,
}

impl MpileupColumn {
    pub fn add_record_to_column(
        &mut self,
        record: &Record,
        cigar_kind: &Kind,
        read_pos: &usize,
        insertion_size: Option<&usize>,
    ) {
        let isize = match insertion_size {
            Some(x) => Some(x.clone()),
            None => None,
        };
        self.records.push(record.clone());
        self.cigar_kind.push(cigar_kind.clone());
        self.read_position.push(read_pos.clone());
        self.insertion_size.push(isize);
    }
}

pub fn build_mpileupcolumn(
    contig: &usize,
    refposition: &Position,
    record: &Record,
    cigar_kind: &Kind,
    read_pos: &usize,
    insertion_size: Option<&usize>,
) -> MpileupColumn {
    let isize = match insertion_size {
        Some(x) => Some(x.clone()),
        None => None,
    };
    MpileupColumn {
        contig: contig.clone(),
        position: refposition.clone(),
        ref_base: None,
        records: vec![record.clone()],
        cigar_kind: vec![cigar_kind.clone()],
        read_position: vec![read_pos.clone()],
        insertion_size: vec![isize],
    }
}

/// FORMAT
fn format_cols(record: &Record, kind: &Kind, read_position: &usize, isize_opt: Option<&usize>, reference_base: Option<u8>) -> Option<(Vec<u8>, Option<u8>)> {
    let strand_val = record.flags().is_reverse_complemented();
    
    let out_seq: Vec<u8> = match kind {
        Kind::Match
        | Kind::SequenceMatch
        | Kind::SequenceMismatch => {
            let seq_val = record.sequence().iter().nth(*read_position);
            if let Some(refbase) = reference_base {
                match seq_val {
                    Some(x) => {
                        if refbase.to_ascii_uppercase() != x {
                            match strand_val {
                                false => vec![x],
                                true => vec![x.to_ascii_lowercase()],
                            }
                        } else {
                            match strand_val {
                                false => vec![b'.'],
                                true => vec![b','],
                            }
                        }

                    } 
                    None => vec![b'N'],
                }
            } else {
                match seq_val {
                    Some(x) => {
                        match strand_val {
                            false => vec![x],
                            true => vec![x.to_ascii_lowercase()],
                        }
                    } 
                    None => vec![b'N'],
                }
            }

        },
        Kind::Insertion => {
            let isize = isize_opt.unwrap();
            let isize_str = format!("{}", isize);
            let mut seq_full_ins = vec![b'+'];
            seq_full_ins.extend(isize_str.into_bytes());
            record.sequence().iter().skip(*read_position).take(*isize).for_each(|x| {
                let stval = match strand_val {
                    false => x,
                    true => x.to_ascii_lowercase(),
                };
                seq_full_ins.push(stval);
            });
            seq_full_ins
        },
        Kind::Deletion => {
            match strand_val {
                false => vec![b'*'],
                true => vec![b'#'],
            }
        },
        _ => {
            return None
        },
    };

    let bqual = match kind {
        Kind::Insertion => {
            None
        }, 
        _ => {
            let qual_val = record.quality_scores();
            let bqual = qual_val.as_ref().get(*read_position).unwrap().checked_add(33).unwrap();     
            Some(bqual)
        }
    };
    Some((out_seq , bqual))
}

fn query_fa_position(fasta_reader: &mut IndexedReader<BufReader<File>>,
                     contig_name: &[u8],
                     pos: Position) -> u8 {
    let region = Region::new(contig_name, pos..=pos);
    let fa_record = fasta_reader.query(&region).unwrap();
    let fa_seq: &Sequence = fa_record.sequence();
    let tt = fa_seq.as_ref().get(0).unwrap();
    tt.clone()
}

impl MpileupColumn {
    pub fn get_default_depth(&self) -> u32 {
        let depth = self.cigar_kind.iter().filter(|x| {
                match x {
                    Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Deletion => true,
                    _ => false,
                }
            }).count() as u32;
        depth
    }
    pub fn add_reference_base(&mut self,
                              header: &Header,
                              fasta_reader: &mut IndexedReader<BufReader<File>>) -> () {
        let (seqname, _extra) = header.reference_sequences().get_index(self.contig).unwrap();
        let contig_name: &[u8] = seqname.as_ref();
        let seq = query_fa_position(fasta_reader, contig_name, self.position);
        self.ref_base = Some(seq);
    }
    pub fn format_column(&self,
                        header: &Header) -> String {
        let (seqname, _extra) = header.reference_sequences().get_index(self.contig).unwrap();

        let ref_base = match self.ref_base {
            Some(x) => char::from(x),
            None => 'N',
        };

        // PARAM: INTRONS COUNT FOR DEPTH????
        let depth = self.get_default_depth();

        let iterator = MpileupColumnIterator::new(self);
        // i am not sure if this use of strings is the most optimal, 
        // a solution is to use a Vec<u8> and then convert it to a string
        // only one time. The thing is that I am not sure if I can represent
        // the + and - for insertion and deletion with u8, I think it should 
        // initial test indicate that has increased the time but it could also
        // be because i am now including insertions in the string.
        let mut seq_arr: Vec<u8> = Vec::new();
        let mut bqual_arr: Vec<u8> = Vec::new();
        iterator.for_each(|(record, kind_ref, read_position, isize_opt)| {
            let seq_opt = format_cols(record, kind_ref, read_position, isize_opt, self.ref_base);
            match seq_opt {
                Some((seq_val, qual_val_opt)) => {
                    seq_arr.extend(seq_val);
                    if let Some(qual_val) =  qual_val_opt {
                        bqual_arr.push(qual_val);
                    }
                },
                None => (),
            }});
        
        let str_out = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            seqname.to_string(),
            self.position.to_string(),
            ref_base,
            depth.to_string(),
            seq_arr.iter().map(|x| *x as char).collect::<String>(),
            bqual_arr.iter().map(|x| *x as char).collect::<String>(),
        );
        str_out
    }
}


