

use std::path::PathBuf;
use serde::{Serialize, Deserialize};
use regex::Regex;

// so I need to import the trait Read here because it is needed for 
// the read_exact method
use std::io::Write;

use noodles::fasta;
use fasta::io::reader::Builder;
use fasta::io::indexed_reader::Builder as IndexedBuilder;
use noodles::core::{Region, Position};


// this is for fast hashing
use rustc_hash::FxHashMap;

// internal dependencies
// use varianth_core::position::Contig;

use log::{error, info};

/*
STRUCTS
*/

#[derive(Serialize, Deserialize)]
struct TotalCount {
    seqnames: Vec<String>,
    counts: FxHashMap<String, FxHashMap<String, usize>>,
}

#[derive(Serialize, Deserialize)]
struct AggregatedCount {
    counts: FxHashMap<String, usize>,
}

/*
RUNS
*/

pub fn run (
    fasta_path: PathBuf, 
    size: usize, 
    regions_str: Option<String>, 
    regions_path: Option<PathBuf>, 
    output: Option<PathBuf>, 
    table_size: Option<usize>,
    verbose: bool
) -> Result<(), std::io::Error>{

    let tbsize = match table_size {
        Some(table_size) => {
            if verbose {
                info!("Setting table size: {}", table_size);
            }
            table_size
        },
        None => {
            4_u32.pow(size as u32) as usize
        }
    };

    let run_result  = match (regions_str, regions_path) {
        (Some(rst), None) => {
            let regions = region_string_to_vec(&rst).expect("Error parsing the regions");
            let res = run_indexed(fasta_path, regions, size, output, tbsize, verbose);
            res
        },
        (None, Some(_)) => {
            panic!("Not implemented");
        },
        (Some(_), Some(_)) => {
            error!("Not possible to use both regions and regions file at the same time");
            std::process::exit(1);
        },
        (None, None) => {
            let res = run_unindexed(fasta_path, size, output, tbsize, verbose);
            res
        }
    };

    run_result
}

fn region_string_to_vec(regions: &str) -> Result<Vec<Region>, std::num::TryFromIntError> {
    // we should set this as static.
    let re = Regex::new(r"([^,:.]*):(\d+)-(\d+)").expect("Error compiling the regex");

    let regions_vec: Result<Vec<Region>,std::num::TryFromIntError>  = re.captures_iter(regions)
        .map(|x| x.extract() )
        .map(|(_, [contig, start_str, end_str])| {
            let start_u32 = start_str.parse::<usize>().expect("Error parsing start of the region");
            let end_u32 = end_str.parse::<usize>().expect("Error parsing end of the region");
            let start = Position::try_from(start_u32)?;
            let end = Position::try_from(end_u32)?;
            Ok(Region::new(contig, start..=end))
        })
        .collect();

    regions_vec
}

fn run_indexed(fasta_path: PathBuf, regions: Vec<Region>, size: usize, output: Option<PathBuf>, table_size: usize, verbose: bool) -> Result<(), std::io::Error> {

    if verbose {
        info!("Running indexed mode, total of {} regions", regions.len());
    }

    // fn cte
    let mut fa = IndexedBuilder::default()
                 .build_from_path(fasta_path)?;
    
    let mut agg_count = AggregatedCount {
        counts: FxHashMap::with_capacity_and_hasher(table_size, Default::default()),
    };

    let mut table_vec: Vec<usize> = vec![0; table_size];
    // i am not sure if this is ideal? 
    // currently a work around
    let mut vec_sequences: Vec<Vec<u8>> = Vec::new();
    for rr in regions.iter() {
        // I am not sure why I can't use fa.read_sequence(); here
        // I don't seem to be able to find a seek method
        // in addition to that this would be a problem because I need
        // to know the length of the sequence to be able to read it,
        // not impossible

        let fasta_record_res = fa.query(rr);

        let fasta_record = match fasta_record_res {
            Ok(ff) => {ff},
            Err(e) => {
                error!("Error reading the region: {}", e);
                std::process::exit(1);
            }
        };
        // is this putting the sequence obj in memory?
        let seq = fasta_record.sequence().as_ref().to_vec();

        vec_sequences.push(seq);
    }

    for seq in vec_sequences.iter() {
        update_table(&mut table_vec, seq, size);
    }

    for (idx, count) in table_vec.iter().enumerate() {
        let kmer_string = index_to_string(idx, size);
        agg_count.counts.insert(kmer_string, *count);
    }

    let _ = serialize_to_json(output, &agg_count);

    Ok(())
}    

fn run_unindexed (fasta_path: PathBuf, size: usize, output: Option<PathBuf>, table_size: usize, verbose: bool) -> Result<(), std::io::Error> {
// fn cte

    let mut total_count = TotalCount {
        seqnames: Vec::new(),
        counts: FxHashMap::with_capacity_and_hasher(table_size, Default::default()),
    };

    let mut fa = Builder::default()
        .build_from_path(fasta_path)?;

    loop {
        let mut string_contig_name: String = Default::default();
        let bytes_read = fa.read_definition(&mut string_contig_name)?;
        if bytes_read == 0 {
            break;
        }
        // this should remove the >, is this always true?
        string_contig_name.remove(0);
        if verbose {
            info!("Processing contig: {}", string_contig_name);
        }
        total_count.seqnames.push(string_contig_name.clone());

        // sequence fun
        let mut table_vec: Vec<usize> = vec![0; table_size];
        let mut sequence_buf = Vec::new();
        let _ = fa.read_sequence(&mut sequence_buf)?;
        update_table(&mut table_vec, &sequence_buf, size);

        let mut hash_table_string = FxHashMap::with_capacity_and_hasher(table_size, Default::default());
        
        for (idx, count) in table_vec.iter().enumerate() {
            let kmer_string = index_to_string(idx, size);
            hash_table_string.insert(kmer_string, *count);
        }

        total_count.counts.insert(string_contig_name.clone(), hash_table_string);
        if verbose {
            info!("Contig {} processed", string_contig_name);
        }
    }

    if verbose {
        info!("Seqnames: {:?}", total_count.seqnames);
    }

    let _ = serialize_to_json(output, &total_count);

    Ok(())
}

fn slice_to_index(kmer: &[u8]) -> usize {
    let mut hash_val = 0;
    let mut count = 0;
    for &byte in kmer.iter() {
        match byte {
            65 => {
                hash_val = hash_val + 4_usize.pow(count) * 0
            }, 
            67 => {
                hash_val = hash_val + 4_usize.pow(count) * 1
            },
            71 => {
                hash_val = hash_val + 4_usize.pow(count) * 2
            },
            84 => {
                hash_val = hash_val + 4_usize.pow(count) * 3
            },
            _ => {
                panic!("Invalid character in sequence");
            }
        }
        count += 1;
    }
    hash_val
}

fn index_to_string(idx: usize, ksize: usize) -> String {
    
    let mut kmer_chars =  Vec::with_capacity(ksize);
    let max_pos = ksize - 1;
    let mut idx = idx;

    for i in (0..=max_pos).rev() {
        // starts with max_pos
        let intdiv = idx / 4_usize.pow(i as u32);
        // this is the remainder
        idx = idx % 4_usize.pow(i as u32);

        let base = match intdiv {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!("Invalid index"),
        };
        kmer_chars.push(base);
    }
    let kmer_string: String = kmer_chars.into_iter().rev().collect();
    kmer_string
}

fn update_table(table: &mut Vec<usize>, sequence_buf: &Vec<u8>, ksize: usize) -> () {
    let mut cursor = 0;
    let mut cend = ksize;
    while cend <= sequence_buf.len() {
        let kmer: &[u8] = &sequence_buf[cursor..cend];
        let seq_idx = slice_to_index(kmer);
        table[seq_idx] += 1;
        cursor += 1;
        cend += 1;
    }
}

fn serialize_to_json<T: Serialize>(json_path: Option<PathBuf>, obj: &T) -> () {
    let total_count_json = serde_json::to_string(obj);

    let json_str = match total_count_json {
        Ok(json) => {
            json
        },
        Err(e) => {
            // There was an error serializing total_count
            error!("Error generating the json output: {}", e);
            std::process::exit(1);
        }
    };

    match json_path {
        Some(output_path) => {
            let mut output_file = std::fs::File::create(output_path).expect("Error creating the output file");
            output_file.write_all(json_str.as_bytes()).expect("kdfjsk");
        },
        None => {
            println!("{}", json_str);
        }
    }
}

