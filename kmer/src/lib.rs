use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::PathBuf;

use fasta::io::indexed_reader::Builder as IndexedBuilder;
use fasta::io::reader::Builder;
use log::info;
use noodles::bed;
use noodles::core::{Position, Region};
use noodles::fasta;
use regex::Regex;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct TotalCount {
    seqnames: Vec<String>,
    counts: FxHashMap<String, FxHashMap<String, usize>>,
}

#[derive(Serialize, Deserialize)]
struct AggregatedCount {
    counts: FxHashMap<String, usize>,
}

pub fn run(
    fasta_path: PathBuf,
    size: usize,
    regions_str: Option<String>,
    regions_path: Option<PathBuf>,
    output: Option<PathBuf>,
    table_size: Option<usize>,
    verbose: bool,
    skip_ambiguous: bool,
) -> io::Result<()> {
    let table_size = match table_size {
        Some(table_size) => {
            if verbose {
                info!("Setting table size: {}", table_size);
            }
            table_size
        }
        None => 4_u32.pow(size as u32) as usize,
    };

    match (regions_str, regions_path) {
        (Some(rst), None) => {
            let regions = region_string_to_vec(&rst)?;
            run_indexed(
                fasta_path,
                regions,
                size,
                output,
                table_size,
                verbose,
                skip_ambiguous,
            )
        }
        (None, Some(regions_path)) => {
            let regions = regions_file_to_vec(regions_path)?;
            run_indexed(
                fasta_path,
                regions,
                size,
                output,
                table_size,
                verbose,
                skip_ambiguous,
            )
        }
        (Some(_), Some(_)) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "not possible to use both regions and regions file at the same time",
        )),
        (None, None) => run_unindexed(
            fasta_path,
            size,
            output,
            table_size,
            verbose,
            skip_ambiguous,
        ),
    }
}

fn region_string_to_vec(regions: &str) -> io::Result<Vec<Region>> {
    let re = Regex::new(r"([^,:.]*):(\d+)-(\d+)").map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("error compiling region regex: {e}"),
        )
    })?;

    let regions_vec: io::Result<Vec<Region>> = re
        .captures_iter(regions)
        .map(|x| x.extract())
        .map(|(_, [contig, start_str, end_str])| {
            let start_u32 = start_str.parse::<usize>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("error parsing start of the region: {e}"),
                )
            })?;
            let end_u32 = end_str.parse::<usize>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("error parsing end of the region: {e}"),
                )
            })?;
            let start = Position::try_from(start_u32).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid region start: {e}"),
                )
            })?;
            let end = Position::try_from(end_u32).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid region end: {e}"),
                )
            })?;
            Ok(Region::new(contig, start..=end))
        })
        .collect();

    let regions_vec = regions_vec?;

    if regions_vec.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no regions found in region string",
        ));
    }

    Ok(regions_vec)
}

fn regions_file_to_vec(regions_path: PathBuf) -> io::Result<Vec<Region>> {
    let mut reader = File::open(regions_path)
        .map(BufReader::new)
        .map(bed::io::Reader::<3, _>::new)?;

    let mut regions = Vec::new();

    loop {
        let mut record = bed::Record::<3>::default();
        let bytes_read = reader.read_record(&mut record)?;

        if bytes_read == 0 {
            break;
        }

        regions.push(bed_record_to_region(record)?);
    }

    if regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no regions found in BED file",
        ));
    }

    Ok(regions)
}

fn bed_record_to_region(record: bed::Record<3>) -> io::Result<Region> {
    let start = record.feature_start()?;
    let end = record.feature_end().transpose()?.ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "BED record has no end position")
    })?;

    Ok(Region::new(
        record.reference_sequence_name().to_owned(),
        start..=end,
    ))
}

fn run_indexed(
    fasta_path: PathBuf,
    regions: Vec<Region>,
    size: usize,
    output: Option<PathBuf>,
    table_size: usize,
    verbose: bool,
    skip_ambiguous: bool,
) -> io::Result<()> {
    if verbose {
        info!("Running indexed mode, total of {} regions", regions.len());
    }

    let mut fa = IndexedBuilder::default().build_from_path(fasta_path)?;

    let mut agg_count = AggregatedCount {
        counts: FxHashMap::with_capacity_and_hasher(table_size, Default::default()),
    };

    let mut table_vec: Vec<usize> = vec![0; table_size];

    for region in regions.iter() {
        let fasta_record = fa.query(region)?;
        let seq = fasta_record.sequence().as_ref().to_vec();
        update_table(&mut table_vec, &seq, size, skip_ambiguous)?;
    }

    for (idx, count) in table_vec.iter().enumerate() {
        let kmer_string = index_to_string(idx, size);
        agg_count.counts.insert(kmer_string, *count);
    }

    serialize_to_json(output, &agg_count)
}

fn run_unindexed(
    fasta_path: PathBuf,
    size: usize,
    output: Option<PathBuf>,
    table_size: usize,
    verbose: bool,
    skip_ambiguous: bool,
) -> io::Result<()> {
    let mut total_count = TotalCount {
        seqnames: Vec::new(),
        counts: FxHashMap::with_capacity_and_hasher(table_size, Default::default()),
    };

    let mut fa = Builder::default().build_from_path(fasta_path)?;

    for result in fa.records() {
        let record = result?;
        let string_contig_name = String::from_utf8(record.name().to_vec()).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FASTA record name is not valid UTF-8: {e}"),
            )
        })?;

        if verbose {
            info!("Processing contig: {}", string_contig_name);
        }
        total_count.seqnames.push(string_contig_name.clone());

        let mut table_vec: Vec<usize> = vec![0; table_size];
        update_table(
            &mut table_vec,
            record.sequence().as_ref(),
            size,
            skip_ambiguous,
        )?;

        let mut hash_table_string =
            FxHashMap::with_capacity_and_hasher(table_size, Default::default());

        for (idx, count) in table_vec.iter().enumerate() {
            let kmer_string = index_to_string(idx, size);
            hash_table_string.insert(kmer_string, *count);
        }

        total_count
            .counts
            .insert(string_contig_name.clone(), hash_table_string);
        if verbose {
            info!("Contig {} processed", string_contig_name);
        }
    }

    if verbose {
        info!("Seqnames: {:?}", total_count.seqnames);
    }

    serialize_to_json(output, &total_count)
}

fn slice_to_index(kmer: &[u8], skip_ambiguous: bool) -> io::Result<Option<usize>> {
    let mut hash_val = 0;
    for (count, &byte) in kmer.iter().enumerate() {
        match byte {
            b'A' | b'a' => {
                hash_val += 4_usize.pow(count as u32) * 0;
            }
            b'C' | b'c' => {
                hash_val += 4_usize.pow(count as u32) * 1;
            }
            b'G' | b'g' => {
                hash_val += 4_usize.pow(count as u32) * 2;
            }
            b'T' | b't' => {
                hash_val += 4_usize.pow(count as u32) * 3;
            }
            b'N' | b'n' => {
                if skip_ambiguous {
                    return Ok(None);
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "ambiguous base 'N' found in k-mer; use --skip-ambiguous to skip such k-mers",
                    ));
                }
            }
            _ => {
                if skip_ambiguous {
                    return Ok(None);
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "invalid or ambiguous character '{}' (byte {}) found in k-mer; use --skip-ambiguous to skip such k-mers",
                            byte as char, byte
                        ),
                    ));
                }
            }
        }
    }
    Ok(Some(hash_val))
}

fn index_to_string(idx: usize, ksize: usize) -> String {
    let mut kmer_chars = Vec::with_capacity(ksize);
    let max_pos = ksize - 1;
    let mut idx = idx;

    for i in (0..=max_pos).rev() {
        let intdiv = idx / 4_usize.pow(i as u32);
        idx %= 4_usize.pow(i as u32);

        let base = match intdiv {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!("Invalid index"),
        };
        kmer_chars.push(base);
    }
    kmer_chars.into_iter().rev().collect()
}

fn update_table(
    table: &mut Vec<usize>,
    sequence_buf: &[u8],
    ksize: usize,
    skip_ambiguous: bool,
) -> io::Result<()> {
    let mut cursor = 0;
    let mut cend = ksize;
    while cend <= sequence_buf.len() {
        let kmer: &[u8] = &sequence_buf[cursor..cend];
        if let Some(seq_idx) = slice_to_index(kmer, skip_ambiguous)? {
            table[seq_idx] += 1;
        }
        cursor += 1;
        cend += 1;
    }
    Ok(())
}

fn serialize_to_json<T: Serialize>(json_path: Option<PathBuf>, obj: &T) -> io::Result<()> {
    let json_str = serde_json::to_string(obj).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("error generating the json output: {e}"),
        )
    })?;

    match json_path {
        Some(output_path) => {
            let mut output_file = File::create(output_path)?;
            output_file.write_all(json_str.as_bytes())?;
        }
        None => {
            println!("{}", json_str);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_bed3_record(src: &str) -> io::Result<bed::Record<3>> {
        let mut reader = bed::io::Reader::<3, _>::new(src.as_bytes());
        let mut record = bed::Record::<3>::default();
        reader.read_record(&mut record)?;
        Ok(record)
    }

    #[test]
    fn bed_record_to_region_converts_bed_coordinates() -> io::Result<()> {
        let record = parse_bed3_record("chr1\t0\t4")?;
        let region = bed_record_to_region(record)?;

        assert_eq!(region.to_string(), "chr1:1-4");

        Ok(())
    }

    #[test]
    fn bed_record_to_region_accepts_extra_fields() -> io::Result<()> {
        let record = parse_bed3_record("chr1\t1\t5\tname\t42")?;
        let region = bed_record_to_region(record)?;

        assert_eq!(region.to_string(), "chr1:2-5");

        Ok(())
    }

    #[test]
    fn invalid_bed_record_returns_invalid_data() {
        let record = parse_bed3_record("chr1\tbad\t5").unwrap();
        let result = bed_record_to_region(record);

        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn slice_to_index_skips_ambiguous_bases() -> io::Result<()> {
        assert_eq!(slice_to_index(b"AN", true)?, None);
        assert!(slice_to_index(b"AN", false).is_err());
        Ok(())
    }
}
