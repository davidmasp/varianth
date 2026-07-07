use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::fasta::fai;
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use bstr::BString;

pub fn reverse_complement(seq: &BString) -> Result<BString, String> {
    let mut rc = Vec::with_capacity(seq.len());
    for (idx, b) in seq.as_slice().iter().copied().enumerate().rev() {
        let comp = match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            _ => {
                return Err(format!(
                    "Non-canonical DNA base '{}' at position {}",
                    b as char, idx
                ));
            }
        };
        rc.push(comp);
    }
    Ok(BString::from(rc))
}

/// Open an indexed FASTA reader.
///
/// If `fai_path` is `None`, the index path is derived as `<fasta>.fai`
/// (e.g. `genome.fa` → `genome.fa.fai`). The function panics if the
/// derived index file does not exist.
pub fn open_indexed_fasta<P: AsRef<Path>>(
    fasta_path: P,
    fai_path: Option<P>,
) -> fasta::io::IndexedReader<BufReader<File>> {
    let fasta_path = fasta_path.as_ref();

    let fai_path: PathBuf = match fai_path {
        Some(p) => p.as_ref().to_path_buf(),
        None => {
            let derived = fasta_path
                .to_str()
                .expect("fasta path is not valid UTF-8")
                .to_string()
                + ".fai";
            let derived = PathBuf::from(derived);
            assert!(
                derived.exists(),
                "expected FAI index at '{}' but it does not exist; run `samtools faidx` first",
                derived.display()
            );
            derived
        }
    };

    let index = File::open(&fai_path)
        .map(BufReader::new)
        .map(fai::io::Reader::new)
        .unwrap_or_else(|e| panic!("failed to open FAI file '{}': {e}", fai_path.display()))
        .read_index()
        .unwrap_or_else(|e| panic!("failed to read FAI index '{}': {e}", fai_path.display()));

    File::open(fasta_path)
        .map(BufReader::new)
        .map(|r| fasta::io::IndexedReader::new(r, index))
        .unwrap_or_else(|e| panic!("failed to open FASTA file '{}': {e}", fasta_path.display()))
}

/// Return the number of sequences described in the index of an [`IndexedReader`].
pub fn sequence_count(reader: &fasta::io::IndexedReader<BufReader<File>>) -> usize {
    reader.index().as_ref().len()
}

/// Return all sequence IDs (names) from the index of an [`IndexedReader`].
pub fn sequence_ids(reader: &fasta::io::IndexedReader<BufReader<File>>) -> Vec<String> {
    reader
        .index()
        .as_ref()
        .iter()
        .map(|r| r.name().to_string())
        .collect()
}

pub fn pull_entire_record(
    reader: &mut fasta::io::IndexedReader<BufReader<File>>,
    index: &fai::Index,
    reference_name: &str,
) -> io::Result<Vec<u8>> {
    let fai_record = index
        .as_ref()
        .iter()
        .find(|record| record.name() == reference_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("reference not found in index: {reference_name}"),
            )
        })?;

    let start = Position::try_from(1usize)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    let sequence_length = usize::try_from(fai_record.length())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    let end = Position::try_from(sequence_length)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;

    let region = Region::new(reference_name, start..=end);
    let record = reader.query(&region)?;

    Ok(record.sequence().as_ref().to_vec())
}
