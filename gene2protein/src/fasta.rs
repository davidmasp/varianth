use anyhow::{Context, Result, bail};
use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::fasta::fai;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use bstr::BString;

pub fn reverse_complement(seq: &BString) -> Result<BString> {
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
                bail!("non-canonical DNA base '{}' at position {}", b as char, idx);
            }
        };
        rc.push(comp);
    }
    Ok(BString::from(rc))
}

/// Open an indexed FASTA reader.
///
/// If `fai_path` is `None`, the index path is derived as `<fasta>.fai`
/// (e.g. `genome.fa` → `genome.fa.fai`).
pub fn open_indexed_fasta<P: AsRef<Path>>(
    fasta_path: P,
    fai_path: Option<P>,
) -> Result<fasta::io::IndexedReader<BufReader<File>>> {
    let fasta_path = fasta_path.as_ref();

    let fai_path: PathBuf = match fai_path {
        Some(p) => p.as_ref().to_path_buf(),
        None => {
            let mut derived = fasta_path.as_os_str().to_os_string();
            derived.push(".fai");
            let derived = PathBuf::from(derived);
            if !derived.exists() {
                bail!(
                    "expected FAI index at '{}' but it does not exist; run `samtools faidx` first",
                    derived.display()
                );
            }
            derived
        }
    };

    let index = File::open(&fai_path)
        .map(BufReader::new)
        .map(fai::io::Reader::new)
        .with_context(|| format!("failed to open FAI file '{}'", fai_path.display()))?
        .read_index()
        .with_context(|| format!("failed to read FAI index '{}'", fai_path.display()))?;

    File::open(fasta_path)
        .map(BufReader::new)
        .map(|r| fasta::io::IndexedReader::new(r, index))
        .with_context(|| format!("failed to open FASTA file '{}'", fasta_path.display()))
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
) -> Result<Vec<u8>> {
    let fai_record = index
        .as_ref()
        .iter()
        .find(|record| record.name() == reference_name)
        .with_context(|| format!("reference not found in FASTA index: {reference_name}"))?;

    let start = Position::try_from(1usize).context("failed to build FASTA start position")?;
    let sequence_length = usize::try_from(fai_record.length()).with_context(|| {
        format!("reference length for {reference_name} does not fit into usize")
    })?;
    let end = Position::try_from(sequence_length).with_context(|| {
        format!("failed to build FASTA end position for reference {reference_name}")
    })?;

    let region = Region::new(reference_name, start..=end);
    let record = reader
        .query(&region)
        .with_context(|| format!("failed to query full FASTA record {reference_name}"))?;

    Ok(record.sequence().as_ref().to_vec())
}
