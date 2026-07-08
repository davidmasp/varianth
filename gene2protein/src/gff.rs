use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

use bstr::BString;

// ── GFF field types ───────────────────────────────────────────────────────────

/// Strand as represented in a GFF file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    fn from_str(s: &str) -> Self {
        match s {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

// ── GffRecord ─────────────────────────────────────────────────────────────────

/// A single feature record parsed from a GFF3 file.
///
/// Fields follow the 9-column GFF3 specification:
/// seqid, source, type, start, end, score, strand, phase, attributes.
#[derive(Debug, Clone)]
pub struct GffRecord {
    /// Column 1 – sequence identifier (e.g. chromosome name).
    pub seqid: BString,
    /// Column 2 – origin annotation source (program or database).
    pub source: BString,
    /// Column 3 – feature type (SO term in GFF3, e.g. "gene", "CDS").
    pub feature_type: BString,
    /// Column 4 – 1-based start coordinate.
    pub start: usize,
    /// Column 5 – 1-based end coordinate (inclusive).
    pub end: usize,
    /// Column 6 – numeric score; `None` when the field contains ".".
    pub score: Option<f64>,
    /// Column 7 – strand orientation.
    pub strand: Strand,
    /// Column 8 – reading frame offset (0, 1, or 2) for CDS features; `None` for ".".
    pub phase: Option<u8>,
    /// Column 9 – parsed `key=value` attribute pairs.
    /// Values that contain commas are split into multiple entries in the `Vec`.
    pub attributes: HashMap<String, Vec<String>>,
}

/// Collection of CDS records grouped by protein ID.
#[derive(Debug, Default, Clone)]
pub struct CdsProteome {
    /// Map of `protein_id` to its CDS records.
    pub records: HashMap<String, Vec<GffRecord>>,
}

impl CdsProteome {
    /// Return CDS records for a protein ID, if present, cloned.
    pub fn get_cloned(&self, protein_id: &str) -> Option<Vec<GffRecord>> {
        self.records.get(protein_id).cloned()
    }

    /// Return CDS records for a protein ID, if present, as a borrowed slice.
    pub fn get_ref(&self, protein_id: &str) -> Option<&[GffRecord]> {
        self.records.get(protein_id).map(Vec::as_slice)
    }

    /// Return the number of protein IDs with CDS records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Return all protein IDs available in this proteome.
    pub fn protein_ids(&self) -> HashSet<String> {
        let hs_out: HashSet<String> = self.records.keys().cloned().collect();
        hs_out
    }

    /// Return an iterator over all protein IDs in this proteome.
    pub fn protein_ids_ref(&self) -> impl Iterator<Item = &String> {
        self.records.keys()
    }
}

impl GffRecord {
    /// Convenience getter for the first value of an attribute key.
    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes.get(key)?.first().map(String::as_str)
    }

    /// Parse a raw tab-delimited GFF line into a `GffRecord`.
    fn parse(line: &str) -> Result<Self, io::Error> {
        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 9 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "expected 9 tab-separated fields, got {}: {line}",
                    fields.len()
                ),
            ));
        }

        let start = fields[3].parse::<usize>().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("invalid start: {e}"))
        })?;
        let end = fields[4]
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid end: {e}")))?;

        let score = match fields[5] {
            "." | "" => None,
            s => Some(s.parse::<f64>().map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("invalid score: {e}"))
            })?),
        };

        let phase = match fields[7] {
            "." | "" => None,
            s => Some(s.parse::<u8>().map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("invalid phase: {e}"))
            })?),
        };

        let attributes = parse_attributes(fields[8]);

        Ok(GffRecord {
            seqid: BString::from(fields[0]),
            source: BString::from(fields[1]),
            feature_type: BString::from(fields[2]),
            start,
            end,
            score,
            strand: Strand::from_str(fields[6]),
            phase,
            attributes,
        })
    }
}

/// Parse the GFF3 attribute string (`key=val1,val2;key2=val`) into a map.
fn parse_attributes(raw: &str) -> HashMap<String, Vec<String>> {
    let mut map = HashMap::new();
    for pair in raw.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        if let Some((key, value)) = pair.split_once('=') {
            let values: Vec<String> = value.split(',').map(|v| v.to_string()).collect();
            map.insert(key.to_string(), values);
        }
    }
    map
}

// ── Reader ────────────────────────────────────────────────────────────────────

/// Inner buffer type used by [`GffReader`].
enum GffBuf {
    Plain(BufReader<File>),
    Gz(BufReader<GzDecoder<File>>),
}

/// Streaming GFF3 reader that works with both plain and gzip-compressed files.
///
/// Skips comment lines (starting with `#`) and blank lines automatically.
///
/// # Example
/// ```no_run
/// use gene2protein::GffReader;
///
/// fn main() -> std::io::Result<()> {
///     let reader = GffReader::open("annotations.gff3.gz")?;
///     for record in reader {
///         let record = record?;
///         println!("{}\t{}\t{}", record.seqid, record.feature_type, record.start);
///     }
///     Ok(())
/// }
/// ```
pub struct GffReader {
    buf: GffBuf,
    line: String,
}

impl GffReader {
    /// Open a GFF file for reading.
    ///
    /// The format (plain or gzip) is detected from the `.gz` file extension.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let buf = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
            GffBuf::Gz(BufReader::new(GzDecoder::new(file)))
        } else {
            GffBuf::Plain(BufReader::new(file))
        };
        Ok(GffReader {
            buf,
            line: String::new(),
        })
    }
}

impl Iterator for GffReader {
    type Item = io::Result<GffRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line.clear();
            let bytes_read = match &mut self.buf {
                GffBuf::Plain(r) => r.read_line(&mut self.line),
                GffBuf::Gz(r) => r.read_line(&mut self.line),
            };

            match bytes_read {
                Err(e) => return Some(Err(e)),
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let trimmed = self.line.trim_end_matches(['\n', '\r']);
                    // skip blank lines and comment / pragma lines
                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        continue;
                    }
                    return Some(GffRecord::parse(trimmed));
                }
            }
        }
    }
}

/// Collect CDS records by protein_id for a set of protein IDs.
pub fn collect_cds_by_protein_id<P: AsRef<Path>>(
    gff_path: P,
    prot_ids: &HashSet<String>,
) -> io::Result<CdsProteome> {
    let gff_reader = GffReader::open(gff_path)?;
    let cds_bstring = BString::from("CDS");
    let mut proteome = CdsProteome {
        records: HashMap::new(),
    };

    for record in gff_reader {
        let record = record?;
        if record.feature_type != cds_bstring {
            continue;
        }

        let pid = record.attribute("protein_id").ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "missing protein_id attribute on CDS record",
            )
        })?;

        if !prot_ids.contains(pid) {
            continue;
        }
        log::debug!("Found CDS with protein_id: {}", pid);
        proteome
            .records
            .entry(pid.to_string())
            .or_insert_with(Vec::new)
            .push(record);
    }

    Ok(proteome)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GFF_LINE: &str = "chr1\tensembl_havana\tgene\t3069168\t3438621\t.\t+\t.\tID=ENSG00000142611.17;gene_id=ENSG00000142611.17;gene_type=protein_coding;gene_name=PRDM16";

    #[test]
    fn parse_basic_record() {
        let rec = GffRecord::parse(GFF_LINE).unwrap();
        assert_eq!(rec.seqid, "chr1");
        assert_eq!(rec.source, "ensembl_havana");
        assert_eq!(rec.feature_type, "gene");
        assert_eq!(rec.start, 3069168);
        assert_eq!(rec.end, 3438621);
        assert_eq!(rec.score, None);
        assert_eq!(rec.strand, Strand::Forward);
        assert_eq!(rec.phase, None);
        assert_eq!(rec.attribute("gene_name"), Some("PRDM16"));
        assert_eq!(rec.attribute("gene_type"), Some("protein_coding"));
    }

    #[test]
    fn parse_attributes_multi_value() {
        let line = "chr1\tsrc\tCDS\t1\t10\t.\t+\t0\ttag=MANE_Select,MANE_Plus_Clinical;ID=cds1";
        let rec = GffRecord::parse(line).unwrap();
        let tags = rec.attributes.get("tag").unwrap();
        assert_eq!(tags, &vec!["MANE_Select", "MANE_Plus_Clinical"]);
    }

    #[test]
    fn parse_strand_variants() {
        let make = |strand: &str| -> GffRecord {
            let line = format!("chr1\t.\tgene\t1\t100\t.\t{strand}\t.\t.");
            GffRecord::parse(&line).unwrap()
        };
        assert_eq!(make("+").strand, Strand::Forward);
        assert_eq!(make("-").strand, Strand::Reverse);
        assert_eq!(make(".").strand, Strand::Unknown);
    }

    #[test]
    fn skip_comments_and_blanks() {
        use std::io::{BufReader, Cursor};
        // Build a fake in-memory GFF
        let data = format!("##gff-version 3\n# comment\n\n{GFF_LINE}\n");
        let cursor = Cursor::new(data.into_bytes());
        // Manually iterate using the same logic as GffReader
        let mut reader = BufReader::new(cursor);
        let mut records: Vec<GffRecord> = Vec::new();
        let mut line = String::new();
        loop {
            line.clear();
            let n = std::io::BufRead::read_line(&mut reader, &mut line).unwrap();
            if n == 0 {
                break;
            }
            let t = line.trim_end_matches(['\n', '\r']);
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            records.push(GffRecord::parse(t).unwrap());
        }
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].feature_type, "gene");
    }
}
