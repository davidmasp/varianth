pub use crate::codons::{
    CodonError, MutationList, NonSynonymousMutation, expand_codons_from_sequence,
};
pub use crate::fasta::{pull_entire_record, reverse_complement};
pub use crate::gff::{GffReader, GffRecord, Strand, collect_cds_by_protein_id};

use anyhow::{Context, Result, bail, ensure};
use bstr::BString;
use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::fasta::fai;
use std::fs::File;
use std::io::BufReader;

pub fn g2pflow(
    pid: &str,
    mut cds_vec: Vec<GffRecord>,
    genome_reader: &mut fasta::io::IndexedReader<BufReader<File>>,
    proteome_reader: &mut fasta::io::IndexedReader<BufReader<File>>,
    proteome_index: &fai::Index,
) -> Result<MutationList> {
    // I am pretty sure this should be sorted already.
    cds_vec.sort_by_key(|cds| cds.start);

    // here we need a first iteration to pull out the strand
    let first_cds = cds_vec
        .first()
        .with_context(|| format!("no CDS records found for protein_id {pid}"))?;
    let unique_strand = first_cds.strand.clone();
    if cds_vec.iter().any(|cds| cds.strand != unique_strand) {
        return Err(CodonError::InconsistentStrand {
            protein_id: pid.to_string(),
        }
        .into());
    }
    log::debug!("{}: unique strand: {:?}", pid, unique_strand);
    // here we need to pull out both dna and genomic positions for each CDS
    let cds_extracted_info: Vec<(Vec<u8>, Vec<usize>)> = cds_vec
        .iter()
        .map(|cds| {
            let start_pos = Position::try_from(cds.start).with_context(|| {
                format!("failed to convert CDS start position {} for protein_id {pid}", cds.start)
            })?;
            let end_pos = Position::try_from(cds.end).with_context(|| {
                format!("failed to convert CDS end position {} for protein_id {pid}", cds.end)
            })?;
            let cds_region = Region::new(cds.seqid.clone(), start_pos..=end_pos);
            let cds_sequence = genome_reader.query(&cds_region).map_err(|_| {
                CodonError::MissingReferenceSequence {
                    protein_id: pid.to_string(),
                    seqid: cds.seqid.to_string(),
                }
            })?;
            let dna_seq = cds_sequence.sequence().as_ref().to_vec();
            let gpos = (cds.start..=cds.end).collect::<Vec<usize>>();
            ensure!(
                dna_seq.len() == gpos.len(),
                "length of extracted DNA sequence ({}) does not match genomic position range ({}) for protein_id {pid}",
                dna_seq.len(),
                gpos.len()
            );
            Ok((dna_seq, gpos))
        })
        .collect::<Result<Vec<_>>>()?;

    // once we have the full sequence we can do the reverse complement if needed.
    let full_cds_sequence_vec = cds_extracted_info
        .iter()
        .flat_map(|(seq, _)| seq.clone())
        .collect::<Vec<u8>>();
    let mut full_cds_sequence = BString::from(full_cds_sequence_vec);
    full_cds_sequence.make_ascii_uppercase();
    if unique_strand == Strand::Reverse {
        full_cds_sequence = reverse_complement(&full_cds_sequence).with_context(|| {
            format!("failed to compute reverse complement for protein_id {pid}")
        })?;
        log::debug!(
            "{}: reverse complemented CDS sequence:\n{}\n",
            pid,
            full_cds_sequence
        );
    }

    // in case we need to do the reverse complement, we also need to reverse the genomic positions.
    let mut full_gpos = cds_extracted_info
        .iter()
        .flat_map(|(_, gpos)| gpos.clone())
        .collect::<Vec<usize>>();
    match unique_strand {
        Strand::Reverse => full_gpos.reverse(),
        Strand::Forward => {}
        Strand::Unknown => bail!("missing strand information for protein_id {pid}"),
    }

    // and we here need to pull the third element of the main function, the protein sequence
    let prot_seq = pull_entire_record(proteome_reader, proteome_index, pid).with_context(|| {
        format!("failed to pull protein sequence from FASTA for protein_id {pid}")
    })?;

    let mutation_list_raw: Vec<NonSynonymousMutation> =
        expand_codons_from_sequence(&full_gpos, full_cds_sequence.as_ref(), &prot_seq)
            .with_context(|| format!("failed to expand codons for protein_id {pid}"))?;

    let mutation_list: MutationList = MutationList::new(
        mutation_list_raw,
        pid.to_string(),
        first_cds.seqid.clone(),
        unique_strand,
    );
    Ok(mutation_list)
}
