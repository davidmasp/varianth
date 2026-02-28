


pub use crate::codons::{CodonError, NonSynonymousMutation, expand_codons_from_sequence, MutationList};
pub use crate::gff::{collect_cds_by_protein_id, GffReader, GffRecord, Strand};
pub use crate::fasta::{open_indexed_fasta, pull_entire_record, sequence_count, sequence_ids, reverse_complement};

use noodles::core::{Region, Position};
use bstr::BString;

pub fn g2pflow(
    pid: &String,
    mut cds_vec: Vec<GffRecord>,
    genome_fasta_path: &str,
    proteome_fasta_path: &str,
) -> Result<MutationList, CodonError> {

    // here i initialize the fasta reader
    let mut fasta_reader_genome = open_indexed_fasta(genome_fasta_path, None::<&str>);

    // I am pretty sure this should be sorted already.
    cds_vec.sort_by_key(|cds| cds.start);

    // here we need a first iteration to pull out the strand
    let unique_strand = cds_vec.first().map(|cds| cds.strand.clone());
    if let Some(ref strand) = unique_strand {
        if cds_vec.iter().any(|cds| cds.strand != *strand) {
            return Err(CodonError::InconsistentStrand {
                protein_id: pid.to_string(),
            });
        }
    }
    log::debug!("{}: unique strand: {:?}", pid, unique_strand);
    // here we need to pull out both dna and genomic positions for each CDS
    let cds_extracted_info: Vec<(Vec<u8>, Vec<usize>)> = cds_vec
        .iter()
        .map(|cds| {
            let start_pos =
                Position::try_from(cds.start).expect("failed to convert start position");
            let end_pos = Position::try_from(cds.end).expect("failed to convert end position");
            let cds_region = Region::new(cds.seqid.clone(), start_pos..=end_pos);
            let cds_sequence = fasta_reader_genome
                .query(&cds_region)
                .map_err(|_| CodonError::MissingReferenceSequence {
                    protein_id: pid.to_string(),
                    seqid: cds.seqid.to_string(),
                })?;
            let dna_seq = cds_sequence.sequence().as_ref().to_vec();
            let gpos = (cds.start..=cds.end)
                .collect::<Vec<usize>>();
            assert_eq!(
                dna_seq.len(),
                gpos.len(),
                "length of extracted DNA sequence does not match genomic position range"
            );
            Ok((dna_seq, gpos))
        })
        .collect::<Result<Vec<_>, CodonError>>()?;

    // once we have the full sequence we can do the reverse complement if needed.
    let full_cds_sequence_vec = cds_extracted_info.iter()
        .flat_map(|(seq, _)| seq.clone())
        .collect::<Vec<u8>>();
    let mut full_cds_sequence = BString::from(full_cds_sequence_vec);
    full_cds_sequence.make_ascii_uppercase();
    if let Some(ref strand) = unique_strand {
        if strand == &Strand::Reverse {
            full_cds_sequence = reverse_complement(&full_cds_sequence).expect("failed to compute reverse complement");
            log::debug!("{}: reverse complemented CDS sequence:\n{}\n", pid, full_cds_sequence);
        }
    } else {
        // id think this is not possible.
        panic!("Missing strand information for protein_id {}", pid);
    }

    // in case we need to do the reverse complement, we also need to reverse the genomic positions.
    let mut full_gpos = cds_extracted_info.iter()
        .flat_map(|(_, gpos)| gpos.clone())
        .collect::<Vec<usize>>();
    match unique_strand {
        Some(Strand::Reverse) => full_gpos.reverse(),
        Some(Strand::Forward) => {}
        Some(Strand::Unknown) => panic!("Missing strand information for protein_id {}", pid),
        None => panic!("Missing strand information for protein_id {}", pid),
    }

    // the reason we do this here again is because we can paralellize faster
    let mut fasta_reader_proteins_instance = open_indexed_fasta(proteome_fasta_path, None::<&str>);
    let protein_index_instance = fasta_reader_proteins_instance.index().clone();

    // and we here need to pull the third element of the main function, the protein sequence
    let prot_seq = pull_entire_record(&mut fasta_reader_proteins_instance, &protein_index_instance, pid)
        .expect("failed to pull protein sequence from FASTA");

    let mutation_list_raw: Vec<NonSynonymousMutation> =
        expand_codons_from_sequence(full_gpos, full_cds_sequence, BString::from(prot_seq))?;

    let mutation_list: MutationList = MutationList::new(
        mutation_list_raw,
        pid.to_string(),
        cds_vec.first().unwrap().seqid.clone(),
        unique_strand.expect("Missing strand information for protein_id"),
    );
    Ok(mutation_list)

}

