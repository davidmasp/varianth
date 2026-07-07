// THIS CODE IS BASED ORIGINALLY IN THE SEQ2MUT MODULE

use crate::gff::Strand;

use bstr::{BString, ByteSlice};

const VALID_DNA: [u8; 4] = [b'A', b'C', b'G', b'T'];

// Dense 64-entry standard genetic code lookup with A/C/G/T => 0/1/2/3 encoding.
const CODON_TABLE: [u8; 64] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S', b'I', b'I', b'M', b'I',
    b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P', b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L',
    b'E', b'D', b'E', b'D', b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C', b'L', b'F', b'L', b'F',
];

#[inline]
fn base_to_bits(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn codon_to_index(c0: u8, c1: u8, c2: u8) -> Option<usize> {
    let b0 = base_to_bits(c0)?;
    let b1 = base_to_bits(c1)?;
    let b2 = base_to_bits(c2)?;
    Some((b0 << 4) | (b1 << 2) | b2)
}

#[inline]
fn is_alt_start_codon(c0: u8, c1: u8, c2: u8) -> bool {
    ((c0 == b'T' || c0 == b't') && (c1 == b'T' || c1 == b't') && (c2 == b'G' || c2 == b'g'))
        || ((c0 == b'C' || c0 == b'c') && (c1 == b'T' || c1 == b't') && (c2 == b'G' || c2 == b'g'))
        || ((c0 == b'G' || c0 == b'g') && (c1 == b'T' || c1 == b't') && (c2 == b'G' || c2 == b'g'))
}

#[inline]
fn translate_codon(c0: u8, c1: u8, c2: u8, prot_position: usize) -> Option<u8> {
    if prot_position == 1 && is_alt_start_codon(c0, c1, c2) {
        return Some(b'M');
    }
    codon_to_index(c0, c1, c2).map(|idx| CODON_TABLE[idx])
}

pub struct MutationList {
    pub mutations: Vec<NonSynonymousMutation>,
    pub protein_id: String,
    pub sequence_id: BString, // this is the chromosome or contig name, e.g. "chr1"
    pub strand: Strand,
}

impl MutationList {
    pub fn new(
        mutations: Vec<NonSynonymousMutation>,
        protein_id: String,
        sequence_id: BString,
        strand: Strand,
    ) -> Self {
        Self {
            mutations,
            protein_id,
            sequence_id,
            strand,
        }
    }

    pub fn format_lines(&self) -> Vec<String> {
        let chr = self.sequence_id.to_str_lossy();

        self.mutations
            .iter()
            .map(|m| {
                let (ref_nucleotide, alt_nucleotide) = match self.strand {
                    Strand::Reverse => (
                        complement_nucleotide(m.ref_nucleotide),
                        complement_nucleotide(m.alt_nucleotide),
                    ),
                    _ => (m.ref_nucleotide, m.alt_nucleotide),
                };

                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    chr,
                    m.genome_position,
                    ref_nucleotide,
                    alt_nucleotide,
                    self.protein_id,
                    m.ref_aa,
                    m.prot_position,
                    m.alt_aa,
                    m.incodon_position,
                )
            })
            .collect()
    }
}

fn complement_nucleotide(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'a' => 't',
        't' => 'a',
        'c' => 'g',
        'g' => 'c',
        _ => base,
    }
}

pub struct NonSynonymousMutation {
    pub prot_position: usize,
    pub incodon_position: usize,
    pub genome_position: usize,
    pub ref_nucleotide: char,
    pub alt_nucleotide: char,
    pub ref_aa: char,
    pub alt_aa: char,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CodonError {
    DnaProteinLengthMismatch {
        dna_len: usize,
        protein_len: usize,
    },
    ReferenceAminoAcidMismatch {
        position: usize,
        expected: char,
        found: char,
    },
    InconsistentStrand {
        protein_id: String,
    },
    MissingReferenceSequence {
        protein_id: String,
        seqid: String,
    },
}

impl std::fmt::Display for CodonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CodonError::DnaProteinLengthMismatch {
                dna_len,
                protein_len,
            } => write!(
                f,
                "DNA length/3 does not match protein length + 1 (dna_len={}, protein_len={})",
                dna_len, protein_len
            ),
            CodonError::ReferenceAminoAcidMismatch {
                position,
                expected,
                found,
            } => write!(
                f,
                "reference amino acid mismatch at position {}: expected {}, got {}",
                position, expected, found
            ),
            CodonError::InconsistentStrand { protein_id } => write!(
                f,
                "Inconsistent strand information for protein_id {}",
                protein_id
            ),
            CodonError::MissingReferenceSequence { protein_id, seqid } => write!(
                f,
                "Reference sequence '{}' not found in genome FASTA for protein_id {}",
                seqid, protein_id
            ),
        }
    }
}

impl std::error::Error for CodonError {}

impl NonSynonymousMutation {
    fn new(
        prot_position: usize,
        incodon_position: usize,
        genome_position: usize,
        ref_nucleotide: char,
        alt_nucleotide: char,
        ref_aa: char,
        alt_aa: char,
    ) -> Self {
        Self {
            prot_position,
            incodon_position,
            genome_position,
            ref_nucleotide,
            alt_nucleotide,
            ref_aa,
            alt_aa,
        }
    }
}

pub fn expand_codons_from_sequence(
    genome_pos: &[usize],
    dna_seq: &[u8],
    protein_seq: &[u8],
) -> Result<Vec<NonSynonymousMutation>, CodonError> {
    let mut generated_mutations: Vec<NonSynonymousMutation> =
        Vec::with_capacity(protein_seq.len().saturating_mul(6));

    if dna_seq.len() / 3 != protein_seq.len() + 1 {
        return Err(CodonError::DnaProteinLengthMismatch {
            dna_len: dna_seq.len(),
            protein_len: protein_seq.len(),
        });
    }

    // This consumes all codons except the trailing stop codon.
    for codon_idx in 0..protein_seq.len() {
        let prot_position = codon_idx + 1; // 1-based protein coordinate
        let offset = codon_idx * 3;

        let c0 = dna_seq[offset];
        let c1 = dna_seq[offset + 1];
        let c2 = dna_seq[offset + 2];
        let aa = protein_seq[codon_idx];

        let wt_cdn_aa =
            translate_codon(c0, c1, c2, prot_position).expect("invalid codon in DNA sequence");
        if wt_cdn_aa != aa {
            return Err(CodonError::ReferenceAminoAcidMismatch {
                position: prot_position,
                expected: aa as char,
                found: wt_cdn_aa as char,
            });
        }

        let codon_positions = [
            genome_pos[offset],
            genome_pos[offset + 1],
            genome_pos[offset + 2],
        ];
        let ref_codon = [c0, c1, c2];

        for i in 0..3 {
            for alt_base in VALID_DNA {
                if ref_codon[i] == alt_base {
                    continue;
                }

                let mut mutated = ref_codon;
                mutated[i] = alt_base;
                let new_aa = translate_codon(mutated[0], mutated[1], mutated[2], prot_position)
                    .expect("invalid codon generated during mutation expansion");
                if new_aa != aa {
                    generated_mutations.push(NonSynonymousMutation::new(
                        prot_position,
                        i,
                        codon_positions[i],
                        ref_codon[i] as char,
                        alt_base as char,
                        wt_cdn_aa as char,
                        new_aa as char,
                    ));
                }
            }
        }
    }
    Ok(generated_mutations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translate_standard_codon_works() {
        assert_eq!(translate_codon(b'A', b'T', b'G', 4), Some(b'M'));
        assert_eq!(translate_codon(b'T', b'A', b'A', 4), Some(b'*'));
    }

    #[test]
    fn translate_alt_start_only_at_position_one() {
        assert_eq!(translate_codon(b'T', b'T', b'G', 1), Some(b'M'));
        assert_eq!(translate_codon(b'T', b'T', b'G', 2), Some(b'L'));
    }

    #[test]
    fn expand_returns_reference_mismatch() {
        let genome_pos = [1, 2, 3, 4, 5, 6];
        let dna_seq = b"ATGTAA";
        let protein_seq = b"A";

        match expand_codons_from_sequence(&genome_pos, dna_seq, protein_seq) {
            Err(CodonError::ReferenceAminoAcidMismatch {
                position,
                expected,
                found,
            }) => {
                assert_eq!(position, 1);
                assert_eq!(expected, 'A');
                assert_eq!(found, 'M');
            }
            Ok(_) => panic!("expected reference amino acid mismatch"),
            Err(_) => panic!("unexpected error variant"),
        }
    }

    #[test]
    fn expand_emits_nonsynonymous_mutations() {
        let genome_pos = [1, 2, 3, 4, 5, 6];
        let dna_seq = b"ATGTAA";
        let protein_seq = b"M";

        let mutations = expand_codons_from_sequence(&genome_pos, dna_seq, protein_seq).unwrap();
        assert!(!mutations.is_empty());
    }
}
