
// THIS CODE IS BASED ORIGINALLY IN THE SEQ2MUT MODULE

use crate::gff::Strand;

use std::collections::HashMap;
use bstr::{BString, ByteSlice};

// this is somehow trying to mimic Biostrings::GENETIC_CODE from R
struct HumanGeneticCode {
    codon_map: HashMap<BString, BString>,
    alt_init_codons: Vec<BString>,
}

impl HumanGeneticCode {
    fn new() -> Self {
        let codon_map = generate_codon_map();
        // ALT start codons are only used if they are in the first position of the CDS, and they code for Methionine (M) instead of their usual amino acid.
        // see https://www.tandfonline.com/doi/full/10.4161/trla.28387#d1e583
        let alt_init_codons = vec![BString::from("TTG"), BString::from("CTG"),  BString::from("GTG")];
        HumanGeneticCode {
            codon_map,
            alt_init_codons,
        }
    }
    fn translate_codon(&self, codon: &BString, prot_position: &usize) -> Option<u8> {
        if *prot_position == 1 {
            if self.alt_init_codons.contains(codon) {
                Some(b'M')
            } else {
                self.codon_map.get(codon).map(|aa| aa[0])
            }
        } else {
            self.codon_map.get(codon).map(|aa| aa[0])
        }
    }
}

fn generate_codon_map() -> HashMap<BString, BString> {
    let mut codon_map = HashMap::new();
    codon_map.insert(BString::from("AAA"), BString::from("K"));
    codon_map.insert(BString::from("AAC"), BString::from("N"));
    codon_map.insert(BString::from("AAG"), BString::from("K"));
    codon_map.insert(BString::from("AAT"), BString::from("N"));
    codon_map.insert(BString::from("ACA"), BString::from("T"));
    codon_map.insert(BString::from("ACC"), BString::from("T"));
    codon_map.insert(BString::from("ACG"), BString::from("T"));
    codon_map.insert(BString::from("ACT"), BString::from("T"));
    codon_map.insert(BString::from("AGA"), BString::from("R"));
    codon_map.insert(BString::from("AGC"), BString::from("S"));
    codon_map.insert(BString::from("AGG"), BString::from("R"));
    codon_map.insert(BString::from("AGT"), BString::from("S"));
    codon_map.insert(BString::from("ATA"), BString::from("I"));
    codon_map.insert(BString::from("ATC"), BString::from("I"));
    codon_map.insert(BString::from("ATG"), BString::from("M"));
    codon_map.insert(BString::from("ATT"), BString::from("I"));
    codon_map.insert(BString::from("CAA"), BString::from("Q"));
    codon_map.insert(BString::from("CAC"), BString::from("H"));
    codon_map.insert(BString::from("CAG"), BString::from("Q"));
    codon_map.insert(BString::from("CAT"), BString::from("H"));
    codon_map.insert(BString::from("CCA"), BString::from("P"));
    codon_map.insert(BString::from("CCC"), BString::from("P"));
    codon_map.insert(BString::from("CCG"), BString::from("P"));
    codon_map.insert(BString::from("CCT"), BString::from("P"));
    codon_map.insert(BString::from("CGA"), BString::from("R"));
    codon_map.insert(BString::from("CGC"), BString::from("R"));
    codon_map.insert(BString::from("CGG"), BString::from("R"));
    codon_map.insert(BString::from("CGT"), BString::from("R"));
    codon_map.insert(BString::from("CTA"), BString::from("L"));
    codon_map.insert(BString::from("CTC"), BString::from("L"));
    codon_map.insert(BString::from("CTG"), BString::from("L"));
    codon_map.insert(BString::from("CTT"), BString::from("L"));
    codon_map.insert(BString::from("GAA"), BString::from("E"));
    codon_map.insert(BString::from("GAC"), BString::from("D"));
    codon_map.insert(BString::from("GAG"), BString::from("E"));
    codon_map.insert(BString::from("GAT"), BString::from("D"));
    codon_map.insert(BString::from("GCA"), BString::from("A"));
    codon_map.insert(BString::from("GCC"), BString::from("A"));
    codon_map.insert(BString::from("GCG"), BString::from("A"));
    codon_map.insert(BString::from("GCT"), BString::from("A"));
    codon_map.insert(BString::from("GGA"), BString::from("G"));
    codon_map.insert(BString::from("GGC"), BString::from("G"));
    codon_map.insert(BString::from("GGG"), BString::from("G"));
    codon_map.insert(BString::from("GGT"), BString::from("G"));
    codon_map.insert(BString::from("GTA"), BString::from("V"));
    codon_map.insert(BString::from("GTC"), BString::from("V"));
    codon_map.insert(BString::from("GTG"), BString::from("V"));
    codon_map.insert(BString::from("GTT"), BString::from("V"));
    codon_map.insert(BString::from("TAA"), BString::from("*"));
    codon_map.insert(BString::from("TAC"), BString::from("Y"));
    codon_map.insert(BString::from("TAG"), BString::from("*"));
    codon_map.insert(BString::from("TAT"), BString::from("Y"));
    codon_map.insert(BString::from("TCA"), BString::from("S"));
    codon_map.insert(BString::from("TCC"), BString::from("S"));
    codon_map.insert(BString::from("TCG"), BString::from("S"));
    codon_map.insert(BString::from("TCT"), BString::from("S"));
    codon_map.insert(BString::from("TGA"), BString::from("*"));
    codon_map.insert(BString::from("TGC"), BString::from("C"));
    codon_map.insert(BString::from("TGG"), BString::from("W"));
    codon_map.insert(BString::from("TGT"), BString::from("C"));
    codon_map.insert(BString::from("TTA"), BString::from("L"));
    codon_map.insert(BString::from("TTC"), BString::from("F"));
    codon_map.insert(BString::from("TTG"), BString::from("L"));
    codon_map.insert(BString::from("TTT"), BString::from("F"));
    codon_map
}


pub struct MutationList {
    pub mutations: Vec<NonSynonymousMutation>,
    pub protein_id: String,
    pub sequence_id: BString, // this is the chromosome or contig name, e.g. "chr1"
    pub strand: Strand,
}

impl MutationList {
    pub fn new(mutations: Vec<NonSynonymousMutation>, protein_id: String, sequence_id: BString, strand: Strand) -> Self {
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
}

impl std::fmt::Display for CodonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CodonError::DnaProteinLengthMismatch { dna_len, protein_len } => write!(
                f,
                "DNA length/3 does not match protein length + 1 (dna_len={}, protein_len={})",
                dna_len,
                protein_len
            ),
            CodonError::ReferenceAminoAcidMismatch {
                position,
                expected,
                found,
            } => write!(
                f,
                "reference amino acid mismatch at position {}: expected {}, got {}",
                position,
                expected,
                found
            ),
            CodonError::InconsistentStrand { protein_id } => write!(
                f,
                "Inconsistent strand information for protein_id {}",
                protein_id
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
    genome_pos: Vec<usize>,
    dna_seq: BString,
    protein_seq: BString,
) -> Result<Vec<NonSynonymousMutation>, CodonError> {

    let mut generated_mutations: Vec<NonSynonymousMutation> = Vec::new();
    // here the protein_position needs to be /3 size the seq length.
    let genetic_code = HumanGeneticCode::new();
    let valid_dna = vec![b'A', b'C', b'G', b'T'];

    if dna_seq.len() / 3 != protein_seq.len() + 1 {
        return Err(CodonError::DnaProteinLengthMismatch {
            dna_len: dna_seq.len(),
            protein_len: protein_seq.len(),
        });
    }

    let dna_chunks = dna_seq.chunks(3).take(protein_seq.len());
    let pos_chunks = genome_pos.chunks(3).take(protein_seq.len());

    let input_triple_zip: Vec<(&[u8], char, &[usize])> = dna_chunks
        .zip(protein_seq.chars())
        .zip(pos_chunks)
        .map(|((codon, aa), pos)| (codon, aa, pos))
        .collect::<Vec<_>>();
    
    // this .take(protein_seq.len() 
    // takes all codons but the last one, which is the stop codon.
    for (pos_idx, (dna_codon, aa, position_codon)) in input_triple_zip.iter().enumerate() {
            let aa = *aa as u8;
            let cdn_bstr = BString::from(*dna_codon);
            let cdn_pos: Vec<usize> = position_codon.to_vec();
            let pos = pos_idx + 1; // 1-based protein coordinate
            let wt_cdn_aa =  genetic_code.translate_codon(&cdn_bstr, &pos).unwrap();
            
            if wt_cdn_aa != aa {
                return Err(CodonError::ReferenceAminoAcidMismatch {
                    position: pos,
                    expected: aa as char,
                    found: wt_cdn_aa as char,
                });
            }
            
            for i in 0..3 {
                for k in &valid_dna {
                    if dna_codon[i] != *k {
                        let mut cdn_to_modify = dna_codon.to_vec();
                        cdn_to_modify[i] = *k;
                        let cdn_mut_bstr = BString::from(cdn_to_modify);
                        let new_aa = genetic_code
                            .translate_codon(&cdn_mut_bstr, &pos)
                            .unwrap();
                        if new_aa != aa {
                            // this includes both the missense and nonsense mutations, but not the synonymous ones.
                            
                            // pos is the position in the protein
                            // i is the position in the codon (0, 1, or 2), aka phase
                            // k is the ALT new nucleotide
                            // new_aa is the new amino acid after the mutation
                            // wt_cdn_aa is the original aa

                            // in codon position of phase is 0 based
                            let mut_ns = NonSynonymousMutation::new(
                                pos,
                                i,
                                cdn_pos[i],
                                dna_codon[i] as char,
                                *k as char,
                                wt_cdn_aa as char,
                                new_aa as char,
                            );
                            generated_mutations.push(mut_ns);
                        }
                    } else {
                        // this does nothing
                    }
                }
            } // end of the loop per codon
        }
        Ok(generated_mutations)
    }

