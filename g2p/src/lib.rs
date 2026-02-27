

use std::{collections::HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::mpsc;

pub mod gff;
pub mod fasta;
pub mod codons;
pub mod flow;

pub use codons::{CodonError, NonSynonymousMutation, expand_codons_from_sequence, MutationList};
pub use gff::{collect_cds_by_protein_id, GffReader, GffRecord, Strand, CdsProteome};
pub use fasta::{open_indexed_fasta, pull_entire_record, sequence_count, sequence_ids, reverse_complement};

pub use flow::g2pflow;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

pub fn g2p_run(
    gff_path: &str,
    genome_fasta_path: &str,
    proteome_fasta_path: &str,
    debug_flag: Option<usize>,
    output_path: &str,
    thread_number: usize,
) {
    let output_path = output_path.to_string();
    let output_dir = Path::new(&output_path).parent();
    if let Some(dir) = output_dir {
        if !dir.is_dir() {
            log::error!("Output directory does not exist: {}", dir.display());
            std::process::exit(1);
        }
    }

    ThreadPoolBuilder::new()
        .num_threads(thread_number)
        .build_global()
        .expect("failed to build global thread pool");

    // fai derived automatically from genome_fasta_path + ".fai"
    let fasta_reader_proteins = open_indexed_fasta(proteome_fasta_path, None::<&str>);
    let _protein_index = fasta_reader_proteins.index().clone();
    let seq_count = sequence_count(&fasta_reader_proteins);
    log::info!("Number of sequences in proteome FASTA: {}", seq_count);

    let prot_ids: HashSet<String> = sequence_ids(&fasta_reader_proteins).into_iter().collect();
    let proteome = collect_cds_by_protein_id(gff_path, &prot_ids)
        .expect("failed to collect CDS records by protein_id");
    log::info!("Proteome length: {}", proteome.len());

    // this is just a sanity check to make sure all proteins are matched
    let proteome_keys = proteome.protein_ids();
    let missing_proteins: Vec<&String> = prot_ids.difference(&proteome_keys).collect();
    if !missing_proteins.is_empty() {
        log::warn!("The following protein IDs were found in the proteome FASTA but are missing from the GFF CDS records: {:?}", missing_proteins);
    } else {
        log::info!("All protein IDs from the proteome FASTA are present in the GFF CDS records.");
    }

    // this is to maintain a small footprint even when using the big files
    let proteome_keys_input = match debug_flag {
        Some(limit) => {
            log::warn!("Debug counter limit set.");
            proteome_keys.into_iter().take(limit).collect::<HashSet<String>>()
        }
        _ => {
            proteome_keys
        }
    };

    let (tx, rx) = mpsc::channel::<Vec<String>>();
    let writer_handle = std::thread::spawn(move || -> std::io::Result<()> {
        let out_file = File::create(output_path)?;
        let mut writer = BufWriter::new(out_file);
        /*
        this is a bit of a mess for soprting the file later so maybe
        better we just get rid of the headers.
        writeln!(
            writer,
            "chr\tgenome_position\tref_dna\talt_dna\tprotein_id\tref_aa\tprotein_position\talt_aa\tincodon_position"
        )?;
        */
        for lines in rx {
            for line in lines {
                writeln!(writer, "{}", line)?;
            }
        }
        writer.flush()?;
        Ok(())
    });

    proteome_keys_input.par_iter().for_each_with(tx, |tx, pid| {
        let cds_vec = proteome.get_cloned(pid).expect("Error in internal GFF object.");
        log::debug!("{}: {} CDS records", pid, cds_vec.len());
        let pid_mutation_list_result =
            g2pflow(pid, cds_vec, genome_fasta_path, proteome_fasta_path);

        let pid_mutation_list = match pid_mutation_list_result {
            Ok(ml) => ml,
            Err(e) => {
                log::error!("Error processing protein_id {}: {}", pid, e);
                panic!("Error processing protein_id {}: {}", pid, e);
            }
        };

        let lines = pid_mutation_list.format_lines();
        tx.send(lines).expect("failed to send TSV lines");
    });

    match writer_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => panic!("failed to write output file: {}", e),
        Err(_) => panic!("writer thread panicked"),
    }
}



