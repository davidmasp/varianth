

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use indicatif::{ProgressBar, ProgressStyle};
use serde::Serialize;

pub mod gff;
pub mod fasta;
pub mod codons;
pub mod flow;

pub use codons::{CodonError, NonSynonymousMutation, expand_codons_from_sequence, MutationList};
pub use gff::{collect_cds_by_protein_id, GffReader, GffRecord, Strand, CdsProteome};
pub use fasta::{open_indexed_fasta, pull_entire_record, sequence_count, sequence_ids, reverse_complement};

pub use flow::g2pflow;

#[derive(Serialize)]
pub struct G2pMetrics {
    pub total_proteins: usize,
    pub successful_count: usize,
    pub failed_count: usize,
    pub failed_ids: Vec<String>,
    pub failed_errors: HashMap<String, String>,
    pub elapsed_seconds: f64,
    pub proteins_per_second: f64,
}

pub fn g2p_run(
    gff_path: &str,
    genome_fasta_path: &str,
    proteome_fasta_path: &str,
    debug_flag: Option<usize>,
    output_prefix: &str,
) {
    let tsv_path = format!("{}.tsv", output_prefix);
    let json_path = format!("{}.json", output_prefix);

    let output_dir = Path::new(output_prefix).parent();
    if let Some(dir) = output_dir {
        if !dir.as_os_str().is_empty() && !dir.is_dir() {
            log::error!("Output directory does not exist: {}", dir.display());
            std::process::exit(1);
        }
    }

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
    let proteome_keys_input: Vec<String> = match debug_flag {
        Some(limit) => {
            log::warn!("Debug counter limit set.");
            proteome_keys.into_iter().take(limit).collect()
        }
        _ => {
            proteome_keys.into_iter().collect()
        }
    };

    let total_proteins = proteome_keys_input.len();
    let mut successful_count: usize = 0;
    let mut failed_ids: Vec<String> = Vec::new();
    let mut failed_errors: HashMap<String, String> = HashMap::new();

    let out_file = File::create(&tsv_path).expect("failed to create output TSV file");
    let mut writer = BufWriter::new(out_file);

    let start = Instant::now();

    let pb = ProgressBar::new(total_proteins as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec}, ETA {eta})"
        )
        .expect("failed to parse progress bar template")
        .progress_chars("#>-"),
    );

    for pid in &proteome_keys_input {
        let cds_vec = proteome.get_cloned(pid).expect("Error in internal GFF object.");
        log::debug!("{}: {} CDS records", pid, cds_vec.len());

        let pid_mutation_list_result =
            g2pflow(pid, cds_vec, genome_fasta_path, proteome_fasta_path);

        match pid_mutation_list_result {
            Ok(ml) => {
                let lines = ml.format_lines();
                for line in lines {
                    writeln!(writer, "{}", line).expect("failed to write TSV line");
                }
                successful_count += 1;
            }
            Err(e) => {
                log::error!("Error processing protein_id {}: {}", pid, e);
                failed_ids.push(pid.clone());
                failed_errors.insert(pid.clone(), e.to_string());
            }
        }
        pb.inc(1);
    }
    pb.finish_with_message("done");

    writer.flush().expect("failed to flush output TSV file");

    let elapsed = start.elapsed().as_secs_f64();
    let proteins_per_second = if elapsed > 0.0 {
        total_proteins as f64 / elapsed
    } else {
        0.0
    };

    let metrics = G2pMetrics {
        total_proteins,
        successful_count,
        failed_count: failed_ids.len(),
        failed_ids,
        failed_errors,
        elapsed_seconds: elapsed,
        proteins_per_second,
    };

    let json = serde_json::to_string_pretty(&metrics).expect("failed to serialize metrics");
    let mut json_file = File::create(&json_path).expect("failed to create metrics JSON file");
    json_file.write_all(json.as_bytes()).expect("failed to write metrics JSON file");

    log::info!(
        "Done. {} proteins processed ({} succeeded, {} failed) in {:.2}s. Metrics written to {}",
        total_proteins,
        metrics.successful_count,
        metrics.failed_count,
        elapsed,
        json_path,
    );
}



