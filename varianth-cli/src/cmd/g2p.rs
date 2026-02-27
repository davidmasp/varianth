

use g2p::g2p_run;

pub fn run(
    gff_path: String,
    genome_fasta_path: String,
    proteome_fasta_path: String,
    debug_flag: Option<usize>,
    output_path: String,
    thread_number: usize,
) {
    g2p_run(
        &gff_path,
        &genome_fasta_path,
        &proteome_fasta_path,
        debug_flag,
        &output_path,
        thread_number,
    );
}

