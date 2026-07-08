mod cmd;

use cmd::gene2protein;
use cmd::kmercount;
use cmd::ms;
use cmd::readinfo;
use cmd::vep2table;

// LOGS
use simplelog;

use log::error;

// ARGUMENTS
use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Adds the mutation subtype (MS) to a VCF file based on the reference genome and variants provided.
    Ms(MsArgs),
    /// Counts k-mers from a FASTA file, optionally restricted to specific genomic regions.
    Kcount(KcountArgs),
    /// Expands VEP INFO/CSQ annotations from a VCF into a flat tabular output.
    Vep2table(Vep2tableArgs),
    /// Computes a histogram of read-start offsets for variants using an indexed BAM file.
    Readinfo(ReadinfoArgs),
    #[command(name = "gene2protein")]
    /// Generates all possible nonsynonymous SNVs by combining GFF, genome, and proteome inputs.
    Gene2protein(Gene2proteinArgs),
}

#[derive(Args)]
struct MsArgs {
    #[arg(short = 'g', long)]
    fasta: PathBuf,
    #[arg(short = 'i', long)]
    variants: PathBuf,
    #[arg(short = 'o', long)]
    output: PathBuf,
    #[arg(short = 'k', long)]
    kval: usize,
    #[arg(short = 'f', long, default_value = "MS")]
    feature: String,
    #[arg(short = 'F', long, default_value = "Mutation Subtype")]
    featuredescription: String,
}

#[derive(Args)]
struct KcountArgs {
    fasta: PathBuf,
    #[arg(short = 'K', long)]
    size: usize,
    #[arg(short = 'S', long)]
    table_size: Option<usize>,
    #[arg(short = 'r', long)]
    regions: Option<String>,
    #[arg(short = 'R', long)]
    regions_file: Option<PathBuf>,
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    /// verbose flag
    #[arg(short = 'v', long)]
    verbose: bool,
    /// Skip k-mers containing ambiguous bases (N or other IUPAC codes) instead of failing
    #[arg(long)]
    skip_ambiguous: bool,
}

#[derive(Args)]
struct Vep2tableArgs {
    #[arg(short = 'i', long)]
    input: PathBuf,
    #[arg(short = 'o', long)]
    output: PathBuf,
}

#[derive(Args)]
struct ReadinfoArgs {
    /// BAM file with read information. Requires a BAM index next to it.
    #[arg(short = 'r', long)]
    reads: PathBuf,
    /// VCF file with variants to query.
    #[arg(short = 'v', long)]
    variants: PathBuf,
    /// Output JSON histogram file.
    #[arg(short = 'o', long, default_value = "out.json")]
    output: PathBuf,
}

#[derive(Args)]
struct Gene2proteinArgs {
    #[arg(long, default_value = "MANE.GRCh38.v1.4.ensembl_genomic.gff.gz")]
    gff_path: String,
    #[arg(long, default_value = "genome.fa")]
    genome_fasta_path: String,
    #[arg(long, default_value = "MANE.GRCh38.v1.4.ensembl_protein.faa")]
    proteome_fasta_path: String,
    #[arg(long)]
    debug_flag: Option<usize>,
    #[arg(long, default_value = "tables/all_mutations")]
    output_prefix: String,
}

fn main() {
    let _ =
        simplelog::SimpleLogger::init(simplelog::LevelFilter::Info, simplelog::Config::default());
    let cli = Cli::parse();

    match cli.command {
        Commands::Ms(args) => {
            if let Err(e) = ms::run(
                args.fasta,
                args.variants,
                args.output,
                args.kval,
                args.feature,
                args.featuredescription,
            ) {
                error!("Error: {:#}", e);
                std::process::exit(1);
            }
        }
        Commands::Kcount(args) => {
            if let Err(e) = kmercount::run(
                args.fasta,
                args.size,
                args.regions,
                args.regions_file,
                args.output,
                args.table_size,
                args.verbose,
                args.skip_ambiguous,
            ) {
                error!("Error: {:#}", e);
                std::process::exit(1);
            }
        }
        Commands::Vep2table(args) => {
            vep2table::run(args.input, args.output);
        }
        Commands::Readinfo(args) => {
            if let Err(e) = readinfo::run(args.reads, args.variants, args.output) {
                error!("Error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Gene2protein(args) => {
            if let Err(e) = gene2protein::run(
                &args.gff_path,
                &args.genome_fasta_path,
                &args.proteome_fasta_path,
                args.debug_flag,
                &args.output_prefix,
            ) {
                error!("Error: {:#}", e);
                std::process::exit(1);
            }
        }
    }
}
