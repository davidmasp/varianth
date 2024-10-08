
mod core;
mod ms;
mod readinfo;
mod getrf;

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(author="DMP", version, about="A compilation of utilities for variant data", long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Adds mutation subtype (MS) into vcf file from a fasta
    Addms(AddmsArgs),
    Readinfo(ReadinfoArgs),
    Readfreq(ReadfreqArgs),
}

// to add a non-positional argument, we need to put this above the arg.
//#[clap(short, long)]

#[derive(Args)]
struct AddmsArgs {
    /// Integer value to define the number of subtype adjacent bases, use 1 for trinucleotide.
    #[clap(short, long, default_value = "1")]
    kval: u8,
    /// Fasta file with the reference genome. (needs to be indexed)
    #[clap(short, long)]
    genome: Option<PathBuf>,
    /// VCF file with to modify. (needs to be indexed)
    #[clap(short, long, conflicts_with = "use_stdin")]
    variants: Option<PathBuf>,
    /// A flag to indicate variants come from stdin instead of a file.
    #[clap(long, action)]
    use_stdin: bool,
    /// Output VCF file. Incompatible with --stdout/-O
    #[clap(short = 'o', long, default_value = "out.vcf.gz", conflicts_with = "use_stdout")]
    outfile: Option<PathBuf>,
    /// A flag to indicate stdout instead of a file.
    #[clap(long, action)]
    use_stdout: bool,
    /// Name for the information field to add into the vcf.
    #[clap(short, long, default_value = "MS")]
    infoname: Option<String>,
    /// Description for the information field to add into the vcf.
    #[clap(short = 'I', long, default_value = "mutation subtype")]
    infodescription: Option<String>,
}

#[derive(Args)]
struct ReadinfoArgs {
    /// BAM file with read information. (needs to be indexed)
    #[clap(short, long)]
    reads: Option<PathBuf>,
    /// VCF file with to modify. (needs to be indexed)
    #[clap(short, long)]
    variants: Option<PathBuf>,
    /// Output json file
    #[clap(short = 'o', long, default_value = "out.json")]
    outfile: Option<PathBuf>,
}

#[derive(Args)]
struct ReadfreqArgs {
    /// BAM file with read information. (needs to be indexed)
    #[clap(short, long)]
    reads: Option<PathBuf>,
    /// BED file with codons.
    #[clap(short, long)]
    variants: Option<PathBuf>,
    /// Output json file
    #[clap(short = 'o', long, default_value = "out.tsv")]
    outfile: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Addms(addmsargs) => {
            // println!("'myapp add' was used, name is: {:?}", addmsargs.infoname);
            let kval_in: usize = TryInto::try_into(addmsargs.kval).unwrap();
            let info_name = addmsargs.infoname.clone().unwrap();
            let info_description = addmsargs.infodescription.clone().unwrap();
            let use_stdin = addmsargs.use_stdin;
            let use_stdout = addmsargs.use_stdout;
            ms::addms(
                addmsargs.genome.clone().unwrap(),
                addmsargs.variants.clone().unwrap(),
                addmsargs.outfile.clone().unwrap(),
                kval_in,
                info_name,
                info_description,
                use_stdin,
                use_stdout,
            );
        },
        Commands::Readinfo(readinfoargs) => {

            let reads_file_result = readinfoargs.reads.clone();
            let reads_file = match reads_file_result {
                Some(reads_file) => reads_file,
                None => panic!("No reads file provided"),
            };

            let variants_file_result = readinfoargs.variants.clone();
            let variants_file = match variants_file_result {
                Some(variants_file) => variants_file,
                None => panic!("No variants file provided"),
            };

            let outfile_result = readinfoargs.outfile.clone();
            let outfile = match outfile_result {
                Some(outfile) => outfile,
                None => panic!("No outfile provided"),
            };

            readinfo::readinfo(
                reads_file,
                variants_file,
                outfile,
            );
        },

        Commands::Readfreq(readfreqargs) => {

            let reads_file_result = readfreqargs.reads.clone();
            let reads_file = match reads_file_result {
                Some(reads_file) => reads_file,
                None => panic!("No reads file provided"),
            };

            let variants_file_result = readfreqargs.variants.clone();
            let variants_file = match variants_file_result {
                Some(variants_file) => variants_file,
                None => panic!("No variants file provided"),
            };

            let outfile_result = readfreqargs.outfile.clone();
            let outfile = match outfile_result {
                Some(outfile) => outfile,
                None => panic!("No outfile provided"),
            };

            getrf::readfreq(
                reads_file,
                variants_file,
                outfile,
            );
        },
    }
}

