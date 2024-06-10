
mod cmd;

use cmd::ms;
use cmd::kmercount;

// LOGS
use simplelog;

use log::error;

// ARGUMENTS
use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Adds files to myapp
    Ms(MsArgs),
    Kcount(KcountArgs),
}

#[derive(Args)]
struct MsArgs {
    fasta: PathBuf,
}

#[derive(Args)]
struct KcountArgs {
    fasta: PathBuf,
    #[arg(short='K', long)]
    size: usize,
    #[arg(short='S', long)]
    table_size: Option<usize>,
    #[arg(short='r', long)]
    regions: Option<String>,
    #[arg(short='R', long)]
    regions_file: Option<PathBuf>,
    #[arg(short='o', long)]
    output: Option<PathBuf>,
    /// verbose flag
    #[arg(short='v', long)]
    verbose: bool,
}

fn main() {

    let _ = simplelog::SimpleLogger::init(simplelog::LevelFilter::Info, simplelog::Config::default());
    let cli = Cli::parse();

    match cli.command {
        Commands::Ms(args) => {
            ms::run();
        },
        Commands::Kcount(args) => {
            let rres = kmercount::run(args.fasta, args.size, args.regions, args.regions_file, args.output, args.table_size, args.verbose);
            match rres {
                Ok(_) => {},
                Err(e) => {
                    error!("Error: {}", e);
                }
            }
        }
    }
}
