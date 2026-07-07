use context::addms;
use std::path::PathBuf;

pub fn run(
    genome: PathBuf,
    variants_in: PathBuf,
    variants_out: PathBuf,
    kval: usize,
    key_name: String,
    key_description: String,
) {
    addms(
        genome,
        variants_in,
        variants_out,
        kval,
        key_name,
        key_description,
        false,
        false,
    );
}
