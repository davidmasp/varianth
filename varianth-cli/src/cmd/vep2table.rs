use std::path::PathBuf;
use vep2table::io::vep2table;

pub fn run(input_file: PathBuf, output_file: PathBuf) {
    vep2table(
        &input_file.to_string_lossy(),
        &output_file.to_string_lossy(),
    );
}
