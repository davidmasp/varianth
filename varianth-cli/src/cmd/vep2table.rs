

use vep2table::io::vep2table;
use std::path::PathBuf;

pub fn run (input_file: PathBuf, output_file: PathBuf) {
    vep2table(&input_file.to_string_lossy(), &output_file.to_string_lossy());
}
