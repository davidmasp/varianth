use std::path::PathBuf;

use readinfo as readinfo_lib;

pub fn run(reads: PathBuf, variants: PathBuf, output: PathBuf) -> readinfo_lib::Result<()> {
    readinfo_lib::readinfo(reads, variants, output).map(|_| ())
}
