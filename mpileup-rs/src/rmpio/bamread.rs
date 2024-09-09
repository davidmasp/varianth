use std::fs::File;

use noodles::bam::io;
use noodles::bam::io::Reader as BamReader;
use noodles::bgzf::Reader as GzReader;
use noodles::sam::Header;

pub fn reader_from_path(
    file_path: &str,
) -> Result<(BamReader<GzReader<File>>, Header), std::io::Error> {
    let mut reader: BamReader<GzReader<File>> =
        io::reader::Builder::default().build_from_path(file_path)?;
    let header = reader.read_header()?;
    Ok((reader, header))
}
