use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use noodles::{bam, core, sam, vcf};

pub type Error = Box<dyn std::error::Error + Send + Sync + 'static>;
pub type Result<T> = std::result::Result<T, Error>;
pub type OffsetHistogram = BTreeMap<usize, u64>;

pub fn readinfo(
    reads: impl AsRef<Path>,
    variants: impl AsRef<Path>,
    json_out: impl AsRef<Path>,
) -> Result<OffsetHistogram> {
    let histogram = read_start_offsets(reads, variants)?;
    write_histogram(json_out, &histogram)?;
    Ok(histogram)
}

pub fn read_start_offsets(
    reads: impl AsRef<Path>,
    variants: impl AsRef<Path>,
) -> Result<OffsetHistogram> {
    let mut variants_reader = vcf::reader::Builder::default().build_from_path(variants)?;
    let variants_header = variants_reader.read_header()?;

    let mut bam_reader = bam::indexed_reader::Builder::default().build_from_path(reads)?;
    let bam_header = bam_reader.read_header()?;

    let mut histogram = OffsetHistogram::new();

    for result in variants_reader.records(&variants_header) {
        let variant = result?;
        update_histogram_for_variant(&variant, &mut bam_reader, &bam_header, &mut histogram)?;
    }

    Ok(histogram)
}

fn update_histogram_for_variant(
    variant: &vcf::Record,
    bam_reader: &mut bam::IndexedReader<noodles::bgzf::Reader<File>>,
    bam_header: &sam::Header,
    histogram: &mut OffsetHistogram,
) -> Result<()> {
    let variant_position = usize::from(variant.position());
    let region = variant_region(variant)?;
    let query = bam_reader.query(bam_header, &region)?;

    for result in query {
        let record = result?;
        let Some(alignment_start) = record.alignment_start() else {
            continue;
        };

        if let Some(offset) = variant_position.checked_sub(usize::from(alignment_start)) {
            *histogram.entry(offset).or_default() += 1;
        }
    }

    Ok(())
}

fn variant_region(variant: &vcf::Record) -> Result<core::Region> {
    let position = core::Position::try_from(usize::from(variant.position()))?;
    Ok(core::Region::new(
        variant.chromosome().to_string(),
        position..=position,
    ))
}

fn write_histogram(json_out: impl AsRef<Path>, histogram: &OffsetHistogram) -> Result<()> {
    let file = File::create(json_out)?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer(&mut writer, histogram)?;
    writer.write_all(b"\n")?;
    Ok(())
}
