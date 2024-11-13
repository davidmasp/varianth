

use std::num::{NonZeroU8, NonZeroUsize};
use bstr::BString;
use bstr::BStr;
use noodles::core::Region;
use noodles::core::Position as NoodlesPosition;
use noodles::vcf::variant::record::info::field::key::END_POSITION;


trait CheckedSub: Sized {
    fn checked_sub(self, other: usize) -> Option<Self>;
}

// here we re implementing a checked sub that returns None, either if the
// result is 0 or if the result is negative
// this is needed because on the default implementation NonZeroUsize this
// is not included

impl CheckedSub for NonZeroUsize {
    fn checked_sub(self, other: usize) -> Option<Self> {
        // LEARN: get returns a u64 (primitive type)
        // this is for the negative overflow
        let val = self.get().checked_sub(other)?;
        // this is for the 0
        match val {
            0 => None,
            _ => NonZeroUsize::new(val)
        }
    }
}


// this aims at representing a position in a genome
// (should allow multi-nucleotide positions)
// (should have into-methods for noodles::core::Position)
// (should have into-methods for rust-bio)

// usize is either u32 or u64 depending on the machine,
// because u32 max is 4_294_967_295u32
// it is "possible" to have a bigger contig
// so we use u64

// for now owning, to some extend this is the same as &[u8] but 
// I don't need to keep the data around

/// `Contig` represents a contig in a genome.
///
/// # Examples
///
/// ```
/// use std::num::NonZeroUsize;
/// use bstr::BString;
/// use varianth_core::position::Contig;
///
/// let length = NonZeroUsize::new(248956422).unwrap();
/// let contig = Contig::new("chr1", length);
///
/// assert_eq!(contig.name, BString::from("chr1"));
/// assert_eq!(contig.length, length);
/// ```
#[derive(Debug)]
pub struct Contig {
    pub name: BString,
    pub length: Option<NonZeroUsize>,
}

impl Contig {
    pub fn new(name: &str, length: Option<NonZeroUsize>) -> Self {
        Self {
            name: BString::from(name),
            length,
        }
    }
    pub fn new_from_noodles(name: &bstr::BStr, length: Option<NonZeroUsize>) -> Self {
        Self {
            name: BString::from(name),
            length,
        }
    }
}

// 1-based
// width needs to be >0, width 1 -> single nucleotide

/// `Position` represents a 1-based position in a genome.
///
/// # Examples
///
/// ```
/// use std::num::{NonZeroUsize, NonZeroU8};
/// use bstr::BString;
/// use varianth_core::position::{Contig, Position};
///
/// let length = NonZeroUsize::new(248956422).unwrap();
/// let contig = Contig::new("chr1", length);
/// let position = NonZeroUsize::new(2489564).unwrap();
/// let width = NonZeroU8::new(10).unwrap();
///
/// let pos = Position::new(contig, position, width);
///
/// assert_eq!(pos.contig.name, BString::from("chr1"));
/// assert_eq!(pos.position, position);
/// assert_eq!(pos.width, width);
/// ```
#[derive(Debug)]
pub struct Position {
    pub contig: Contig,
    pub position: NonZeroUsize,
    pub width: NonZeroU8 // is this gonna be enough?
}

impl Position {
    pub fn new(contig: Contig, position: NonZeroUsize, width: NonZeroU8) -> Self {
        Self {
            contig,
            position,
            width
        }
    }
}

impl From<&Position> for Region {
    fn from(position: &Position) -> Self {
        let start_position = position.position.get();
        let start: NoodlesPosition = NoodlesPosition::try_from(start_position).unwrap();
        let end: NoodlesPosition = NoodlesPosition::try_from(start_position).unwrap();
        let contig_name = position.contig.name.to_string();
        let region: Region = Region::new(contig_name, start..=end);
        region
    }
}

pub struct CenteredPosition {
    pub contig: Contig,
    pub position: NonZeroUsize,
    pub k: usize
}

impl CenteredPosition {
    pub fn new(contig: Contig, position: NonZeroUsize, k: usize) -> Self {
        Self {
            contig,
            position,
            k
        }
    }
    pub fn get_region(&self) -> Option<Region> {
        let start_val = self.position.checked_sub(self.k)?;
        let start = NoodlesPosition::new(start_val.get())?;
        let end_val = self.position.checked_add(self.k)?;
        let end = NoodlesPosition::new(end_val.get())?;
        Some(Region::new(self.contig.name.clone(), start..=end))
    }
}

// this should be able to represent SNVs and indels
// for SNVs the width should be 1, 
// for insertions the widht should be equal to the (delta sequence)+1
// for deleteions the widht should be equal to the (delta sequence)+1
// to classify as insertion or deletion we check the ref-alt length
// SNV ref A alt T
// insertion ref A alt AT (this inserts a T after the A at position)
// deletion ref AT alt A (this deletes the T after the A at position)
// see 5.2 of vcf spec 
/// here https://samtools.github.io/hts-specs/VCFv4.2.pdf
pub struct VariantPosition {
    pub position: Position,
    pub ref_allele: BString,
    pub alt_allele: BString
}

impl VariantPosition {
    pub fn new(position: Position, ref_allele: &str, alt_allele: &str) -> Self {
        Self {
            position,
            ref_allele: BString::from(ref_allele),
            alt_allele: BString::from(alt_allele)
        }
    }
    pub fn new_from_mut(position: Position, mutstring: &BStr) -> Self {
        // convert to string
        // 62 should be the ascii code for >
        let mut parts = mutstring.splitn(2, |num| *num == 62);
        let ref_allele = parts.next().unwrap_or_else(|| panic!("Invalid mutation string"));
        let alt_allele = parts.next().unwrap_or_else(|| panic!("Invalid mutation string"));
        Self {
            position,
            ref_allele: BString::from(ref_allele),
            alt_allele: BString::from(alt_allele)
        }
    }
}
