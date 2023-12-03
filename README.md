# Helpers for VARIANT analysis (hvariant)


This is a set of small rust tools that can be helpful in the analysis
of genomic variants, written in rust.

:warning: This is experimental development

## Tools

### Readinfo

Run with:

```bash
# bam needs to be indexed
havariant readinfo --reads path/to/file.bam --variants path/to/variants.vcf.gz
```

It returns a json file (`output.json`, see `-o`) that contains a histogram of position of the variants
in the read.

The current estimate for 1Mb (_germline_) performance is:

* Number of records: `11281`
* Time (seconds): `1260`
* Memory: `< 1G`

Note: this estimates can change dramatically with different conditions.

### Readfreq

It returns a table file (`out.tsv`, see `-o`) that contains the position in
bed format (0-based) and two extra columns that indicate the sequence and
number of reads. Sequence is extracted from the bam file, not the
reference sequence.

The input needs to be a bed (N=3) file and a bam file.

:warning: **note** that because I am not sure how to treat the
hard clipped and pan cigar operations, reads with 
any of such bases are not considered.

#### Example

```bash
bedfile="devdata/small.bed"
bamfile="devdata/wgs1kg/HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"
hvariant readfreq -r ${bedfile} -v ${bamfile}  
```

```bash
$ more out.tsv 
20      47000001        47000003        CAA     4
20      47100001        47100003        CTG     5
20      47099956        47099958        TCG     1
20      47099956        47099958        TAG     4
```

