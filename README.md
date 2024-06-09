# collection for VARIANT Helpers (varianth)


Set of rust tools and functions that aim at helping
the analysis of genomic variant sites.

pronounced like tenth (with final &Theta;)

:warning: This is experimental development

This crate relies massively in noodles, although the structs might
look sometimes duplicated here we aim at generating a unique set of
types that is sharable across the package while in noodles this is not
necessarly the case.

```

hyperfine -m 5 --parameter-scan KMER 5 8 --warmup 2 -n "jelly_5" -n "jelly_6" -n "jelly_7" -n "jelly_8" "jellyfish count -m {KMER} -s 100M -t 1 Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa && jellyfish dump mer_counts.jf > mer_counts_dumps.fa" --export-json jelly.json

hyperfine -m 5 --parameter-scan KMER 5 8 --warmup 2 -n "var_5" -n "var_6" -n "var_7" -n "var_8" "~/projects/varianth/target/release/varianth kcount -K {KMER} Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa -o test.json" --export-json var.json


```


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

