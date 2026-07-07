# varianth - collection for VARIANT Helpers

This is a growing, actively developed and opinionated collection of
small rust tools and libraries aiming at *helping* in the analysis of genomic
variant sites.

Pronounced like "*tenth*" (with final θ).

⚠️ **This software is in experimental development and distributed as-is**.

This project relies heavily on [noodles](https://github.com/zaeleus/noodles) for genomic file format handling.

* [kmer count](#kcount---k-mer-counting)
* [ms](#ms---add-mutation-subtype-annotations)
* [vep2table](#vep2table---expand-vep-csq-annotations-to-table)
* [gene2protein](#gene2protein---genome-to-protein-nonsynonymous-mutation-expansion)


## Overview

This workspace contains multiple crates providing utilities for variant analysis:

- **varianth-cli**: Main command-line interface

Distinct commands:

- **ms**: Adds mutation subtype annotation to vcf (only for SNVs)
- **vep2table**: Utility to transform the CSQ annotation from VEP into a table
- **kmercounts**: counts kmers in selected genomic region.
- **gene2protein**: From an annotation file, creates all possible DNA mutations and its associated protein equivalents.

Other legacy components:

- **context**: Mutation subtype annotation functionality
- **hvariant**: Legacy BAM/VCF analysis tools (deprecated, being migrated)
- **mpileup-rs**: Rust implementation of mpileup functionality (in development)

---

## Installation

**Install directly from GitHub:**

```bash
cargo install --git https://github.com/davidmasp/varianth --bin varianth
```

**Or build from source:**

```bash
git clone https://github.com/davidmasp/varianth
cd varianth/varianth-cli
cargo install -path .
```

---

## Usage (via `varianth` CLI)

These commands are fully functional and available through the `varianth` CLI:

### `ms` - Add Mutation Subtype Annotations

Annotates VCF files with mutation subtype information (e.g., trinucleotide context) by extracting k-mer context around each variant from a reference genome.

**Usage:**

```bash
varianth ms \
  --fasta genome.fa \
  --variants input.vcf.gz \
  --output output.vcf.gz \
  --kval 1 \
  --feature MS \
  --featuredescription "Mutation Subtype"
```

**Parameters:**

- `-g, --fasta`: Reference genome FASTA file (must be indexed)
- `-i, --variants`: Input VCF file (must be indexed)
- `-o, --output`: Output VCF file
- `-k, --kval`: Number of adjacent bases (1 = trinucleotide context)
- `-f, --feature`: INFO field name (default: "MS")
- `-F, --featuredescription`: INFO field description (default: "Mutation Subtype")

### `kcount` - K-mer Counting

Fast k-mer counting from FASTA files with optional region-based filtering.

**Usage:**

```bash
varianth kcount \
  --size 7 \
  genome.fa \
  --output counts.json \
  --verbose
```

**Parameters:**

- `-K, --size`: K-mer size
- `-S, --table-size`: Hash table size (optional, for optimization)
- `-r, --regions`: Region string for filtering (e.g., "chr1:1000-2000")
- `-R, --regions-file`: File containing regions (🚨 not implemented yet)
- `-o, --output`: Output JSON file
- `-v, --verbose`: Enable verbose output

**Benchmarking:**

```bash
# Example comparison with jellyfish
hyperfine -m 5 --parameter-scan KMER 5 8 --warmup 2 \
  -n "jellyfish" "jellyfish count -m {KMER} -s 100M -t 1 genome.fa" \
  -n "varianth" "varianth kcount -K {KMER} genome.fa -o test.json"
```

Also see a particular usage example in the [docs/examples/cpgs](docs/counting_cpgs_in_genome.md).

### `vep2table` - Expand VEP CSQ Annotations to Table

Converts VCF records annotated with VEP (`INFO/CSQ`) into a pipe-delimited flat table, creating one output row per CSQ entry.

**Usage:**

```bash
varianth vep2table \
  --input input.vep.vcf.gz \
  --output output.tsv
```

**Parameters:**

- `-i, --input`: Input VCF/VCF.GZ file containing `CSQ` annotations
- `-o, --output`: Output table file

**Output:**

- Header starts with `chrom|pos|ref|alt|` followed by the VEP CSQ `Format:` fields from the VCF header
- One row is emitted for each CSQ annotation entry in a variant record

**Current assumptions/limitations:**

- Expects `INFO/CSQ` to be present and declared with a `Format:` section in the VCF header
- Expects exactly one ALT allele per record

### `gene2protein` - Genome-to-Protein Nonsynonymous Mutation Expansion

Generates all possible nonsynonymous single-nucleotide substitutions for coding sequences by combining genome FASTA, proteome FASTA, and CDS records from GFF.

**Usage:**

```bash
varianth gene2protein \
  --gff-path MANE.GRCh38.v1.4.ensembl_genomic.gff.gz \
  --genome-fasta-path genome.fa \
  --proteome-fasta-path MANE.GRCh38.v1.4.ensembl_protein.faa \
  --output-prefix tables/all_mutations
```

**Parameters:**

- `--gff-path`: Input GFF with CDS/protein mappings (default: `MANE.GRCh38.v1.4.ensembl_genomic.gff.gz`)
- `--genome-fasta-path`: Indexed genome FASTA path (default: `genome.fa`)
- `--proteome-fasta-path`: Indexed proteome FASTA path (default: `MANE.GRCh38.v1.4.ensembl_protein.faa`)
- `--debug-flag`: Optional limit on number of proteins processed (useful for debugging)
- `--output-prefix`: Prefix for outputs (default: `tables/all_mutations`)

**Output:**

- `<output-prefix>.tsv`: Nonsynonymous mutation table (columns: chromosome, genomic position, ref nt, alt nt, protein ID, ref aa, protein position, alt aa, codon position)
- `<output-prefix>.json`: Run metrics (`total_proteins`, `successful_count`, `failed_count`, elapsed time, failed IDs/errors)

**Notes:**

- Output directory in `--output-prefix` must already exist
- Reverse-strand CDS entries are handled with reverse complement logic before mutation expansion

Also see a particular usage example for the [MANE](https://ncbiinsights.ncbi.nlm.nih.gov/2024/10/28/mane-v1-4-mane-select-non-coding-genes/) dataset in the [docs/examples/MANE](docs/gene2protein_run_in_mane.md).

---

## 🚧 In Development

These tools are under active development and may not be fully functional or integrated into the main CLI:

### `mpileup-rs` - Rust Mpileup Implementation

A Rust reimplementation of samtools mpileup for generating pileup format from BAM files.

**Status:** Core functionality implemented but standalone binary only. Performance is currently ~24% slower than samtools (1.12s vs 0.90s on test data).

**Current capabilities:**

- Basic pileup generation with reference base
- Quality filtering (mapping quality, base quality)
- Flag-based read filtering
- Compatible output format with samtools mpileup

**Usage (standalone binary):**

## Development Roadmap

**High Priority:**

- [ ] Integrate `mpileup-rs` into main CLI
- [ ] Migration of [breadth](https://github.com/davidmasp/breadth)
- [ ] Migration of [tabix unique](https://github.com/davidmasp/tabixunique)
- [ ] Migration of [matchseq](https://github.com/davidmasp/matchseq)

---

## License

See LICENSE file for details.
