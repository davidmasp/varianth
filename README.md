# varianth - collection for VARIANT Helpers

This is a growing, actively developed and opinionated collection of
Rust tools and libraries aiming at helping the analyzing genomic variant sites.

_Pronounced like "tenth" (with final θ)_

⚠️ **This is experimental development** - APIs and functionality may change.

This project relies heavily on [noodles](https://github.com/zaeleus/noodles) for genomic file format handling.

## Overview

This workspace contains multiple crates providing utilities for variant analysis:
- **varianth-cli**: Main command-line interface (active development)
- **varianth-core**: Core data structures and utilities
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
cd varianth
cargo build --release
```

The main binary will be available at `target/release/varianth`.

---

## ✅ Ready to Use (via `varianth` CLI)

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

🚨 TO INTEGRATE INTO VARIANTH

## 📚 Legacy Tools (Being Migrated)

The following tools exist in the `hvariant` crate but are deprecated and being migrated to the main `varianth` CLI:

### `readinfo` - Variant Position in Reads Histogram

**Note:** Currently only available in legacy `hvariant` binary (not built by default).

Analyzes BAM files to generate histograms showing where variants appear within reads.

**Legacy Usage:**
```bash
hvariant readinfo \
  --reads sample.bam \
  --variants variants.vcf.gz \
  --outfile output.json
```

**Output:** JSON file with position histograms

**Performance (1Mb germline data):**
- Records: 11,281
- Time: 1,260 seconds
- Memory: < 1GB

### `readfreq` - Read Frequency at Positions

**Note:** Currently only available in legacy `hvariant` binary (not built by default).

Extracts sequences and read counts from BAM files at specified positions.

**Legacy Usage:**
```bash
hvariant readfreq \
  --reads sample.bam \
  --variants positions.bed \
  --outfile output.tsv
```

**Input:** BED file (3 columns) + indexed BAM file

**Output:** TSV with format:
```
chr     start   end     sequence        count
20      47000001        47000003        CAA     4
20      47100001        47100003        CTG     5
```

**Limitations:**
- Reads with hard-clipped or pan CIGAR operations are excluded
- Sequences extracted from reads, not reference


## Development Roadmap

**High Priority:**
- [ ] Integrate `mpileup-rs` into main CLI
- [ ] Migrate `readinfo` and `readfreq` to `varianth` CLI
- [ ] Add comprehensive tests

- [ ] Migration of [breadth](https://github.com/davidmasp/breadth)
- [ ] Migration of [tabix unique](https://github.com/davidmasp/tabixunique)
- [ ] Migration of [matchseq](https://github.com/davidmasp/matchseq)

**Future Enhancements:**
- [ ] Parallel processing support for multiple chromosomes
- [ ] Streaming VCF processing
- [ ] Additional variant annotation types
- [ ] Multi-sample support

---

## License

See LICENSE file for details.

