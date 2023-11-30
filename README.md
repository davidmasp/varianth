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

Note: this estimates can change dramatically.

