# README

To render this README.md file run:

    Rscript -e "rmarkdown::render('README.Rmd')"

## Feature parity and benchmarks

### Pairity test A

This test

    /Users/dmas/bin/hyperfine --warmup 2 \
        --command-name samtoolsA \
        --output ./pairitytest_basicpileupA/samtools.mp \
        --export-json pairitytest_basicpileupA/samtools.json \
        "/Users/dmas/bin/samtools mpileup -B -f ./hg38.fa --no-output-del --no-output-del --min-MQ 24 --max-depth 0 -Q 0 -x --reverse-del --no-output-ends minsample.bam"

    /Users/dmas/bin/hyperfine --warmup 2 \
        --command-name mprsA \
        --max-runs 5 \
        --output ./pairitytest_basicpileupA/mprs.mp \
        --export-json pairitytest_basicpileupA/mprs.json \
        "./target/release/mprs"

    ## Benchmark 1: samtoolsA
    ##   Time (mean ± σ):     901.2 ms ±  34.1 ms    [User: 808.1 ms, System: 68.9 ms]
    ##   Range (min … max):   872.4 ms … 972.3 ms    10 runs
    ##  
    ## Benchmark 1: mprsA
    ##   Time (mean ± σ):      1.121 s ±  0.015 s    [User: 0.880 s, System: 0.215 s]
    ##   Range (min … max):    1.103 s …  1.142 s    5 runs
    ## 

### Benchmarks results

<table>
<thead>
<tr class="header">
<th style="text-align: left;">name</th>
<th style="text-align: left;">mean</th>
<th style="text-align: left;">median</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">samtools_testA</td>
<td style="text-align: left;">0.90</td>
<td style="text-align: left;">0.89</td>
</tr>
<tr class="even">
<td style="text-align: left;">mprs_testA</td>
<td style="text-align: left;">1.12</td>
<td style="text-align: left;">1.12</td>
</tr>
</tbody>
</table>

### Parity checks

    ## Rows: 3 Columns: 3
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): test, type
    ## dbl (1): diff_error_code
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

<table>
<thead>
<tr class="header">
<th style="text-align: left;">test</th>
<th style="text-align: left;">positions covered</th>
<th style="text-align: left;">depth</th>
<th style="text-align: left;">whole file</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">A</td>
<td style="text-align: left;">✅</td>
<td style="text-align: left;">✅</td>
<td style="text-align: left;">✅</td>
</tr>
</tbody>
</table>
