# CALCULATE COVERAGE OF EXOME ALIGNMENT

This is an Rscript to plot and print out the coverage metrics of an exome alignment.
The input file is a gzip compressed per-base coverage file produced by bedtools:

`coverageBed -abam exome.bam -b intervals.bed -d | gzip > sample.per.base.coverage.gz`

The output is a histogram plot showing the mean and median coverage, and percent covered at a minimum read depth of 10, 20 and 50 reads.
The sex is calculated by comparing the ratio of average coverage on Y chromosome:X chromosome. A ratio of less than 0.3 is a female.

A table of the metrics is also printed out.

## Usage:

`Rscript exome.coverage.R sample.per.base.coverage.gz`
