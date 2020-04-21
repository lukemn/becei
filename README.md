# Genetic map validation for the *C. becei* reference genome

This repository contains a collection of scripts for assessing concordance between a genome assembly and genetic data from a sequenced recombinant inbred line panel, run by [snakemake](https://snakemake.readthedocs.io/en/stable/).

It infers RIL ancestries by [HMM](http://genomics.princeton.edu/AndolfattoLab/MSG.html) from parental SNPs called for each RIL, and genetic mapping functions from [R/qtl](https://rqtl.org/) to test the dominant 6 linkage groups for purity and consistency.

