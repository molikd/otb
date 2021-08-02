#!/usr/bin/env nextflow

params.outdir = "$baseDir/otb"
params.bam = "$baseDir/hic.bam"
params.threads = 8

process stats.sh {}
process HiFiAdapterFilt {}
process HiFiASM {}
process gfa2fasta {}
process ragtag.py {}
process hicstuff {}
process Shhquis.jl {}
