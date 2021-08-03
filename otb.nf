#!/usr/bin/env nextflow

params.outdir = "$baseDir/otb"
params.bam = "$baseDir/hic.bam"
params.outfasta = "genome.reorinted.fasta"
params.threads = 1

bam_ch = Channel.fromPath(params.bam)

/*
process stats.sh {

}
*/

process HiFiAdapterFilt {
  input:
    file bam from bam_ch
  output:
    file '*.fasta' into filt_fasta_ch
  """
    pbadapterfilt.sh ${bam} -t ${params.threads}
  """
}

process HiFiASM {}

process gfa2fasta {}

process ragtag.py {}

process hicstuff {}

process Shhquis.jl {
 input:
 file unorganized

  """
    shh.jl --reorient "Hemileuca_maia.reorient.fasta" --genome "genome.fasta" --fai "genome.fasta.fai" --bg2 "abs_fragments_contacts_weighted.bg2" --contig "info_contigs.txt" --hclust-linkage "average"  """
}
*/
