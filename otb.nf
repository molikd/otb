#!/usr/bin/env nextflow

params.assembly = "hmai"
params.bam = "$baseDir/data/*.bam"
params.readf = "$baseDir/data/*_R1.fastq.gz"
params.readr = "$baseDir/data/*_R2.fastq.gz"
params.outfasta = "genome.reorinted.fasta"
params.outdir = 'results'
params.threads = 40

bam_ch = Channel.fromPath(params.bam)
reads_ch = Channel.fromFilePairs(params.reads)

process stats.sh {
  publishDir params.outdir, mode: 'copy'

  input:
    file fasta from fasta_res_ch
  output:
    file '*.stats'
  """
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
  """
}

process HiFiAdapterFilt {
  input:
    file bam from bam_ch
  output:
    file '*.fasta' into filt_fasta_ch
  """
    pbadapterfilt.sh ${bam} -t ${params.threads}
  """
}

process HiFiASM {
  input:
    file fasta from filt_fasta_ch
  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' in fasta_ec_ch
  """
    hifiasm -o ${params.assembly} -t ${params.threads} --write-paf --write-ec --h1 ${params.readf} --h2 ${params.readr} ${fasta}
  """
}

process gfa2fasta {
  input:
    file gfa from gfa_ch
  output:
    file '*.fasta' into fasta_res_ch
    file '*.hic.p_ctg.fasta' optional true into fasta_unoriented_ch
  """
    any2fasta ${gfa} > $(echo ${gfa} | sed 's/\.gfa/\.fasta/g')
  """
}

process ragtag.py {
  input:
    file fasta from fasta_unoriented_ch
    file fasta_ec from fasta_ec_ch
  output:
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_res_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_genome_ch
  """
    ragtag.py patch --aligner unimap -t ${params.threads} -o ./${params.assembly}_ragtag_ec_patch ${fasta} ${fasta_ec}
  """
}

process faidx {
  input:
   file genome from fasta_genome_ch
  ouput:
   file '*.fai' info fai_ch
  """
    samtools faidx -o ${genome}.fai ${genome}
  """
}

process hicstuff {
  publishDir params.outdir, mode: 'rellink'

  input:
    file genome from fasta_genome_ch
  output:
    file 'hicstuff_out/abs_fragments_contacts_weighted.bg2' into abs_ch
    file 'hicstuff_out/fragments_list.txt'
    file 'hicstuff_out/info_contigs.txt' into contigs_ch
    file 'hicstuff_out/plots/frags_hist.pdf'
  """
    hicstuff pipeline -t ${params.threads} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${params.readf} ${params.readr}
  """
}

process Shhquis.jl {
  publishDir params.outdir, mode: 'rellink'

  input:
    file abs from abs_ch
    file contig from contigs_ch
    file genome from fasta_genome_ch
    file fai from fai_ch

  """
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage "average"
  """
}
