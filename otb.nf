#!/usr/bin/env nextflow

params.assembly = "hmai"
params.readbam = "$baseDir/data/*.bam"
params.readf = "$baseDir/data/*_R1.fastq.gz"
params.readr = "$baseDir/data/*_R2.fastq.gz"
params.outfasta = "genome.reorinted.fasta"
params.outdir = 'results'
params.mode = 'homozygous'

/*TODO gotta be a better way to handle cluster stuff
*/
params.threads = 40

bam_ch = Channel.fromPath(params.readbam)

process HiFiAdapterFilt {
  input:
    file bam from bam_ch
  output:
    file '*.fasta' into filt_fasta_ch
  """
    pbadapterfilt.sh ${bam} -t ${params.threads}
    echo "finished adapter filtering"
    exit 0;
  """
}

process HiFiASM {
  input:
    file fasta from filt_fasta_ch
  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch
  script:

  if( params.mode == 'phasing' )
  """
    hifiasm -o ${params.assembly} -t ${params.threads} --write-paf --write-ec --h1 ${params.readf} --h2 ${params.readr} ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if( params.mode == 'homozygous' )
  """
    hifiasm -o ${params.assembly} -t ${params.threads} --write-paf --write-ec -l0 ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if( params.mode == 'heterozygous')
  """
    hifiasm -o ${params.assembly} -t ${params.threads} --write-paf --write-ec ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if ( params.mode == 'trio')
  """
    yak count -b37 -t${params.threads} -o pat.yak <(cat ${params.readf}) <(cat ${params.readf})
    yak count -b37 -t${params.threads} -o mat.yak <(cat ${params.readr}) <(cat ${params.readr})
    hifiasm -o ${params.assembly} -t ${params.threads} --write-paf --write-ec 1 pat.yak -2 mat.yak ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else
    error "Invalid alignment mode: ${params.mode}"
}

process gfa2fasta {
  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.fasta' into gfa2fasta_fasta_res_ch
    file '*.hic.p_ctg.fasta' optional true into fasta_unoriented_ch
  """
    any2fasta ${gfa} > ${gfa}.fasta
    echo "finished gfa to fasta conversion"
    exit 0;
  """
}

process ragtag_dot_py {
  input:
    file fasta from fasta_unoriented_ch
    file fasta_ec from fasta_ec_ch
  output:
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into ragtag_fasta_res_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_genome_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_fai_genome_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_sshquis_genome_ch
  """
    ragtag.py patch --aligner unimap -t ${params.threads} -o ./${params.assembly}_ragtag_ec_patch ${fasta} ${fasta_ec}
    echo "finished patching"
    exit 0;
  """
}

process faidx {
  input:
   file genome from fasta_fai_genome_ch
  output:
   file '*.fai' into fai_ch
  """
    samtools faidx -o ${genome}.fai ${genome}
    echo "finished indexing"
    exit 0;
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
    echo "finished fragment calculations"
    exit 0;
  """
}

process Shhquis_dot_jl {
  publishDir params.outdir, mode: 'rellink'

  input:
    file abs from abs_ch
    file contig from contigs_ch
    file genome from fasta_sshquis_genome_ch
    file fai from fai_ch
  output:
    file '${params.outfasta}' into shhquis_fasta_res_ch

  """
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage "average"
    echo "finished reorientation"
    exit 0;
  """
}

process gfa2fasta_stats_dot_sh {
  publishDir params.outdir, mode: 'copy'

  input:
    file fasta from gfa2fasta_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process ragtag_stats_dot_sh {
  publishDir params.outdir, mode: 'copy'

  input:
    file fasta from ragtag_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process sshquis_stats_do_sh {
  publishDir params.outdir, mode: 'copy'

  input:
    file fasta from shhquis_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}
