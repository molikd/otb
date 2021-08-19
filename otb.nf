#!/usr/bin/env nextflow

params.assembly = "hmai"
params.readbam = "$baseDir/data/*.bam"
params.readf = "$baseDir/data/*_R1.fastq.gz"
params.readr = "$baseDir/data/*_R2.fastq.gz"
params.outfasta = "genome.reorinted.fasta"
params.outdir = 'results'
params.mode = 'heterozygous'
params.threads = '40'

bam_ch = Channel.fromPath(params.readbam)

process HiFiAdapterFilt {
  container = 'dmolik/pbadapterfilt'
  cpus = params.threads

  input:
    file bam from bam_ch
  output:
    file '*.fasta' into filt_fasta_ch
  """
    pbadapterfilt.sh ${bam} -t ${task.cpus}
    echo "finished adapter filtering"
    exit 0;
  """
}

process HiFiASM {
  publishDir params.outdir, mode: 'rellink'
  container = 'dmolik/hifiasm'
  cpus = params.threads

  input:
    file fasta from filt_fasta_ch
  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch
    file '*hap[12].p_ctg.gfa.fasta' optional true
  script:

  if( params.mode == 'phasing' )
  """
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec --h1 ${params.readf} --h2 ${params.readr} ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if( params.mode == 'homozygous' )
  """
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec -l0 ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if( params.mode == 'heterozygous')
  """
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if ( params.mode == 'trio')
  """
    yak count -b37 -t${task.cpus} -o pat.yak <(cat ${params.readf}) <(cat ${params.readf})
    yak count -b37 -t${task.cpus} -o mat.yak <(cat ${params.readr}) <(cat ${params.readr})
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec 1 pat.yak -2 mat.yak ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else
    error "Invalid alignment mode: ${params.mode}"
}

process gfa2fasta {
  container = 'pvstodghill/any2fasta'
  cpus 1

  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.fasta' into gfa2fasta_fasta_res_ch
    file '*.p_ctg.gfa.fasta' optional true into fasta_unoriented_ch
  """
    any2fasta ${gfa} > ${gfa}.fasta
    echo "finished gfa to fasta conversion"
    exit 0;
  """
}

process ragtag_dot_py {
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta from fasta_unoriented_ch
    file fasta_ec from fasta_ec_ch
  output:
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into ragtag_fasta_res_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_genome_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_fai_genome_ch
    file './${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta' into fasta_sshquis_genome_ch
  """
    ragtag.py patch --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_ec_patch ${fasta} ${fasta_ec}
    echo "finished patching"
    exit 0;
  """
}

process faidx {
  container = 'dmolik/samtools'
  cpus 1

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
  container = 'koszullab/hicstuff'
  cpus = params.threads

  input:
    file genome from fasta_genome_ch
  output:
    file 'hicstuff_out/abs_fragments_contacts_weighted.bg2' into abs_ch
    file 'hicstuff_out/fragments_list.txt'
    file 'hicstuff_out/info_contigs.txt' into contigs_ch
    file 'hicstuff_out/plots/frags_hist.pdf'
  """
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${params.readf} ${params.readr}
    echo "finished fragment calculations"
    exit 0;
  """
}

process Shhquis_dot_jl {
  publishDir params.outdir, mode: 'rellink'
  container = 'dmolik/shhquis'
  cpus 1

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
  container = 'bryce911/bbtools'
  cpus 1

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
  container = 'bryce911/bbtools'
  cpus 1

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
  container = 'bryce911/bbtools'
  cpus 1

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

process HiFiASM_Version {
  container = 'dmolik/hifiasm'
  cpus 1

  output:
    stdout hifiasm_version

  """
    hifiasm --version
    exit 0;
  """
}

process any2fasta_Version {
  container = 'pvstodghill/any2fasta'
  cpus 1

  output:
    stdout any2fasta_version

  """
    any2fasta -v
    exit 0;
  """
}

process ragtag_Version {
  container = 'dmolik/ragtag'
  cpus 1

  output:
    stdout ragtag_version

  """
    ragtag.py --version
    exit 0;
  """
}

process samtools_Version {
  container = 'dmolik/samtools'
  cpus 1

  output:
    stdout samtools_version

  """
    samtools --version
    exit 0;
  """
}

process hicstuff_Version {
  container = 'koszullab/hicstuff'
  cpus 1

  output:
    stdout hicstuff_version

  """
    hicstuff --version
    exit 0;
  """
}

process bbtools_Version {
  container = 'bryce911/bbtools'
  cpus 1

  output:
    stdout bbtools_version

  """
    stats.sh --version
    exit 0;
  """
}

process Other_Version {
  cpus 1

  output:
    stdout other_version

  """
    echo "Shhquis.jl - - - - - 0.1.0
    echo "HiFiAdapterFilt  - - v1.0.0
  """
}

pbadapterfilt_version.subscribe {
  println "HiFiAdapterFilt Version"
  println "$it"
}

hifiasm_version.subscribe {
  println "HiFiASM Version"
  println "$it"
}

any2fasta_version.subscribe {
  println "any2fasta Version"
  println "$it"
}

ragtag_version.subscribe {
  println "RagTag Version"
  println "$it"
}

samtools_version.subscribe {
  println "Samtools Version"
  println "$it"
}

hicstuff_version.subscribe {
  println "hicstuff Version"
  println "$it"
}

bbtools_version.subscribe {
  println "BBMap Version"
  println "$it"
}

other_version.subscribe {
  println "Other Versions"
  ptintln "$it"
}
