#!/usr/bin/env nextflow

params.assembly = "hmai"
params.readbam = "$baseDir/data/*.bam"
params.readf = " "
params.readr = " "
params.outfasta = "genome.reorinted.fasta"
params.outdir = 'results'
params.mode = 'heterozygous'
params.threads = '40'
params.linreage = 'insecta'
params.busco = 'false'
params.polish = 'true'

bam_ch = Channel.fromPath(params.readbam)

process check_bam {
  container = 'mgibio/samtools:1.9'
  cpus = 1

  """
   stat ${params.readbam}
   samtools flagstat ${params.readbam}
   exit 0;
  """
}

process check_R1 {
  cpus = 1
 
  """
   stat ${params.readf}
   grep -q "^>" -m1 ${params.readf} && echo true || echo false
   exit 0;
  """
}

process check_R2 {
  cpus = 1

  """
   stat ${params.readr}
   grep -q "^>" -m1 ${params.readf} && echo true || echo false
   exit 0;
  """
}


process HiFiAdapterFilt {
  container = 'dmolik/pbadapterfilt'
  cpus = params.threads

  input:
    file bam from bam_ch
  output:
    file '*.fasta' into filt_fasta_ch
    stdout pbadapterfilt_output
  """
    pbadapterfilt.sh ${bam} -t ${task.cpus}
    echo "finished adapter filtering"
    exit 0;
  """
}

process HiFiASM {
  container = 'dmolik/hifiasm'
  cpus = params.threads

  input:
    file fasta from filt_fasta_ch
  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch
    file 'R1.fastq.gz' optional true into fastqR1_ch
    file 'R2.fastq.gz' optional true into fastqR2_ch

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
    yak count -b37 -t${task.cpus} -o pat.yak <(zcat ${params.readf}) <(zcat ${params.readf})
    yak count -b37 -t${task.cpus} -o mat.yak <(zcat ${params.readr}) <(zcat ${params.readr})
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec 1 pat.yak -2 mat.yak ${fasta}
    ln -s ${params.readf} R1.fastq.gz
    ln -s ${params.readr} R2.fastq.gz
    echo "finished alignment"
    exit 0;
  """
  else
    error "Invalid alignment mode: ${params.mode}"
}

process gzip_fastqR1_ch {
  cpus 1

  input:
    file R1 from fastqR1_ch
  output:
    file 'R1.fastq' into fastqR1_d_ch

  """
    gzip -d -c $R1 >'R1.fastq'
  """
}

process gzip_fastqR2_ch {
  cpus 1

  input:
    file R2 from fastqR2_ch
  output:
    file 'R2.fastq' into fastqR2_d_ch

  """
    gzip -d -c $R2 > 'R2.fastq'
  """
}

process gfa2fasta {
  publishDir params.outdir, mode: 'rellink'
  container = 'pvstodghill/any2fasta'
  cpus 1

  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.fasta' into gfa2fasta_fasta_res_ch
    file '*.bp.p_ctg.gfa.fasta' optional true into fasta_unoriented_ch
    file '*.bp.p_ctg.gfa.fasta' optional true into fasta_genome_ch
    file '*.p_ctg.gfa.fasta' optional true into fasta_busco_ch
    file '*hap[12].p_ctg.gfa.fasta' optional true
  """
    any2fasta ${gfa} > ${gfa}.fasta
    echo "finished gfa to fasta conversion"
    exit 0;
  """
}

process busco_gfa {
  publishDir params.outdir, mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  input:
    file fasta from fasta_busco_ch.flatten()
  output:
    file '*'

  when:
  params.busco == 'true'

  script:

  if( params.linreage == 'auto-lineage')
  """
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok')
  """
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk')
  """
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else
  """
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
}

process ragtag_dot_py {
  container = 'dmolik/ragtag'
  cpus = params.threads

  when:
  params.polish == 'true' 

  input:
    file fasta from fasta_unoriented_ch
    file fasta_ec from fasta_ec_ch
  output:
    file "${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta" into ragtag_fasta_res_ch
    file "${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta" into ragtag_fasta_genome_ch
    file "${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta" into fasta_fai_genome_ch
    file "${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta" into fasta_sshquis_genome_ch
  """
    ragtag.py patch --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_ec_patch ${fasta} ${fasta_ec}
    echo "finished patching"
    exit 0;
  """
}

process faidx {
  container = 'mgibio/samtools:1.9'
  cpus 1

  when:
  params.polish == 'true'

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
    file 'hicstuff_out/fragments_list.txt'
    file 'hicstuff_out/plots/frags_hist.pdf'
  """
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${params.readf} ${params.readr}
    echo "finished fragment calculations"
    exit 0;
  """
}

process hicstuff_polish {
  publishDir params.outdir, mode: 'rellink'
  container = 'koszullab/hicstuff'
  cpus = params.threads

  when:
  params.polish == 'true'

  input:
    file genome from ragtag_fasta_genome_ch
  output:
    file 'hicstuff_out/abs_fragments_contacts_weighted.bg2' into abs_ch
    file 'hicstuff_out/polish_fragments_list.txt'
    file 'hicstuff_out/info_contigs.txt' into contigs_ch
    file 'hicstuff_out/plots/polish_frags_hist.pdf'
  """
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${params.readf} ${params.readr}
    mv hicstuff_out/fragments_list.txt hicstuff_out/polish_fragments_list.txt
    mv hicstuff_out/plots/frags_hist.pdf hicstuff_out/plots/polish_frags_hist.pdf
    echo "finished fragment calculations"
    exit 0;
  """
}

process Shhquis_dot_jl {
  publishDir params.outdir, mode: 'rellink'
  container = 'dmolik/shhquis'
  cpus 1

  when:
  params.polish == 'true'

  input:
    file abs from abs_ch
    file contig from contigs_ch
    file genome from fasta_sshquis_genome_ch
    file fai from fai_ch
  output:
    file "${params.outfasta}" into shhquis_fasta_res_ch
    file "${params.outfasta}" into shhquis_genome_ch
    file "${params.outfasta}"

  """
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage "average"
    echo "finished reorientation"
    exit 0;
  """
}

process busco_fasta {
  publishDir params.outdir, mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  when:
  params.polish == 'true' && params.busco == 'true'

  input:
    file fasta from shhquis_genome_ch
  output:
    file '*'

  script:

  if( params.linreage == 'auto-lineage')
  """
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok')
  """
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk')
  """
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else
  """
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
}

process jellyfish {
  container = 'dmolik/jellyfish'
  cpus = params.threads

  input:
    file fastqr from fastqR1_d_ch
    file fastqf from fastqR2_d_ch
  output:
    file '*.histo' into jellyfish_histo_ch
    file 'version.txt' into jellyfish_ver_ch

  """
    jellyfish count -C -m 21 -s 1000000000 -t ${task.cpus} *.fastq -o reads.jf
    jellyfish histo -t ${task.cpus} reads.jf > ${params.assembly}.histo
    jellyfish cite > version.txt
  """
}

process genomescope2 {
  publishDir params.outdir, mode: 'rellink'
  container = 'dmolik/genomescope2'
  cpus = params.threads

  input:
    file histo from jellyfish_histo_ch
  output:
    file "${params.assembly}/*"
    file 'version.txt' into genomescope_ver_ch

  """
    genomescope.R -i ${histo} -o ${params.assembly} -k 21
    genomescope.R --version > version.txt
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

  when:
  params.polish == 'true'

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

  when:
  params.polish == 'true'

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
    echo "HiFiASM Version:"
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
    echo "any2fasta Version:"
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
    echo "Ragtag Version:"
    ragtag.py --version
    exit 0;
  """
}

process samtools_Version {
  echo "Samtools Version"
  container = 'mgibio/samtools:1.9'
  cpus 1

  output:
    stdout samtools_version

  """
    echo "Samtools Version:"
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
    echo "HiCStuff Version:"
    hicstuff --version
    exit 0;
  """
}

process bbtools_Version {
  container = 'bryce911/bbtools'
  cpus 1

  output:
    echo "BBTools Version:"
    stdout bbtools_version

  """
    stats.sh --version 2>&1
    exit 0;
  """
}

process jellyfish_Version {
  cpus 1

  input:
    file version from jellyfish_ver_ch

  output:
    stdout jellyfish_version

  """
    echo "Jellyfish Version:"
    cat $version
  """
}

process genomescope_Version {
  cpus 1

  input:
    file version from genomescope_ver_ch

  output:
    stdout genomescope_version

  """
    echo "GenomeScope Version:"
    cat $version
  """
}


process BUSCO_Version {
  cpus 1

  output:
  stdout busco_version 

  when:
  params.busco == 'true'

  """
   echo "BUSCO  - - - - - - - v5.2.2_cv1"
  """   

}

process shhquis_Version {
  cpus 1

  output:
  stdout shhquis_version

  when:
  params.polish == 'true'

  """
   echo "Shhquis.jl - - - - - 0.1.0"
  """
}

process HiFiAdapterFilt_Version {
  cpus 1

  output:
  stdout pbadapterfilt_version  

  """
   echo "HiFiAdapterFilt  - - v1.0.0"
  """
}

pbadapterfilt_output
   .collectFile(name:'filtering_information.txt', newLine: true, storeDir:"${params.outdir}/software_versions")

hifiasm_version
   .collectFile(name:'hifiasm_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

any2fasta_version
   .collectFile(name:'any2fasta_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

ragtag_version
   .collectFile(name:'ragtag_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

samtools_version
   .collectFile(name:'samtools_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

hicstuff_version
   .collectFile(name:'hicstuff_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

bbtools_version
   .collectFile(name:'bbtools_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

jellyfish_version
   .collectFile(name:'jellyfish_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

genomescope_version
   .collectFile(name:'genomescope_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

busco_version
   .collectFile(name:'busco_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

shhquis_version
   .collectFile(name:'shhquis_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

pbadapterfilt_version
   .collectFile(name:'pbadapterfilt_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }
