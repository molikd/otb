#!/usr/bin/env nextflow

params.assembly = "an_assembly"
params.readbam = "$baseDir/data/*.bam"
params.readf = "$baseDir/data/*.R1.fastq.gz"
params.readr = "$baseDir/data/*.R2.fastq.gz"
params.outfasta = "genome.out.fasta"
params.outdir = 'results'
params.mode = 'heterozygous'
params.threads = '20'
params.linreage = 'insecta'
params.busco = false
params.polish = false

bam_ch = Channel.fromPath(params.readbam)
right_fastq_check = Channel.fromPath(params.readr)
left_fastq_check = Channel.fromPath(params.readf)

bam_ch.into{ bam_check_ch; bam_Hifi_ch}

process check_bam {
  container = 'mgibio/samtools:1.9'
  cpus = 1

  input:
    file bam from bam_check_ch.flatten()
  output:
    stdout check_bam_output
  """
   touch check_bam.flag.txt
   stat ${bam}
   samtools flagstat ${bam}
   exit 0;
  """
}

process check_fastq {
  cpus = 1

  input:
    file right_fastq from right_fastq_check
    file left_fastq from left_fastq_check 
  output:
    file 'right.fastq.gz' into right_fastq_HiFiASM, right_fastq_hicstuff, right_fastq_hicstuff_polish, right_fastq_jellyfish
    file 'left.fastq.gz' into left_fastq_HiFiASM, left_fastq_hicstuff, left_fastq_hicstuff_polish, left_fastq_jellyfish
    stdout check_fastq_output
  shell:
  '''
   touch check_fastq.flag.txt
   stat !{right_fastq}
   stat !{left_fastq}

   [[ !{right_fastq}  =~ ".gz" ]] && first=$(zcat !{right_fastq} | awk '{ print $1; exit }') || first=$( cat !{right_fastq} | awk '{ print $1; exit }')
   [[ !{left_fastq} =~ ".gz" ]] && second=$(zcat !{left_fastq} | awk '{ print $1; exit }') || second=$( cat !{left_fastq} | awk '{ print $1; exit }')

   [[ $first =~ '@' ]] || exit 1;
   [[ $second =~ '@' ]] || exit 1;

   [[ !{right_fastq}  =~ ".gz" ]] && ln -s !{right_fastq} right.fastq.gz || (zcat !{right_fastq} > right.fastq.gz)
   [[ !{left_fastq}  =~ ".gz" ]] && ln -s !{left_fastq} left.fastq.gz || (zcat !{left_fastq} > left.fastq.gz)

   exit 0;
  '''
}

process HiFiAdapterFilt {
  container = 'dmolik/pbadapterfilt'
  cpus = params.threads

  input:
    file bam from bam_Hifi_ch
  output:
    file '*.fasta' into filt_fasta_ch
    stdout pbadapterfilt_output
  """
    touch pbadapterfilt.flag.txt
    pbadapterfilt.sh ${bam} -t ${task.cpus}
    echo "finished adapter filtering"
    exit 0;
  """
}

process HiFiASM {
  container = 'dmolik/hifiasm'
  cpus = params.threads

  input:
    file fasta from filt_fasta_ch.toList()
    file left from left_fastq_HiFiASM
    file right from right_fastq_HiFiASM

  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch
    stdout HiFiASM_output
  script:

  if( params.mode == 'phasing' )
  """
    touch hifiasm.flag.txt
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec --h1 ${left} --h2 ${right} ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if( params.mode == 'homozygous' )
  """
    touch hifiasm.flag.txt
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec -l0 ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if( params.mode == 'heterozygous')
  """
    touch hifiasm.flag.txt
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else if ( params.mode == 'trio')
  """
    touch hifiasm.flag.txt
    yak count -b37 -t${task.cpus} -o pat.yak <(zcat ${left}) <(zcat ${left})
    yak count -b37 -t${task.cpus} -o mat.yak <(zcat ${right}) <(zcat ${right})
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec 1 pat.yak -2 mat.yak ${fasta}
    echo "finished alignment"
    exit 0;
  """
  else
    error "Invalid alignment mode: ${params.mode}"
}

process gfa2fasta {
  publishDir "${params.outdir}/genomes", mode: 'rellink'
  container = 'pvstodghill/any2fasta'
  cpus 1

  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.fasta' into gfa2fasta_fasta_res_ch
    file '*.bp.p_ctg.gfa.fasta' optional true into fasta_unoriented_ch, fasta_genome_ch
    file '*.p_ctg.gfa.fasta' optional true into fasta_busco_ch
    file '*hap1.p_ctg.gfa.fasta' optional true into fasta_hap1_ch
    file '*hap1.p_ctg.gfa.fasta' optional true into fasta_hap2_ch
    file '*hap[12].p_ctg.gfa.fasta' optional true
    stdout gfa2fasta_output
  """
    touch any2fasta.flag.txt
    any2fasta ${gfa} > ${gfa}.fasta
    echo "finished gfa to fasta conversion"
    exit 0;
  """
}

process busco_gfa {
  publishDir "${params.outdir}/busco_no_polish", mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  input:
    file fasta from fasta_busco_ch.flatten()
  output:
    file '*'
    stdout busco_gfa_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage')
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok')
  """
    touch busco.flag.txt 
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk')
  """
    touch busco.flag.txt 
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else
  """
    touch busco.flag.txt 
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
}

process ragtag_dot_py {
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta from fasta_unoriented_ch
    file fasta_ec from fasta_ec_ch
  output:
    file "${params.assembly}_ragtag_ec_patch/ragtag.patch.fasta" into ragtag_fasta_res_ch, ragtag_fasta_genome_ch, fasta_fai_genome_ch, fasta_sshquis_genome_ch
    stdout ragtag_dot_py_output
   when:
    params.polish
  """
    touch ragtag.flag.txt 
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
   stdout faidx_output
  """
    touch faidx.flag.txt
    samtools faidx -o ${genome}.fai ${genome}
    echo "finished indexing"
    exit 0;
  """
}

process hicstuff {
  publishDir "${params.outdir}/hicstuff", mode: 'rellink'
  container = 'koszullab/hicstuff'
  cpus = params.threads

  input:
    file genome from fasta_genome_ch
    file left from left_fastq_hicstuff
    file right from right_fastq_hicstuff
  output:
    file 'hicstuff_out/fragments_list.txt'
    file 'hicstuff_out/plots/frags_hist.pdf'
    stdout hicstuff_output
  """
    touch hicstuff.flag.txt
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${left} ${right}
    echo "finished fragment calculations"
    exit 0;
  """
}

process hicstuff_polish {
  publishDir "${params.outdir}/hicstuff", mode: 'rellink'
  container = 'koszullab/hicstuff'
  cpus = params.threads

  when:
  params.polish == 'true'

  input:
    file genome from ragtag_fasta_genome_ch
    file left from left_fastq_hicstuff_polish
    file right from right_fastq_hicstuff_polish
  output:
    file 'hicstuff_out/abs_fragments_contacts_weighted.bg2' into abs_ch
    file 'hicstuff_out/polish_fragments_list.txt'
    file 'hicstuff_out/info_contigs.txt' into contigs_ch
    file 'hicstuff_out/plots/polish_frags_hist.pdf'
    stdout hicstuff_polish_output
  """
    touch hicstuff_for_polished.flag.txt 
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${left} ${right}
    mv hicstuff_out/fragments_list.txt hicstuff_out/polish_fragments_list.txt
    mv hicstuff_out/plots/frags_hist.pdf hicstuff_out/plots/polish_frags_hist.pdf
    echo "finished fragment calculations"
    exit 0;
  """
}

process Shhquis_dot_jl {
  publishDir "${params.outdir}/genomes", mode: 'rellink'
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
    file "${params.outfasta}" into shhquis_fasta_res_ch, shhquis_genome_ch, shhquis_genome_hap1_ch, shhquis_genome_hap2_ch
    file "${params.outfasta}"
    stdout Shhquis_dot_jl_output
  """
    touch shhquis.flag.txt
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage "average"
    echo "finished reorientation"
    exit 0;
  """
}

process ragtag_dot_py_hap1 {
  publishDir "${params.outdir}/genomes", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  when:
  params.polish == 'true' 

  input:
    file fasta_hap1 from fasta_hap1_ch
    file fasta_genome from shhquis_genome_hap1_ch
  output:
    file "${params.assembly}_ragtag_scaffold/hap1.ragtag.scaffold.fasta"
    file "${params.assembly}_ragtag_scaffold/hap1.ragtag.scaffold.fasta" into hap1_res_ch
    stdout ragtag_dot_py_hap1_output
  """
    touch ragtag.hap1.flag.txt
    ragtag.py scaffold --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap1}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/hap1.ragtag.scaffold.fasta
    echo "finished patching"
    exit 0;
  """
}

process ragtag_dot_py_hap2 {
  publishDir "${params.outdir}/genomes", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  when:
  params.polish == 'true' 

  input:
    file fasta_hap2 from fasta_hap2_ch
    file fasta_genome from shhquis_genome_hap2_ch
  output:
    file "${params.assembly}_ragtag_scaffold/hap2.ragtag.scaffold.fasta"
    file "${params.assembly}_ragtag_scaffold/hap2.ragtag.scaffold.fasta" into hap2_res_ch
    stdout ragtag_dot_py_hap2_output
  """
    touch ragtag.hap2.flag.txt
    ragtag.py scaffold --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap2}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/hap2.ragtag.scaffold.fasta
    echo "finished patching"
    exit 0;
  """
}

process busco_fasta {
  publishDir "${params.outdir}/busco_polish", mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  when:
  params.polish == 'true' && params.busco == 'true'

  input:
    file fasta from shhquis_genome_ch
  output:
    file '*'
    stdout busco_fasta_output
  script:

  if( params.linreage == 'auto-lineage')
  """
    touch busco_for_polished.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok')
  """
    touch busco_for_polished.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk')
  """
    touch busco_for_polished.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else
  """
    touch busco_for_polished.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_polish_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
}

process jellyfish {
  container = 'dmolik/jellyfish'
  cpus = params.threads

  input:
    file fastqr from right_fastq_jellyfish
    file fastqf from left_fastq_jellyfish
  output:
    file '*.histo' into jellyfish_histo_ch
    file 'version.txt' into jellyfish_ver_ch
    stdout jellyfish_output
  """
    touch jellyfish.flag.txt
    jellyfish count -C -m 21 -s 1000000000 -t ${task.cpus} -o reads.jf <(zcat ${fastqr}) <(zcat ${fastqf}) 
    jellyfish histo -t ${task.cpus} reads.jf > ${params.assembly}.histo
    jellyfish cite > version.txt
  """
}

process genomescope2 {
  publishDir "${params.outdir}/genomescope", mode: 'rellink'
  container = 'dmolik/genomescope2'
  cpus = params.threads

  input:
    file histo from jellyfish_histo_ch
  output:
    file "${params.assembly}/*"
    file 'version.txt' into genomescope_ver_ch
    stdout genomescope2_output
  """
    touch genomescope.flag.txt
    xvfb-run genomescope.R -i ${histo} -o ${params.assembly} -k 21
    genomescope.R --version > version.txt
  """
}

process gfa2fasta_stats_dot_sh {
  publishDir "${params.outdir}/genomes", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from gfa2fasta_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch any2fasta_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process ragtag_stats_dot_sh {
  publishDir "${params.outdir}/genomes", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  when:
  params.polish == 'true'

  input:
    file fasta from ragtag_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch ragtag_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process sshquis_stats_do_sh {
  publishDir "${params.outdir}/genomes", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  when:
  params.polish == 'true'

  input:
    file fasta from shhquis_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch shhquis_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process ragtag_stats_dot_sh_hap1 {
  publishDir "${params.outdir}/genomes", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  when:
  params.polish == 'true'

  input:
    file fasta from hap1_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch ragtag_hap1_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process ragtag_stats_dot_sh_hap2 {
  publishDir "${params.outdir}/genomes", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  when:
  params.polish == 'true'

  input:
    file fasta from hap2_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch ragtag_hap2_stats.flag.txt
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
    touch HiFiASM_version.flag.txt
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
    touch any2fasta_version.flag.txt
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
    touch ragtag_version.flag.txt
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
    touch samtools_version.flag.txt
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
    touch hicstuff_version.flag.txt
    echo "HiCStuff Version:"
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
    touch bbtools_version.flag.txt
    echo "BBTools Version"
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
    touch jellyfish_version.flag.txt
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
    touch genomescope_version.flag.txt
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
   touch busco_version.flag.txt
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
   touch shhquis_version.flag.txt
   echo "Shhquis.jl - - - - - 0.1.0"
  """
}

process HiFiAdapterFilt_Version {
  cpus 1

  output:
  stdout pbadapterfilt_version  

  """
   touch hifiadapterfilt_version.flag.txt
   echo "HiFiAdapterFilt  - - v1.0.0"
  """
}

pbadapterfilt_output
   .collectFile(name:'filtering_information.log.txt', newLine: true, storeDir:"${params.outdir}/filtering")

check_fastq_output
   .collectFile(name:'fastq_check.log.txt', newLine: true, storeDir:"${params.outdir}/filtering")

check_bam_output
   .collectFile(name:'bam_check.log.txt', newLine: true, storeDir:"${params.outdir}/filtering")

HiFiASM_output
   .collectFile(name:'HiFiASM.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

gfa2fasta_output
   .collectFile(name:'gfa2fasta.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

hicstuff_output
   .collectFile(name:'hicstuff.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

hicstuff_polish_output
   .collectFile(name:'hicstuff_polish.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

ragtag_dot_py_output
   .collectFile(name:'ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

faidx_output
   .collectFile(name:'faidx.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

Shhquis_dot_jl_output
   .collectFile(name:'shhquis.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

ragtag_dot_py_hap1_output
   .collectFile(name:'ragtag_hap1.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

ragtag_dot_py_hap2_output
   .collectFile(name:'ragtag_hap2.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

busco_gfa_output
   .collectFile(name:'busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_no_polish")

busco_fasta_output
   .collectFile(name:'busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_polish" )

jellyfish_output
   .collectFile(name:'jellyfish.log.txt', newLine: true, storeDir:"${params.outdir}/genomescope" )

genomescope2_output
   .collectFile(name:'genomescope2.log.txt', newLine: true, storeDir:"${params.outdir}/genomescope" )

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
