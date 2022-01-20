#!/usr/bin/env nextflow

params.assembly = "an_assembly"
params.readf = "$baseDir/data/*.R1.fastq.gz"
params.readr = "$baseDir/data/*.R2.fastq.gz"
params.outfasta = "genome.out.fasta"
params.outdir = 'results'
params.threads = '21'
params.linreage = 'insecta_odb10'
params.busco = false
params.buscooffline = false
params.buscodb = "/work/busco"

right_fastq_check = Channel.fromPath(params.readr)
left_fastq_check = Channel.fromPath(params.readf)

process check_fastq {
  cpus = 1

  input:
    file right_fastq from right_fastq_check
    file left_fastq from left_fastq_check 
  output:
    file 'out/right.fastq.gz' into right_fastq_HiFiASM, right_fastq_hicstuff_polish, right_fastq_jellyfish
    file 'out/left.fastq.gz' into left_fastq_HiFiASM, left_fastq_hicstuff_polish, left_fastq_jellyfish
    file 'out/*.fastq.gz' into fasta_in_ch
    stdout check_fastq_output
  shell:
  '''
   state () { printf "%b\n" "[$(date)]: $*" 2>&1; }
   error () { printf "%b\n" "[$(date)]: $*" >&2; exit 1; }

   state "use touch to create a flag, so that I can be found easily"
   touch check_fastq.flag.txt
   state "stat !{right_fastq}..."
   stat !{right_fastq}
   state "stat !{left_fastq}..."
   stat !{left_fastq}

   state "if !{right_fastq} ends in gz zcat, and if it does not, cat the first line, save."
   [[ !{right_fastq}  =~ ".gz" ]] && first=$(zcat !{right_fastq} | awk '{ print $1; exit }') || first=$( cat !{right_fastq} | awk '{ print $1; exit }')
   state "if !{left_fastq} end in gz zcat, and if it does not, cat the first line, save"
   [[ !{left_fastq} =~ ".gz" ]] && second=$(zcat !{left_fastq} | awk '{ print $1; exit }') || second=$( cat !{left_fastq} | awk '{ print $1; exit }')

   state "if the first line of the fastqs doesn't start with an @, error out, otherwise continue"
   [[ $first =~ '@' ]] || error "!{right_fastq} doesn't start with an @";
   [[ $second =~ '@' ]] || errror "!{left_fastq} doesn't start with an @";

   state "check to make sure that fastqs are divisable by 4"
   [[ !{right_fastq}  =~ ".gz" ]] && first=$(zcat !{right_fastq} | wc -l) || first=$( cat !{right_fastq} | wc -l)
   [[ !{left_fastq} =~ ".gz" ]] && second=$(zcat !{left_fastq} | wc -l) || second=$( cat !{left_fastq} | wc -l )

   [[ $(( $first % 4 )) -eq 0 ]] || error "number of lines in !{right_fastq} not divisable by four"
   [[ $(( $second % 4 )) -eq 0 ]] || error "number of lines in !{left_fastq} not divisable by four"

   state "make softlinks for both files"
   mkdir out

   if [[ !{right_fastq}  =~ ".gz" ]]; then
     cd out
     ln -s ../!{right_fastq} right.fastq.gz
     cd ..
   else
     gzip -c !{right_fastq} > out/right.fastq.gz
   fi

   if [[ !{left_fastq}  =~ ".gz" ]]; then
     cd out
     ln -s ../!{left_fastq} left.fastq.gz
     cd ..
   else
     gzip -c !{left_fastq} > out/left.fastq.gz
   fi

   state "successful completion"
   sleep 120;
   exit 0;
  '''
}

process HiFiASM {
  container = 'dmolik/hifiasm'
  cpus = params.threads

  input:
    file left from left_fastq_HiFiASM
    file right from right_fastq_HiFiASM

  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch
    stdout HiFiASM_output
  script:
  """
    touch hifiasm.flag.txt
    hifiasm -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec ${left} ${right} 2>&1
    echo "finished alignment"
    sleep 120;
    exit 0;
  """
}

process gfa2fasta {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'pvstodghill/any2fasta'
  cpus 1

  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.p_ctg.gfa.fasta' optional true into gfa2fasta_fasta_res_ch
    file '*.bp.p_ctg.gfa.fasta' optional true into fasta_unoriented_ch, fasta_busco_ch
    file '*hap[12].p_ctg.gfa.fasta' optional true into fasta_hap_ch, simple_fasta_hap_polish_ch
    stdout gfa2fasta_output
  """
    touch any2fasta.flag.txt
    any2fasta ${gfa} > ${gfa}.fasta
    echo "finished gfa to fasta conversion"
    sleep 120;
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

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
    exit 0;
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt 
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
    exit 0;
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt 
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
    exit 0;
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt 
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
    exit 0;
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
    exit 0;
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
    exit 0;
  """
  else
  """
    touch busco.flag.txt
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
  """
    touch ragtag.flag.txt 
    ragtag.py patch --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_ec_patch ${fasta} ${fasta_ec}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process faidx {
  container = 'mgibio/samtools:1.9'
  cpus 1

  input:
    file genome from fasta_fai_genome_ch
  output:
    file '*.fai' into fai_ch
    stdout faidx_output
  """
    touch faidx.flag.txt
    samtools faidx -o ${genome}.fai ${genome}
    echo "finished indexing"
    sleep 120;
    exit 0;
  """
}

process hicstuff_polish {
  publishDir "${params.outdir}/hicstuff", mode: 'rellink'
  container = 'koszullab/hicstuff'
  cpus = params.threads

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
    sleep 120;
    exit 0;
  """
}

process Shhquis_dot_jl {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'dmolik/shhquis'
  cpus 1

  input:
    file abs from abs_ch
    file contig from contigs_ch
    file genome from fasta_sshquis_genome_ch
    file fai from fai_ch
  output:
    file "${params.outfasta}" into shhquis_fasta_res_ch, polish_haps_genome_ch, shhquis_simple_ch
    file "${params.outfasta}"
    stdout Shhquis_dot_jl_output
  """
    touch shhquis.flag.txt
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage "average"
    echo "finished reorientation"
    sleep 120;
    exit 0;
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
    sleep 120;
    exit 0;
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
    file "kcov.txt" into kcov_ch
    file "${params.assembly}/lookup_table.txt" into lookup_table_ch
    file 'version.txt' into genomescope_ver_ch
    stdout genomescope2_output
  """
    touch genomescope.flag.txt
    xvfb-run genomescope.R -i ${histo} -o ${params.assembly} -k 21 -p ${params.ploidy} --fitted_hist 
    genomescope.R --version > version.txt
    awk '/kmercov [0-9]/ { print \$2 }' ${params.assembly}/model.txt >> kcov.txt
    echo "finished genomescope"
    sleep 120;
    exit 0;
  """
}

process simple_polish {
  cpus = 1

  input:
    file genome from shhquis_simple_ch
  output:
    file "${params.assembly}.polished.genome.fasta" into simple_polished_genome_ch, simple_polished_genome_busco_ch
  """
    touch simple_polish.flag.txt
    ln -s ${genome} ${params.assembly}.polished.genome.fasta
    echo "finished softlink"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_hap_simple_polish {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from simple_fasta_hap_polish_ch.flatten()
    file fasta_genome from simple_polished_genome_ch
  output:
    file "polished*fasta"
    file "polished*fasta" into simple_hap_patch_res_ch
    stdout simple_ragtag_dot_py_hap_output
  """
    touch ragtag.hap.flag.txt
    ragtag.py scaffold --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process simple_busco_fasta {
  publishDir "${params.outdir}/busco_polish", mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  input:
    file fasta from simple_polished_genome_busco_ch
  output:
    file '*'
    stdout simple_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process fasta_in_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from  fasta_in_ch.flatten()
  output:
    file '*.stats'
  """
    touch any2fasta_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process gfa2fasta_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from  gfa2fasta_fasta_res_ch.flatten()
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
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

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

process shhquis_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

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

process simple_hap_patch_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from simple_hap_patch_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch hap_patch_stats.flag.txt
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

process bcftools_Version {
  container = 'mgibio/bcftools:1.9'
  cpus 1

  output:
    stdout bcftools_version

  """
    touch bcftools_version.flag.txt
    echo "Bcftools Version:"
    bcftools --version
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
    params.busco

  """
    touch busco_version.flag.txt
    echo "BUSCO  - - - - - - - v5.2.2_cv1"
    exit 0;
  """   

}

process shhquis_Version {
  cpus 1

  output:
    stdout shhquis_version

  """
    touch shhquis_version.flag.txt
    echo "Shhquis.jl - - - - - 0.1.0"
    exit 0;
  """
}

pbadapterfilt_output
   .collectFile(name:'filtering_information.log.txt', newLine: true, storeDir:"${params.outdir}/filtering")

check_fastq_output
   .collectFile(name:'fastq_check.log.txt', newLine: true, storeDir:"${params.outdir}/filtering")

check_in_file_output
   .collectFile(name:'bam_check.log.txt', newLine: true, storeDir:"${params.outdir}/filtering")

HiFiASM_output
   .collectFile(name:'HiFiASM.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

gfa2fasta_output
   .collectFile(name:'gfa2fasta.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

hicstuff_polish_output
   .collectFile(name:'hicstuff_polish.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

ragtag_dot_py_output
   .collectFile(name:'ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

faidx_output
   .collectFile(name:'faidx.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

Shhquis_dot_jl_output
   .collectFile(name:'shhquis.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

simple_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

busco_gfa_output
   .collectFile(name:'busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_no_polish")

simple_busco_fasta_output
   .collectFile(name:'simple.busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_polish" )

jellyfish_output
   .collectFile(name:'jellyfish.log.txt', newLine: true, storeDir:"${params.outdir}/genomescope" )

genomescope2_output
   .collectFile(name:'genomescope2.log.txt', newLine: true, storeDir:"${params.outdir}/genomescope" )

minimap_dot_sh_output
   .collectFile(name:'minimap.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log" )

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

bcftools_version
   .collectFile(name:'bcftools_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
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
