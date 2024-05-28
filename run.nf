#!/usr/bin/env nextflow
//Global Parameters
params.assembly = "an_assembly"

params.readin = "NO_FILE"
params.hicreadf = "NO_FILE_HIC_R1"
params.hicreadr = "NO_FILE_HIC_R2"
params.matreadf = "NO_FILE_MAT_R1"
params.matreadr = "NO_FILE_MAT_R2"
params.patreadf = "NO_FILE_PAT_R1"
params.patreadr = "NO_FILE_PAT_R2"

params.outfasta = "genome.out.fasta"
params.outdir = 'results'
params.threads = '21'
params.lite = false
//Runtype Parameters
params.scaffold = false
params.kmer = 'kmc'
params.polish = false
params.polishtype = 'simple'
params.hapscaffold = false
params.yahs = false
params.patch = false
//HiFIASM Parameters
params.l = 0
params.mode = 'default'
params.ploidy = '2'
//Busco Parameters
params.busco = false
params.buscooffline = false
params.buscodb = "/work/busco"
params.linreage = 'insecta_odb10'
params.buscoevalue = '0.001'
//Shhquis.jl Prameters
params.hclustlinkage = "average"

bam_ch = Channel.fromPath(params.readin)
right_hicfastq_check = Channel.fromPath(params.hicreadr)
left_hicfastq_check = Channel.fromPath(params.hicreadf)

//right_matfastq_check = Channel.fromPath(params.matreadf)
//left_matfastq_check = Channel.fromPath(params.matreadr)
//right_patfastq_check = Channel.fromPath(params.patreadf)
//left_patfastq_check = Channel.fromPath(params.patreadr)

bam_ch.into {
  in_check_ch
  in_Hifi_ch
}

process check_in_file {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = 1

  input:
    file in_file from in_check_ch.flatten()
  output:
    stdout check_in_file_output
  shell:
  '''
   state () { printf "%b\n" "[$(date)]: $*" 2>&1; }
   error () { printf "%b\n" "[$(date)]: $*" >&2; exit 1; }
   warn () { printf "%b\n" "[$(date)]: $*" >&2; }

   touch check_bam.flag.txt
   state "checking file type"
   file=!{in_file}
   stat $file

   if [[ $file == *.bam ]]; then
     state "   ...file type is bam file type"
     samtools flagstat $file
   elif [[ $file == *.fastq.gz || $file == *.fq.gz ]]; then
     state "   ...file type is fastq gz file type"
     state "check if file can be opened, and it starts with @"
     at_check=$(zcat $file | awk '{ print $1; exit }')
     [[ $at_check =~ '@' ]] || error "$file doesn't start with an @";
     state "check if file can be divided by four"
     modulo_four_check=$(zcat $file | grep -v "^=.*"  | wc -l)
     [[ $(( $modulo_four_check % 4 )) -eq 0 ]] || warn "number of lines in $file not divisible by four, continuing anyway"
   elif [[ $file == *.fastq || *.fq ]]; then
     state "   ...file type is fastq file type"
     state "check if file can be opened, and it starts with @"
     at_check=$( head -n 1 $file )
     [[ $at_check =~ '@' ]] || error "$file doesn't start with an @";
     state "check if file can be divided by four"
     modulo_four_check=$(cat $file | grep -v "^=.*" | wc -l)
     [[ $(( $modulo_four_check % 4 )) -eq 0 ]] || warn "number of lines in $file not divisible by four, continuing anyway"
   else
     error "trying to run otb with somthing that does not end with the correct file type"
   fi

   state "check file on $file passed"
   sleep 120;
   exit 0;
  '''
}

process check_fastq {
  label 'shortq'
  cpus = 1

  input:
    file right_fastq from right_hicfastq_check
    file left_fastq from left_hicfastq_check
  output:
    file 'out/right.fastq.gz' into right_fastq_hicstuff, right_fastq_hicstuff_polish, right_yahs, simple_right_yahs, merfin_right_yahs, dv_right_yahs
    file 'out/left.fastq.gz' into left_fastq_hicstuff, left_fastq_hicstuff_polish, left_yahs, simple_left_yahs, merfin_left_yahs, dv_left_yahs
    file 'out/*.fastq.gz' into fasta_in_ch
    stdout check_fastq_output
 when:
    !params.lite
  shell:
  '''
   state () { printf "%b\n" "[$(date)]: $*" 2>&1; }
   error () { printf "%b\n" "[$(date)]: $*" >&2; exit 1; }
   warn () { printf "%b\n" "[$(date)]: $*" >&2; }

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
   [[ !{right_fastq}  =~ ".gz" ]] && first=$(zcat !{right_fastq} | grep -v "^=.*" | wc -l) || first=$( cat !{right_fastq} | grep -v "^=.*" | wc -l)
   [[ !{left_fastq} =~ ".gz" ]] && second=$(zcat !{left_fastq} | grep -v "^=.*"  | wc -l) || second=$( cat !{left_fastq} | grep -v "^=.*"  | wc -l )

   [[ $(( $first % 4 )) -eq 0 ]] || warn "number of lines in !{right_fastq} not divisable by four, continuing anyway"
   [[ $(( $second % 4 )) -eq 0 ]] || warn "number of lines in !{left_fastq} not divisable by four, continuing anyway"

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

process HiFiAdapterFilt {
  label 'shortq'
  container = 'dmolik/pbadapterfilt'
  cpus = params.threads

  input:
    file in_file from in_Hifi_ch.flatten()
  output:
    file '*.filt.fastq' into hifiasm_filt_fastq_ch, filt_fastq_ch, minimap_dv_filt_ch,  minimap_merfin_filt_ch, meryl_filt_ch
    stdout pbadapterfilt_output
  """
    touch pbadapterfilt.flag.txt
    pbadapterfilt.sh ${in_file} -t ${task.cpus}
    gzip -d *.filt.fastq.gz
    echo "finished adapter filtering"
    sleep 120;
    exit 0;
  """
}

process HiFiASM {
  label 'longq'
  container = 'dmolik/hifiasm'
  cpus = params.threads

  input:
    file fasta from hifiasm_filt_fastq_ch.collect()
  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch
    stdout HiFiASM_output

  script:
    if( params.mode == 'phasing' && params.hicreadr != 'NO_FILE_HIC_R1' && params.hicreadr != 'NO_FILE_HIC_R2' )
    """
      touch hifiasm.flag.txt
      hifiasm -l${params.l} -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec --h1 ${params.hicreadf} --h2 ${params.hicreadr} ${fasta} 2>&1
      echo "finished alignment"
      sleep 120;
      exit 0;
    """
    else if( params.mode == 'default')
    """
      touch hifiasm.flag.txt
      hifiasm -l${params.l} -o ${params.assembly} -t ${task.cpus} --write-paf --write-ec ${fasta} 2>&1
      echo "finished alignment"
      sleep 120;
      exit 0;
    """
    else if( params.mode == 'primary')
    """
      touch hifiasm.flag.txt
      hifiasm -l${parms.l} -o ${params.assembly} --primary -t ${task.cpus} --write-paf --write-ec ${fasta} 2>&1
      echo "finished alignment"
      sleep 120;
      exit 0;
    """
    else if( params.mode == 'trio' && params.patreadf != 'NO_FILE_PAT_R1'  && params.patreadr != 'NO_FILE_PAT_R2'  && params.matreadf != 'NO_FILE_MAT_R1'  && params.matreadr != 'NO_FILE_MAT_R2')
    """
      touch hifiasm.flag.txt
      yak count -b37 -t${task.cpus} -o pat.yak <(zcat ${params.patreadf}) <(zcat ${params.patreadr})
      yak count -b37 -t${task.cpus} -o mat.yak <(zcat ${params.matreadf}) <(zcat ${params.matreadr})
      hifiasm -l${parms.l} -o ${params.assembly} --primary -t ${task.cpus} --write-paf --write-ec -1 pat.yak -2 mat.yak${fasta} 2>&1
      echo "finished alignment"
      sleep 120;
      exit 0;
    """
    else
      error "Invalid alignment mode: ${params.mode}"

}

process gfa2fasta {
  label 'shortq'
  publishDir "${params.outdir}/01_hifiasm", mode: 'rellink'
  container = 'pvstodghill/any2fasta'
  cpus 1

  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.p_ctg.gfa.fasta' optional true into gfa2fasta_fasta_res_ch
    file '*[a-z].p_ctg.gfa.fasta' optional true into fasta_unoriented_ch, fasta_genome_ch, fasta_busco_ch, no_polish_yahs_align_genome_ch, fasta_fai_yahs_genome_ch
    file '*hap[12].p_ctg.gfa.fasta' optional true into fasta_hap_ch, simple_fasta_hap_polish_ch, merfin_fasta_hap_polish_ch, dv_fasta_hap_polish_ch, yahs_fasta_hap_polish_ch, yahs_simple_fasta_hap_polish_ch, yahs_merfin_fasta_hap_polish_ch, yahs_dv_fasta_hap_polish_ch
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
  label 'longq'
  publishDir "${params.outdir}/01_hifiasm/busco", mode: 'rellink'
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
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
    exit 0;
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
    exit 0;
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
    exit 0;
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
    exit 0;
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
    exit 0;
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
    exit 0;
  """
  else
  """
    touch busco.flag.txt
  """
}

process ragtag_dot_py {
  label 'longq'
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

  script:

  if ( params.patch )
  """
    touch ragtag.flag.txt
    ragtag.py patch -u --aligner minimap2 -t ${task.cpus} --mm2-params '-x map-hifi -t ${task.cpus} -I 64GB -2 -K 4G -w 21 -k 21' -o ./${params.assembly}_ragtag_ec_patch ${fasta} ${fasta_ec}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
  else
  """
    touch ratag.flag.txt
    echo "ignoring ragtag patch"
    mkdir ${params.assembly}_ragtag_ec_patch
    cd ${params.assembly}_ragtag_ec_patch
    ln -s ../${fasta} ragtag.patch.fasta
    sleep 120;
    exit 0;
  """
}

process faidx {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus 1

  input:
    file genome from fasta_fai_genome_ch
  output:
    file "*.fai" into fai_ch
    stdout faidx_output
  when:
    params.polish
  """
    touch faidx.flag.txt
    samtools faidx -o ${genome}.fai ${genome}
    echo "finished indexing"
    sleep 120;
    exit 0;
  """
}

process yahs_faidx {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus 1

  input:
    file genome from fasta_fai_yahs_genome_ch
  output:
    file "${params.assembly}.yahs.fasta.fai" into yahs_fai_ch
    file "${params.assembly}.yahs.fasta" into yahs_genome_ch
  when:
    params.yahs && params.scaffold
  """
    touch faidx.yahs.flag.txt
    ln -s ${genome} "${params.assembly}.yahs.fasta"
    samtools faidx -o "${params.assembly}.yahs.fasta.fai" "${params.assembly}.yahs.fasta"
    echo "finished indexing"
    sleep 120;
    exit 0;
  """
}

process hicstuff {
  label 'longq'
  publishDir "${params.outdir}/02_hicstuff", mode: 'rellink'
  container = 'koszullab/hicstuff'
  cpus = params.threads

  input:
    file genome from fasta_genome_ch
    file left from left_fastq_hicstuff
    file right from right_fastq_hicstuff
  output:
    file 'hicstuff_out/abs_fragments_contacts_weighted.bg2' into abs_ch
    file 'hicstuff_out/info_contigs.txt' into contigs_ch
    file 'hicstuff_out/fragments_list.txt'
    file 'hicstuff_out/plots/frags_hist.pdf'
    stdout hicstuff_output
  when:
    params.scaffold
  """
    touch hicstuff.flag.txt
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${left} ${right}
    echo "finished fragment calculations"
    sleep 120;
    exit 0;
  """
}

/*process hicstuff_polish {
  label 'longq'
  publishDir "${params.outdir}/02_hicstuff", mode: 'rellink'
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
  when:
    params.polish && params.scaffold
  """
    touch hicstuff_for_polished.flag.txt
    hicstuff pipeline -t ${task.cpus} -a minimap2 --no-cleanup -e 10000000 --force --out hicstuff_out --duplicates --matfmt=bg2 --plot -g ${genome} ${left} ${right}
    mv hicstuff_out/fragments_list.txt hicstuff_out/polish_fragments_list.txt
    mv hicstuff_out/plots/frags_hist.pdf hicstuff_out/plots/polish_frags_hist.pdf
    echo "finished fragment calculations"
    sleep 120;
    exit 0;
  """
}*/

process Shhquis_dot_jl {
  label 'mediumq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'dmolik/shhquis'
  cpus params.threads

  input:
    file abs from abs_ch
    file contig from contigs_ch
    file genome from fasta_sshquis_genome_ch
    file fai from fai_ch
  output:
    file "${params.outfasta}" into shhquis_fasta_res_ch, polish_haps_genome_ch, shhquis_simple_ch, shhquis_merfin_ch, shhquis_dv_fai_ch, shhquis_dv_ch, shhquis_bcftools_dv_ch, shhquis_bcftools_merfin_ch, shhquis_dv_minimap_ch, shquis_minimap_merfin_ch, shhquis_mpileup_ch
    file "${params.outfasta}"
    stdout Shhquis_dot_jl_output
  when:
    params.polish
  """
    touch shhquis.flag.txt
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage ${params.hclustlinkage} --threads ${task.cpus}
    echo "finished reorientation"
    sleep 120;
    exit 0;
  """
}

process K_mer_counting {
  label 'mediumq'
  container = 'dmolik/k-mer-counting-tools'
  cpus = params.threads

  input:
    file filt_reads from filt_fastq_ch.collect()
  output:
    file '*.histo' into histo_ch
    file 'version.txt' into jellyfish_ver_ch
    stdout jellyfish_output
  script:

  if( params.kmer == 'jellyfish' )
  """
    touch jellyfish.flag.txt
    jellyfish count -C -m 21 -s 1000000000 -t ${task.cpus} -o reads.jf ${filt_reads}
    jellyfish histo -t ${task.cpus} reads.jf > ${params.assembly}.histo
    jellyfish cite > version.txt
    sleep 120;
    exit 0;
  """
  else if( params.kmer == 'kmc' )
  """
    mkdir tmp
    ls ${filt_reads} > FILES.lst
    kmc -v -k21 -t${task.cpus} -ci1 -cs10000 @FILES.lst reads tmp/
    kmc_tools transform reads histogram ${params.assembly}.histo -cx10000
    kmc | head -n 1 > version.txt
  """
  else
    error "Invalid k-mer tool: ${params.kmer}"
}

process genomescope2 {
  label 'mediumq'
  publishDir "${params.outdir}/00_ordination/genomescope", mode: 'rellink'
  container = 'dmolik/genomescope2'
  cpus = params.threads

  input:
    file histo from histo_ch
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


process simple_fcs_adaptor {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'ncbi/fcs-adaptor'
  cpus = 1

  input:
   file genome from shhquis_simple_ch
  output:
   file 'cleaned_sequences/*' into simple_fcs_adaptor_ch
   file '*'
   stdout simple_fcs_adaptor_output
  when:
   params.polishtype == "simple"
  """
    touch dv_fcs_adaptor.flag.txt
    /app/fcs/bin/av_screen_x -o . --euk ${genome}
    echo "finished simple fcs adaptor"
    sleep 120;
    exit 0;
  """
}

process simple_polish {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = 1

  input:
    file genome from simple_fcs_adaptor_ch
  output:
    file "${params.assembly}.polished.genome.fasta" into simple_polished_genome_ch, simple_polished_genome_busco_ch, yahs_simple_genome_ch, yahs_simple_align_genome_ch
    file "${params.assembly}.polished.genome.fasta.fai" into yahs_simple_fai_ch
  when:
    params.polishtype == "simple"
  """
    touch simple_polish.flag.txt
    ln -s ${genome} ${params.assembly}.polished.genome.fasta
    samtools faidx -o ${params.assembly}.polished.genome.fasta.fai ${params.assembly}.polished.genome.fasta
    echo "finished softlink"
    sleep 120;
    exit 0;
  """
}

process minimap_for_merfin {
  label 'mediumq'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file filt_reads from minimap_merfin_filt_ch
    file genome from shquis_minimap_merfin_ch
  output:
    file "mapped.sam" into sam_for_merfin_ch
    stdout minimap_dot_sh_output
  when:
    params.polishtype == "merfin"
  """
    touch minimap.flag.sh
    minimap2 -a ${genome} ${filt_reads} > mapped.sam
    echo "finished minimap"
    sleep 120;
    exit 0;
  """
}

process samtools_mpileup_merfin {
  label 'mediumq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = params.threads

  input:
    file sam_file from sam_for_merfin_ch
    file genome from shhquis_mpileup_ch
  output:
    file 'out.mpileup' into bcf_for_merfin_ch
    stdout samtools_mpileup_output
  when:
    params.polishtype == "merfin"
  shell:
  '''
    touch samtools.mpileup.flag.txt
    samtools sort -@ !{task.cpus} -o aln.bam !{sam_file}
    samtools index -@ !{task.cpus} aln.bam
    samtools view -H aln.bam | grep '@SQ' | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P !{task.cpus} sh -c "samtools mpileup -BQ0 -d 100000 -uf !{genome} -r {} aln.bam > tmp.{}.mpileup"
    cat tmp.*.mpileup > out.mpileup
    echo "finished mpileup"
    sleep 120;
    exit 0;
  '''
}

process bcftools_refmt {
  label 'shortq'
  container = 'mgibio/bcftools:1.9'
  cpus = params.threads

  input:
    file mpileup from bcf_for_merfin_ch
  output:
    file "final.reshaped.vcf.gz" into vcf_for_merfin_ch
    stdout bcftools_refmt_output
  when:
    params.polishtype == "merfin"
  """
    touch bcftools.qual.flag.txt
    echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">' > merfin_header.vcf
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> merfin_header.vcf
    echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">' >> merfin_header.vcf
    cat ${mpileup} | bcftools call --threads ${task.cpus} -mv > var.raw.vcf
    bcftools filter --threads ${task.cpus} -s LowQual -e '%QUAL<20 || DP>100' var.raw.vcf  > var.flt.vcf
    grep -v "#" var.flt.vcf | sed 's/,/;/g' > var.temp.reshaped.vcf
    bcftools view --threads ${task.cpus} -h var.flt.vcf > var.temp.reshaped.header.vcf
    cat var.temp.reshaped.header.vcf var.temp.reshaped.vcf > var.temp.reshaped.combined.vcf
    rm var.temp.reshaped.header.vcf var.temp.reshaped.vcf
    bcftools annotate --threads ${task.cpus} -h merfin_header.vcf var.temp.reshaped.combined.vcf > var.temp.reshaped.vcf
    bcftools view --threads ${task.cpus} -h var.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > var.reshaped.vcf
    rm var.temp.reshaped.vcf
    bcftools view --threads ${task.cpus} -H var.temp.reshaped.combined.vcf | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,".\tGT:DP\t1/1:",\$8);print \$0}' >> var.reshaped.vcf
    bcftools view --threads ${task.cpus} var.reshaped.vcf -Oz > final.reshaped.vcf.gz
    rm var.reshaped.vcf
    rm var.temp.reshaped.combined.vcf
    echo "finished bcftools reformat"
    sleep 120;
    exit 0;
  """
}

process merfin {
  label 'longq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'dmolik/merfin'
  cpus = params.threads

  input:
    file genome from shhquis_merfin_ch
    file kcov_file from kcov_ch
    file lookup_table from lookup_table_ch
    file filt_reads from meryl_filt_ch
    file vcf_file from vcf_for_merfin_ch
  output:
    file 'merfin.polish.vcf' into merfin_vcf_ch
    stdout merfin_output
  when:
    params.polishtype == "merfin"

  shell:
    '''
      echo "warning merfin is experimental"
      touch merfin.flag.txt
      meryl count k=21 !{filt_reads} output reads.meryl
      meryl greater-than 1 reads.meryl output reads.gt1.meryl
      merfin -polish -threads !{task.cpus} -sequence !{genome} -peak $( cat !{kcov_file} | sed 's/e.*//g' ) -prob !{lookup_table} -readmers reads.gt1.meryl -vcf !{vcf_file} -output merfin
      echo "finished merfin"
      sleep 120;
      exit 0;
    '''
}

process minimap_for_deep_variant {
  label 'mediumq'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file filt_reads from minimap_dv_filt_ch
    file genome from shhquis_dv_minimap_ch
  output:
    file "mapped.sam" into sam_for_dv_ch
    stdout minimap_dv_output
  when:
      params.polishtype == "dv"
    """
      touch minimap.dv.flag.sh
      minimap2 -t ${task.cpus} -a ${genome} ${filt_reads} > mapped.sam
      echo "finished minimap"
      sleep 120;
      exit 0;
   """
}


process samtools_index_for_deep_variant {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = params.threads

  input:
    file sam from sam_for_dv_ch
    file genome from shhquis_dv_fai_ch
  output:
    file 'mapped.sort.bam' into bam_dv_index_ch
    file '*.bai' into bai_dv_index_ch
    file '*.fai' into fai_dv_index_ch
  when:
    params.polishtype == "dv"
  """
    touch samtools.index.flag.txt
    samtools view -S -b ${sam} > mapped.bam
    samtools sort -@ ${task.cpus} mapped.bam -o mapped.sort.bam
    samtools index -@ ${task.cpus} mapped.sort.bam
    samtools faidx ${genome}
    echo "finished indexing"
    sleep 120;
    exit 0;
  """
}

process deep_variant {
  label 'longq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'google/deepvariant'
  cpus = params.threads

  input:
    file genome from shhquis_dv_ch
    file genome_fai from fai_dv_index_ch
    file bam_read from bam_dv_index_ch
    file bai_read from bai_dv_index_ch
  output:
    file 'google_dv.vcf' into dv_vcf_ch
    file '*'
    stdout deep_variant_output
  when:
    params.polishtype == "dv"
  """
    echo "warning deep variant is experimental"
    touch deep_variant.flag.txt
    /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=${genome} --reads=${bam_read} --output_vcf=google_dv.vcf --output_gvcf=google_dv.gvcf --num_shards=${task.cpus}
    echo "finished deep variant"
    sleep 120;
    exit 0;
  """
}

process dv_bcftools {
  label 'mediumq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'mgibio/bcftools:1.9'
  cpus = params.threads

  input:
    file genome from shhquis_bcftools_dv_ch
    file vcf from dv_vcf_ch
  output:
    file "${params.assembly}.vcf_polished_assembly.fasta" into dv_fcs_adaptor_ch
    file "${params.assembly}.vcf_polished_assembly.fasta"
    stdout dv_bcftools_output
  when:
    params.polishtype == "dv"
  """
    touch dv.bcftools.flag.txt
    bcftools view --threads ${task.cpus} -Oz ${vcf} > ${vcf}.gz
    bcftools index --threads ${task.cpus} ${vcf}.gz
    bcftools consensus ${vcf}.gz -f ${genome} -H 1 > ${params.assembly}.vcf_polished_assembly.fasta
    echo "finished bcftools from deep variant"
    sleep 120;
    exit 0;
  """
}

process dv_fcs_adaptor {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'ncbi/fcs-adaptor'
  cpus = 1

  input:
   file genome from dv_fcs_adaptor_ch
  output:
   file 'cleaned_sequences/*.fa' into dv_vcf_polished_genome_ch, dv_vcf_res_ch, dv_vcf_polished_busco_genome_ch, yahs_dv_genome_ch, yahs_dv_align_genome_ch
   file '*'
   stdout dv_fcs_adaptor_output
  when:
   params.polishtype == "dv"
  """
    touch dv_fcs_adaptor.flag.txt
    /app/fcs/bin/av_screen_x -o . --euk ${genome}
    gzip -d cleaned_sequences/*.fa.gz
    echo "finished dv fcs adaptor"
    sleep 120;
    exit 0;
  """
}

process merfin_bcftools {
  label 'mediumq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'mgibio/bcftools:1.9'
  cpus = params.threads

  input:
    file genome from shhquis_bcftools_merfin_ch
    file vcf from merfin_vcf_ch
  output:
    file "${params.assembly}.vcf_polished_assembly.fasta" into merfin_fcs_adaptor_ch
    file "${params.assembly}.vcf_polished_assembly.fasta"
    stdout merfin_bcftools_output
  when:
    params.polishtype == "merfin"
  """
    touch merfin.bcftools.flag.txt
    bcftools view --threads ${task.cpus} -Oz ${vcf} > ${vcf}.gz
    bcftools index --threads ${task.cpus} ${vcf}.gz
    bcftools consensus ${vcf}.gz -f ${genome} > ${params.assembly}.vcf_polished_assembly.fasta
    echo "finished bcftools from merfin"
    sleep 120;
    exit 0;
  """
}

process merfin_fcs_adaptor {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'ncbi/fcs-adaptor'
  cpus = 1

  input:
   file genome from merfin_fcs_adaptor_ch
  output:
   file 'cleaned_sequences/*.fa' into merfin_vcf_polished_genome_ch, merfin_vcf_res_ch, merfin_vcf_polished_busco_genome_ch, yahs_merfin_genome_ch, yahs_merfin_align_genome_ch
   file '*'
   stdout merfin_fcs_adaptor_output
  when:
   params.polishtype == "merfin"
  """
    touch merfin_fcs_adaptor.flag.txt
    /app/fcs/bin/av_screen_x -o . --euk ${genome}
    gzip -d cleaned_sequences/*.fa.gz
    echo "finished merfin fcs adaptor"
    sleep 120;
    exit 0;
  """
}

process bwa_for_yahs {
  label 'mediumq'
  container = 'staphb/bwa'
  cpus = params.threads

  input:
    file left_reads from left_yahs
    file right_reads from right_yahs
    file genome from no_polish_yahs_align_genome_ch
  output:
    file "mapped.sam" into yahs_sam_ch
    stdout bwa_for_yahs_output
  when:
      params.yahs
    """
      touch bwa.yahs.flag.sh
      bwa index ${genome}
      bwa mem -5SP -t ${task.cpus} ${genome} ${left_reads} ${right_reads} > mapped.sam
      echo "finished bwa"
      sleep 120;
      exit 0;
   """
}

process bwa_for_simple_yahs {
  label 'mediumq'
  container = 'staphb/bwa'
  cpus = params.threads

  input:
    file left_reads from simple_left_yahs
    file right_reads from simple_right_yahs
    file genome from yahs_simple_align_genome_ch
  output:
    file "mapped.sam" into yahs_simple_sam_ch
    stdout bwa_for_yahs_simple_output
  when:
      params.yahs
    """
      touch minimap.yahs.simple.flag.sh
      bwa index ${genome}
      bwa mem -5SP -t ${task.cpus} ${genome} ${left_reads} ${right_reads} > mapped.sam
      echo "finished bwa"
      sleep 120;
      exit 0;
   """
}

process bwa_for_merfin_yahs {
  label 'mediumq'
  container = 'staphb/bwa'
  cpus = params.threads

  input:
    file left_reads from merfin_left_yahs
    file right_reads from merfin_right_yahs
    file genome from yahs_merfin_align_genome_ch
  output:
    file "mapped.sam" into yahs_merfin_sam_ch
    stdout bwa_for_yahs_merfin_output
  when:
      params.yahs
    """
      touch bwa.yahs.merfin.flag.sh
      bwa index ${genome}
      bwa mem -5SP -t ${task.cpus} ${genome} ${left_reads} ${right_reads} > mapped.sam
      echo "finished bwa"
      sleep 120;
      exit 0;
   """
}

process bwa_for_dv_yahs {
  label 'mediumq'
  container = 'staphb/bwa'
  cpus = params.threads

  input:
    file left_reads from dv_left_yahs
    file right_reads from dv_right_yahs
    file genome from yahs_dv_align_genome_ch
  output:
    file "mapped.sam" into yahs_dv_sam_ch
    stdout bwa_for_yahs_dv_output
  when:
      params.yahs
    """
      touch bwa.yahs.merfin.flag.sh
      bwa index ${genome}
      bwa mem -5SP -t ${task.cpus} ${genome} ${left_reads} ${right_reads} > mapped.sam
      echo "finished bwa"
      sleep 120;
      exit 0;
   """
}

process bam_sort_for_yahs {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = params.threads

  input:
    file sam_file from yahs_sam_ch
  output:
    file 'aln.bam' into bam_for_yahs_ch
    stdout bam_sort_for_yahs_output
  when:
     params.yahs
  shell:
  '''
    touch bam.sort.yahs.flag.txt
    samtools view -S -h -F 2316 !{sam_file} | samblaster | samtools sort -n -@ !{task.cpus} -o aln.bam
    echo "finished sort"
    sleep 120;
    exit 0;
  '''
}

process bam_sort_for_simple_yahs {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = params.threads

  input:
    file sam_file from yahs_simple_sam_ch
  output:
    file 'aln.bam' into bam_for_simple_yahs_ch
    stdout bam_sort_for_simple_yahs_output
  when:
     params.yahs
  shell:
  '''
    touch bam.sort.simple.yahs.flag.txt
    samtools view -S -h -F 2316 !{sam_file} | samblaster | samtools sort -n -@ !{task.cpus} -o aln.bam
    echo "finished sort"
    sleep 120;
    exit 0;
  '''
}

process bam_sort_for_merfin_yahs {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = params.threads

  input:
    file sam_file from yahs_merfin_sam_ch
  output:
    file 'aln.bam' into bam_for_merfin_yahs_ch
    stdout bam_sort_for_merfin_yahs_output
  when:
     params.yahs
  shell:
  '''
    touch bam.sort.merfin.yahs.flag.txt
    samtools view -S -h -F 2316 !{sam_file} | samblaster | samtools sort -n -@ !{task.cpus} -o aln.bam
    echo "finished sort"
    sleep 120;
    exit 0;
  '''
}

process bam_sort_for_dv_yahs {
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus = params.threads

  input:
    file sam_file from yahs_dv_sam_ch
  output:
    file 'aln.bam' into bam_for_dv_yahs_ch
    stdout bam_sort_for_dv_yahs_output
  when:
     params.yahs
  shell:
  '''
    touch bam.sort.dv.yahs.flag.txt
    samtools view -S -h -F 2316 !{sam_file} | samblaster | samtools sort -n -@ !{task.cpus} -o aln.bam
    echo "finished sort"
    sleep 120;
    exit 0;
  '''
}

process yahs {
  label 'mediumq'
  publishDir "${params.outdir}/04_yahs", mode: 'rellink'
  container = 'dmolik/yahs'

  input:
    file input_bam from bam_for_yahs_ch
    file input_genome from yahs_genome_ch
    file input_fai from yahs_fai_ch
  output:
    file "yahs.no_polish*"
    file "yahs.no_polish_scaffolds_final.fa" into yahs_no_polish_stats_ch, yahs_no_polish_haps_genome_ch, yahs_no_polish_busco_ch
    file "*JBAT.txt" into yahs_JBAT_txt_ch
    file "yahs.no_polish.JBAT/tmp_juicer_pre_JBAT.log" into yahs_JBAT_ch
    stdout yahs_output
  when:
     params.yahs
  """
    touch yahs.flag.txt
    yahs -o yahs.no_polish --no-contig-ec ${input_genome} ${input_bam}
    mkdir yahs.no_polish.JBAT
    juicer_pre -a -o yahs.no_polish.JBAT yahs.no_polish.bin yahs.no_polish_scaffolds_final.agp ${input_fai} 2> yahs.no_polish.JBAT/tmp_juicer_pre_JBAT.log
    sleep 120;
    exit 0;
  """
}

process simple_yahs {
  label 'mediumq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/yahs'

  input:
    file input_bam from bam_for_simple_yahs_ch
    file input_genome from yahs_simple_genome_ch
    file input_fai from yahs_simple_fai_ch
  output:
    file "yahs.simple*"
    file "yahs.simple_scaffolds_final.fa" into yahs_simple_polish_stats_ch, yahs_simple_polish_haps_genome_ch, yahs_simple_polish_busco_ch
    file "*JBAT.txt" into yahs_simple_JBAT_txt_ch
    file "yahs.simple.JBAT/tmp_juicer_pre_JBAT.log" into yahs_simple_JBAT_ch
    stdout yahs_simple_output
  when:
     params.yahs
  """
    touch simple_yahs.flag.txt
    yahs -o yahs.simple --no-contig-ec ${input_genome} ${input_bam}
    mkdir yahs.simple.JBAT
    juicer_pre -a -o yahs.simple.JBAT yahs.simple.bin yahs.simple_scaffolds_final.agp ${input_fai} 2> yahs.simple.JBAT/tmp_juicer_pre_JBAT.log
    sleep 120;
    exit 0;
  """
}

process merfin_yahs_faidx {
  label 'mediumq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus 1

  input:
    file genome from yahs_merfin_genome_ch
  output:
    file "${params.assembly}.yahs.fasta.fai" into yahs_merfin_fai_genome_fai_ch
    file "${params.assembly}.yahs.fasta" into yahs_merfin_fai_genome_ch
  when:
    params.yahs
  """
    touch faidx.yahs.flag.txt
    ln -s ${genome} "${params.assembly}.yahs.fasta"
    samtools faidx -o "${params.assembly}.yahs.fasta.fai" "${params.assembly}.yahs.fasta"
    echo "finished indexing"
    sleep 120;
    exit 0;
  """
}

process dv_yahs_faidx {
  label 'mediumq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
  cpus 1

  input:
    file genome from yahs_dv_genome_ch
  output:
    file "${params.assembly}.yahs.fasta.fai" into yahs_dv_fai_genome_fai_ch
    file "${params.assembly}.yahs.fasta" into yahs_dv_fai_genome_ch
  when:
    params.yahs
  """
    touch faidx.yahs.flag.txt
    ln -s ${genome} "${params.assembly}.yahs.fasta"
    samtools faidx -o "${params.assembly}.yahs.fasta.fai" "${params.assembly}.yahs.fasta"
    echo "finished indexing"
    sleep 120;
    exit 0;
  """
}

process merfin_yahs {
  label 'mediumq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/yahs'

  input:
    file input_bam from bam_for_merfin_yahs_ch
    file input_genome from yahs_merfin_fai_genome_ch
    file input_fai from yahs_merfin_fai_genome_fai_ch
  output:
    file "yahs.merfin*"
    file "yahs.merfin_scaffolds_final.fa" into yahs_merfin_polish_stats_ch, yahs_merfin_polish_haps_genome_ch, yahs_merfin_polish_busco_ch
    file "*JBAT.txt" into yahs_merfin_JBAT_txt_ch
    file "yahs.merfin.JBAT/tmp_juicer_pre_JBAT.log" into yahs_merfin_JBAT_ch
    stdout yahs_merfin_output
  when:
     params.yahs
  """
    touch merfin_yahs.flag.txt
    yahs -o yahs.merfin --no-contig-ec ${input_genome} ${input_bam}
    mkdir yahs.merfin.JBAT
    juicer_pre -a -o yahs.merfin.JBAT yahs.merfin.bin yahs.merfin_scaffolds_final.agp ${input_fai} 2> yahs.merfin.JBAT/tmp_juicer_pre_JBAT.log
    sleep 120;
    exit 0;
  """
}


process dv_yahs {
  label 'mediumq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/yahs'

  input:
    file input_bam from bam_for_dv_yahs_ch
    file input_genome from yahs_dv_fai_genome_ch
    file input_fai from yahs_dv_fai_genome_fai_ch
  output:
    file "yahs.dv*"
    file "yahs.dv_scaffolds_final.fa" into yahs_dv_polish_stats_ch, yahs_dv_polish_haps_genome_ch, yahs_dv_polish_busco_ch
    file "*JBAT.txt" into yahs_dv_JBAT_txt_ch
    file "yahs.dv.JBAT/tmp_juicer_pre_JBAT.log" into yahs_dv_JBAT_ch
    stdout yahs_dv_output
  when:
     params.yahs
  """
    touch dv_yahs.flag.txt
    yahs -o yahs.dv --no-contig-ec ${input_genome} ${input_bam}
    mkdir yahs.dv.JBAT
    juicer_pre -a -o yahs.dv.JBAT yahs.dv.bin yahs.dv_scaffolds_final.agp ${input_fai} 2> yahs.dv.JBAT/tmp_juicer_pre_JBAT.log
    sleep 120;
    exit 0;
  """
}

process juicer_tools_pre_yahs {
  label 'shortq'
  publishDir "${params.outdir}/04_yahs", mode: 'rellink'
  container = 'dmolik/juicer-tools'
  cpus = params.threads

  input:
    file yahs_JBAT from yahs_JBAT_ch
    file yahs_JBAT_txt from yahs_JBAT_txt_ch
  output:
    file "*_JBAT.hic"
    stdout juicer_tools_pre_yahs_output

  """
    touch juicer_tools_pre_yahs.flag.txt
    java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre ${yahs_JBAT_txt} ${yahs_JBAT_txt}.hic.part <(cat ${yahs_JBAT}  | grep PRE_C_SIZE | awk '{print \$2" "\$3}') && (mv ${yahs_JBAT_txt}.hic.part no_polish_JBAT.hic)
    exit 0;
  """
}

process juicer_tools_pre_yahs_simple {
  label 'shortq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/juicer-tools'
  cpus = params.threads

  input:
    file yahs_JBAT from yahs_simple_JBAT_ch
    file yahs_JBAT_txt from yahs_simple_JBAT_txt_ch
  output:
    file "*_JBAT.hic"
    stdout juicer_tools_pre_yahs_simple_output

  """
    touch juicer_tools_pre_yahs_simple.flag.txt
    java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre ${yahs_JBAT_txt} ${yahs_JBAT_txt}.hic.part <(cat ${yahs_JBAT}  | grep PRE_C_SIZE | awk '{print \$2" "\$3}') && (mv ${yahs_JBAT_txt}.hic.part simple_JBAT.hic)
    exit 0;
  """
}

process juicer_tools_pre_yahs_merfin {
  label 'shortq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/juicer-tools'
  cpus = params.threads

  input:
    file yahs_JBAT from yahs_merfin_JBAT_ch
    file yahs_JBAT_txt from yahs_merfin_JBAT_txt_ch
  output:
    file "*_JBAT.hic"
    stdout juicer_tools_pre_yahs_merfin_output

  """
    touch juicer_tools_pre_yahs_merfin.flag.txt
    java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre ${yahs_JBAT_txt} ${yahs_JBAT_txt}.hic.part <(cat ${yahs_JBAT}  | grep PRE_C_SIZE | awk '{print \$2" "\$3}') && (mv ${yahs_JBAT_txt}.hic.part merfin_JBAT.hic)
    exit 0;
  """
}

process juicer_tools_pre_yahs_dv {
  label 'shortq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/juicer-tools'
  cpus = params.threads

  input:
    file yahs_JBAT from yahs_dv_JBAT_ch
    file yahs_JBAT_txt from yahs_dv_JBAT_txt_ch
  output:
    file "*_JBAT.hic"
    stdout juicer_tools_pre_yahs_dv_output

  """
    touch juicer_tools_pre_yahs_dv.flag.txt
    java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre ${yahs_JBAT_txt} ${yahs_JBAT_txt}.hic.part <(cat ${yahs_JBAT}  | grep PRE_C_SIZE | awk '{print \$2" "\$3}') && (mv ${yahs_JBAT_txt}.hic.part dv_JBAT.hic)
    exit 0;
  """
}

process ragtag_dot_py_hap_simple_polish {
  label 'longq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from simple_fasta_hap_polish_ch.flatten()
    file fasta_genome from simple_polished_genome_ch
  output:
    file "polished*fasta"
    file "polished*fasta" into simple_hap_patch_res_ch
    stdout simple_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.polishtype == "simple"
  """
    touch ragtag.hap.flag.txt
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_hap_merfin_polish {
  label 'longq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from merfin_fasta_hap_polish_ch.flatten()
    file fasta_genome from merfin_vcf_polished_genome_ch
  output:
    file "polished*fasta"
    file "polished*fasta" into merfin_hap_patch_res_ch
    stdout merfin_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.polishtype == "merfin"
  """
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_hap_deep_variant_polish {
  label 'longq'
  publishDir "${params.outdir}/03_polish", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from dv_fasta_hap_polish_ch.flatten()
    file fasta_genome from dv_vcf_polished_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/polished*fasta"
    file "${params.assembly}_ragtag_scaffold/polished*fasta" into dv_hap_patch_res_ch
    stdout dv_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.polishtype == "dv"
  """
    touch ragtag.hap.flag.txt
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_yahs {
  label 'longq'
  publishDir "${params.outdir}/04_yahs", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from yahs_fasta_hap_polish_ch.flatten()
    file fasta_genome from yahs_no_polish_haps_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/polished*fasta"
    file "${params.assembly}_ragtag_scaffold/polished*fasta" into yahs_hap_patch_res_ch
    stdout yahs_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.yahs
  """
    touch ragtag.hap.yahs.flag.txt
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_simple_yahs {
  label 'longq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from yahs_simple_fasta_hap_polish_ch.flatten()
    file fasta_genome from yahs_simple_polish_haps_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/polished*fasta"
    file "${params.assembly}_ragtag_scaffold/polished*fasta" into yahs_simple_hap_patch_res_ch
    stdout yahs_simple_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.polishtype == "simple"
  """
    touch ragtag.hap.yahs.simple.flag.txt
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_merfin_yahs {
  label 'longq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from yahs_merfin_fasta_hap_polish_ch.flatten()
    file fasta_genome from yahs_merfin_polish_haps_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/polished*fasta"
    file "${params.assembly}_ragtag_scaffold/polished*fasta" into yahs_merfin_hap_patch_res_ch
    stdout yahs_merfin_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.polishtype == "merfin"
  """
    touch ragtag.hap.yahs.simple.flag.txt
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process ragtag_dot_py_dv_yahs {
  label 'longq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from yahs_dv_fasta_hap_polish_ch.flatten()
    file fasta_genome from yahs_dv_polish_haps_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/polished*fasta"
    file "${params.assembly}_ragtag_scaffold/polished*fasta" into yahs_dv_hap_patch_res_ch
    stdout yahs_dv_ragtag_dot_py_hap_output
  when:
    params.hapscaffold && params.polishtype == "dv"
  """
    touch ragtag.hap.yahs.simple.flag.txt
    ragtag.py scaffold --aligner minimap2 -t ${task.cpus} --mm2-params '-x asm10 -t ${task.cpus} -I 8GB -2 -K 2G' -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    sleep 120;
    exit 0;
  """
}

process simple_busco_fasta {
  label 'longq'
  publishDir "${params.outdir}/03_polish/busco", mode: 'rellink'
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
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process merfin_busco_fasta {
  label 'longq'
  publishDir "${params.outdir}/03_polish/busco", mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  input:
    file fasta from merfin_vcf_polished_busco_genome_ch
  output:
    file '*'
    stdout merfin_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process dv_busco_fasta {
  label 'longq'
  publishDir "${params.outdir}/03_polish/busco", mode: 'rellink'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  cpus = params.threads

  input:
    file fasta from dv_vcf_polished_busco_genome_ch
  output:
    file '*'
    stdout dv_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process yahs_busco_fasta {
  label 'longq'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  publishDir "${params.outdir}/04_yahs/busco", mode: 'rellink'
  input:
    file fasta from yahs_no_polish_busco_ch
  output:
    file '*'
    stdout yahs_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process yahs_simple_busco_fasta {
  label 'longq'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  publishDir "${params.outdir}/05_yahs_on_polish/busco", mode: 'rellink'
  input:
    file fasta from yahs_simple_polish_busco_ch
  output:
    file '*'
    stdout yahs_simple_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process yahs_merfin_busco_fasta {
  label 'longq'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  publishDir "${params.outdir}/05_yahs_on_polish/busco", mode: 'rellink'
  input:
    file fasta from yahs_merfin_polish_busco_ch
  output:
    file '*'
    stdout yahs_merfin_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process yahs_merfin_dv_fasta {
  label 'longq'
  container = 'ezlabgva/busco:v5.2.2_cv1'
  publishDir "${params.outdir}/05_yahs_on_polish/busco", mode: 'rellink'
  input:
    file fasta from yahs_dv_polish_busco_ch
  output:
    file '*'
    stdout yahs_dv_busco_fasta_output
  when:
    params.busco

  script:

  if( params.linreage == 'auto-lineage' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage
  """
  else if( params.linreage == 'auto-lineage-prok' && params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-prok
  """
  else if( params.linreage == 'auto-lineage-euk'&& params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} --auto-lineage-euk
  """
  else if( params.buscooffline == false)
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage}
  """
  else if( params.buscooffline == true && params.buscodb == 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$baseDir/work/busco"
  """
  else if( params.buscooffline == true && params.buscodb != 'work/busco')
  """
    touch busco.flag.txt
    busco -q --long -e ${params.buscoevalue} -i ${fasta} -o "${params.assembly}_${fasta}_busco" -m genome -c ${task.cpus} -l ${params.linreage} --offline --download_path "$params.buscodb"
  """
  else
  """
    touch busco.flag.txt
  """
}

process fasta_in_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/00_ordination", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from  fasta_in_ch.flatten()
  output:
    file '*.stats'
  """
    touch any2fasta_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process gfa2fasta_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/01_hifiasm", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from  gfa2fasta_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch any2fasta_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process ragtag_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/02_hicstuff", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from ragtag_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch ragtag_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process shhquis_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from shhquis_fasta_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch shhquis_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process merfin_vcf_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from merfin_vcf_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch vcf_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process dv_vcf_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from dv_vcf_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch vcf_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}


process simple_hap_patch_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from simple_hap_patch_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch hap_patch_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}


process merfin_hap_patch_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from merfin_hap_patch_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch hap_patch_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process dv_hap_patch_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/03_polish", mode: 'copy'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from dv_hap_patch_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch hap_patch_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process yahs_no_polish_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/04_yahs", mode: 'rellink'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from yahs_no_polish_stats_ch
  output:
    file '*.stats'
  """
    touch yahs_no_polish_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process yahs_simple_polish_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from yahs_simple_polish_stats_ch
  output:
    file '*.stats'
  """
    touch yahs_simple_polish_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process yahs_merfin_polish_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from yahs_merfin_polish_stats_ch
  output:
    file '*.stats'
  """
    touch yahs_merfin_polish_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process yahs_dv_polish_stats_dot_sh {
  label 'shortq'
  publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
  container = 'dmolik/gfastats'
  cpus 1

  input:
    file fasta from yahs_dv_polish_stats_ch
  output:
    file '*.stats'
  """
    touch yahs_dv_polish_stats.flag.txt
    gfastats ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process yahs_hap_patch_stats_dot_sh {
   label 'shortq'
   publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
   container = 'dmolik/gfastats'
   cpus 1

   input:
     file fasta from yahs_hap_patch_res_ch
   output:
     file '*.stats'
   """
     touch yahs_hap_no_polish_stats.flag.txt
     gfastats ${fasta} > ${fasta}.stats
     echo "finished stats"
     exit 0;
   """
}

process yahs_hap_patch_simple_polish_stats_dot_sh {
   label 'shortq'
   publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
   container = 'dmolik/gfastats'
   cpus 1

   input:
     file fasta from yahs_simple_hap_patch_res_ch
   output:
     file '*.stats'
   """
     touch yahs_hap_simple_polish_stats.flag.txt
     gfastats ${fasta} > ${fasta}.stats
     echo "finished stats"
     exit 0;
   """
}

process yahs_hap_patch_merfin_polish_stats_dot_sh {
   label 'shortq'
   publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
   container = 'dmolik/gfastats'
   cpus 1

   input:
     file fasta from yahs_merfin_hap_patch_res_ch
   output:
     file '*.stats'
   """
     touch yahs_hap_merfin_polish_stats.flag.txt
     gfastats ${fasta} > ${fasta}.stats
     echo "finished stats"
     exit 0;
   """
}

process yahs_hap_patch_dv_polish_stats_dot_sh {
   label 'shortq'
   publishDir "${params.outdir}/05_yahs_on_polish", mode: 'rellink'
   container = 'dmolik/gfastats'
   cpus 1

   input:
     file fasta from yahs_dv_hap_patch_res_ch
   output:
     file '*.stats'
   """
     touch yahs_hap_dv_polish_stats.flag.txt
     gfastats ${fasta} > ${fasta}.stats
     echo "finished stats"
     exit 0;
   """
}

process HiFiASM_Version {
  label 'shortq'
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
  label 'shortq'
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
  label 'shortq'
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
  label 'shortq'
  container = 'mgibio/alignment_helper-cwl:2.2.1'
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

process YaHS_Version {
  label 'shortq'
  container = 'dmolik/yahs'
  cpus 1

  output:
    stdout yahs_version
  when:
    params.yahs
  """
    touch yahs_version.flag.txt
    echo "YaHS Version:"
    yahs --version
    exit 0;
  """
}

process bwa_Version {
  label 'shortq'
  container = 'staphb/bwa'
  cpus 1

  output:
    stdout bwa_mem_version
  when:
    params.yahs
  """
    touch bwa_version.flag.txt
    echo "bwa Version:"
    bwa 2>&1 | head -n 5
    exit 0;
  """
}


process bcftools_Version {
  label 'shortq'
  container = 'mgibio/bcftools:1.9'
  cpus 1

  output:
    stdout bcftools_version
  when:
    params.polish

  """
    touch bcftools_version.flag.txt
    echo "Bcftools Version:"
    bcftools --version
    exit 0;
  """
}

process hicstuff_Version {
  label 'shortq'
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

process gfastats_Version {
  label 'shortq'
  container = 'dmolik/gfastats'
  cpus 1

  output:
    stdout gfastats_version

  """
    touch gfastats_version.flag.txt
    echo "gfastats Version"
    gfastats --version 2>&1
    exit 0;
  """
}

process deepvariant_Version {
   label 'shortq'
   container = 'google/deepvariant'
   cpus 1

   output:
     stdout deepvariant_version
   when:
     params.polishtype == 'dv'

   """
     touch deepvariant_version.flag.txt
     echo "DeepVariant Version"
     /opt/deepvariant/bin/run_deepvariant --version
     exit 0;
   """
}

process fcs_adaptor_Version {
   label 'shortq'
   container = 'ncbi/fcs-adaptor'
   cpus 1

   output:
     stdout fcs_adaptor_version
   when:
     params.polish == true

   """
     touch fcs_adaptor_version.flag.txt
     echo "fcs-adaptor Version:"
     /app/fcs/bin/av_screen_x --help
     echo "note, some versions of fcs-adaptor do not display version number"
     exit 0;
   """
}

process k_mer_Version {
  label 'shortq'
  cpus 1

  input:
    file version from jellyfish_ver_ch
  output:
    stdout jellyfish_version

  """
    touch K_mer_counting_tool_version.flag.txt
    echo "K-mer counting tool version:"
    cat $version
  """
}

process genomescope_Version {
  label 'shortq'
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

process merfin_Version {
  label 'shortq'
  cpus 1

  output:
    stdout merfin_version
  when:
    params.polishtype == 'merfin'

  """
    touch merfin_version.flag.txt
    echo "merfin  Version:"
    echo "beta?"
    exit 0;
  """

}

process BUSCO_Version {
  label 'shortq'
  cpus 1

  output:
    stdout busco_version
  when:
    params.busco

  """
    touch busco_version.flag.txt
    echo "BUSCO  Version:"
    echo "v5.2.2_cv1"
    exit 0;
  """

}

process shhquis_Version {
  label 'shortq'
  cpus 1

  output:
    stdout shhquis_version
  when:
    params.polish

  """
    touch shhquis_version.flag.txt
    echo "Shhquis.jl Version:"
    echo "0.1.0"
    exit 0;
  """
}

process HiFiAdapterFilt_Version {
  label 'shortq'
  container = 'dmolik/pbadapterfilt'
  cpus 1

  output:
    stdout pbadapterfilt_version

  """
    touch hifiadapterfilt_version.flag.txt
    echo "HiFiAdapterFilt Version:"
    pbadapterfilt.sh -version
    exit 0;
  """
}

process Juicer_Tools_Version {
  label 'shortq'
  container = 'dmolik/juicer-tools'
  cpus 1

  output:
    stdout Juicer_Tools_version

  when:
    params.yahs
  """
    touch juicer_tools_version.flag.txt
    java -Xms512m -Xmx2048m -jar /home/genomics/juicer_tools_1.22.01.jar --version
    exit 0;
  """
}

pbadapterfilt_output
   .collectFile(name:'filtering_information.log.txt', newLine: true, storeDir:"${params.outdir}/00_ordination/log/filtering")

check_fastq_output
   .collectFile(name:'hic_check.log.txt', newLine: true, storeDir:"${params.outdir}/00_ordination/log/filtering")

check_in_file_output
   .collectFile(name:'hifi_check.log.txt', newLine: true, storeDir:"${params.outdir}/00_ordination/log/filtering")

HiFiASM_output
   .collectFile(name:'HiFiASM.log.txt', newLine: true, storeDir:"${params.outdir}/01_hifiasm/log")

gfa2fasta_output
   .collectFile(name:'gfa2fasta.log.txt', newLine: true, storeDir:"${params.outdir}/01_hifiasm/log")

hicstuff_output
   .collectFile(name:'hicstuff.log.txt', newLine: true, storeDir:"${params.outdir}/02_hicstuff/log")

/*hicstuff_polish_output
   .collectFile(name:'hicstuff_polish.log.txt', newLine: true, storeDir:"${params.outdir}/02_hicstuff/log")
*/

ragtag_dot_py_output
   .collectFile(name:'ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/02_hicstuff/log")

faidx_output
   .collectFile(name:'faidx.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

Shhquis_dot_jl_output
   .collectFile(name:'shhquis.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

simple_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

merfin_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

dv_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

busco_gfa_output
   .collectFile(name:'busco.log.txt', newLine: true, storeDir:"${params.outdir}/01_hifiasm/log")

simple_busco_fasta_output
   .collectFile(name:'polished.busco.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

merfin_busco_fasta_output
   .collectFile(name:'polished.busco.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

dv_busco_fasta_output
   .collectFile(name:'polished.busco.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

jellyfish_output
   .collectFile(name:'kmer_counter.log.txt', newLine: true, storeDir:"${params.outdir}/00_ordination/log/genomescope" )

genomescope2_output
   .collectFile(name:'genomescope2.log.txt', newLine: true, storeDir:"${params.outdir}/00_ordination/log/genomescope" )

simple_fcs_adaptor_output
   .collectFile(name:'fcs_adaptor.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

minimap_dot_sh_output
   .collectFile(name:'minimap.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

samtools_mpileup_output
   .collectFile(name:'mpileup.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

bcftools_refmt_output
   .collectFile(name:'bcftools_refmt.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

merfin_output
   .collectFile(name:'merfin.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

minimap_dv_output
   .collectFile(name:'bbmap_dot_sh.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log" )

deep_variant_output
   .collectFile(name:'deepvariant.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

merfin_bcftools_output
   .collectFile(name:'bcftools.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

merfin_fcs_adaptor_output
   .collectFile(name:'fcs_adaptor.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

dv_bcftools_output
   .collectFile(name:'bcftools.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

dv_fcs_adaptor_output
   .collectFile(name:'fcs_adaptor.log.txt', newLine: true, storeDir:"${params.outdir}/03_polish/log")

bwa_for_yahs_output
   .collectFile(name:'bwa.log.txt', newLine: true, storeDir:"${params.outdir}/04_yahs/log")

bwa_for_yahs_simple_output
   .collectFile(name:'bwa.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

bwa_for_yahs_merfin_output
   .collectFile(name:'bwa.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

bwa_for_yahs_dv_output
   .collectFile(name:'bwa.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

bam_sort_for_yahs_output
   .collectFile(name:'bam.sort.log.txt', newLine: true, storeDir:"${params.outdir}/04_yahs/log")

bam_sort_for_simple_yahs_output
   .collectFile(name:'bam.sort.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

bam_sort_for_merfin_yahs_output
   .collectFile(name:'bam.sort.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

bam_sort_for_dv_yahs_output
   .collectFile(name:'bam.sort.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_output
   .collectFile(name:'yahs.log.txt', newLine: true, storeDir:"${params.outdir}/04_yahs/log")

yahs_simple_output
   .collectFile(name:'yahs.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_merfin_output
   .collectFile(name:'yahs.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_dv_output
   .collectFile(name:'yahs.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

juicer_tools_pre_yahs_output
   .collectFile(name:'juicer_tools.log.txt', newLine: true, storeDir:"${params.outdir}/04_yahs/log")

juicer_tools_pre_yahs_simple_output
   .collectFile(name:'juicer_tools.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

juicer_tools_pre_yahs_merfin_output
   .collectFile(name:'juicer_tools.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

juicer_tools_pre_yahs_dv_output
   .collectFile(name:'juicer_tools.polished.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_ragtag_dot_py_hap_output
   .collectFile(name:'yahs.ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_simple_ragtag_dot_py_hap_output
   .collectFile(name:'yahs.simple.ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_merfin_ragtag_dot_py_hap_output
   .collectFile(name:'yahs.merfin.ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_dv_ragtag_dot_py_hap_output
   .collectFile(name:'yahs.dv.ragtag.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log")

yahs_busco_fasta_output
   .collectFile(name:'yahs.busco.log.txt', newLine: true, storeDir:"${params.outdir}/04_yahs/log")

yahs_simple_busco_fasta_output
   .collectFile(name:'yahs.polished.busco.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log" )

yahs_merfin_busco_fasta_output
   .collectFile(name:'yahs.polished.busco.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log" )

yahs_dv_busco_fasta_output
   .collectFile(name:'yahs.polished.busco.log.txt', newLine: true, storeDir:"${params.outdir}/05_yahs_on_polish/log" )

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

yahs_version
   .collectFile(name:'yahs_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

bwa_mem_version
    .collectFile(name:'bwa_mem_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
    .view{ it.text }

bcftools_version
   .collectFile(name:'bcftools_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

hicstuff_version
   .collectFile(name:'hicstuff_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

gfastats_version
   .collectFile(name:'gfastats_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

fcs_adaptor_version
   .collectFile(name:'fcs-adaptor_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

deepvariant_version
   .collectFile(name:'deepvariant_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

jellyfish_version
   .collectFile(name:'jellyfish_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

genomescope_version
   .collectFile(name:'genomescope_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }

merfin_version
   .collectFile(name:'merfin_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
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

Juicer_Tools_version
   .collectFile(name:'juicer_tools_version.txt', newLine: true, storeDir: "${params.outdir}/software_versions")
   .view{ it.text }
