#!/usr/bin/env nextflow

params.assembly = "an_assembly"
params.readbam = "$baseDir/data/*.bam"
params.readf = "$baseDir/data/*.R1.fastq.gz"
params.readr = "$baseDir/data/*.R2.fastq.gz"
params.outfasta = "genome.out.fasta"
params.outdir = 'results'
params.mode = 'heterozygous'
params.ploidy = '2'
params.threads = '20'
params.linreage = 'insecta_odb10'
params.busco = false
params.polish = false
params.polishtype = 'simple'

bam_ch = Channel.fromPath(params.readbam)
right_fastq_check = Channel.fromPath(params.readr)
left_fastq_check = Channel.fromPath(params.readf)

bam_ch.into{ bam_check_ch; bam_Hifi_ch; bam_dv_ch}

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
   [[ !{right_fastq}  =~ ".gz" ]] && ln -s !{right_fastq} right.fastq.gz || (zcat !{right_fastq} > right.fastq.gz)
   [[ !{left_fastq}  =~ ".gz" ]] && ln -s !{left_fastq} left.fastq.gz || (zcat !{left_fastq} > left.fastq.gz)

   state "successful completion"

   exit 0;
  '''
}

process HiFiAdapterFilt {
  container = 'dmolik/pbadapterfilt'
  cpus = params.threads

  input:
    file bam from bam_Hifi_ch.flatten()
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
    file fasta from filt_fasta_ch.collect()
    file left from left_fastq_HiFiASM
    file right from right_fastq_HiFiASM

  output:
    file '*.gfa' into gfa_ch
    file '*.ec.fa' into fasta_ec_ch, meryl_ec_ch, bbmap_ec_ch
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
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'pvstodghill/any2fasta'
  cpus 1

  input:
    file gfa from gfa_ch.flatten()
  output:
    file '*.p_ctg.gfa.fasta' optional true into gfa2fasta_fasta_res_ch
    file '*.bp.p_ctg.gfa.fasta' optional true into fasta_unoriented_ch, fasta_genome_ch, fasta_busco_ch
    file '*hap[12].p_ctg.gfa.fasta' optional true into fasta_hap_ch, simple_fasta_hap_polish_ch, merfin_fasta_hap_polish_ch, dv_fasta_hap_polish_ch
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

  input:
    file genome from fasta_fai_genome_ch
  output:
    file '*.fai' into fai_ch
    stdout faidx_output
  when:
    params.polish
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
    params.polish
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
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'dmolik/shhquis'
  cpus 1

  input:
    file abs from abs_ch
    file contig from contigs_ch
    file genome from fasta_sshquis_genome_ch
    file fai from fai_ch
  output:
    file "${params.outfasta}" into shhquis_fasta_res_ch, polish_haps_genome_ch, shhquis_simple_ch, shhquis_merfin_ch, shhquis_dv_ch, shhquis_bcftools_dv_ch, shhquis_bcftools_merfin_ch, shhquis_bbmap_ch, shhquis_mpileup_ch
    file "${params.outfasta}"
    stdout Shhquis_dot_jl_output
  when:
    params.polish
  """
    touch shhquis.flag.txt
    shh.jl --reorient ${params.outfasta} --genome ${genome} --fai ${fai} --bg2 ${abs} --contig ${contig} --hclust-linkage "average"
    echo "finished reorientation"
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
    exit 0;
  """
}

process simple_polish {
  cpus = 1

  input:
    file genome from shhquis_simple_ch
  output:
    file "${params.assembly}.polished.genome.fasta" into simple_polished_genome_ch, simple_polished_genome_busco_ch
  when:
    params.polishtype == "simple"
  """
    touch simple_polish.flag.txt
    ln -s ${genome} ${params.assembly}.polished.genome.fasta
    echo "finished softlink"
    exit 0;
  """
}

process bbmap_dot_sh {
  container = 'bryce911/bbtools'
  cpus = params.threads

  input:
    file ec_reads from bbmap_ec_ch
    file genome from shhquis_bbmap_ch  
  output:
    file "mapped.sam" into sam_for_merfin_ch
    stdout bbmap_dot_sh_output
  when:
    params.polishtype == "merfin"
  """
    touch bbmap.sh.flag.sh
    bbmap.sh t=${task.cpus} in=${ec_reads} out=mapped.sam ref=${genome}
    echo "finished bbmap.sh"
    exit 0;
  """
}

process samtools_mpileup {
  container = 'mgibio/samtools:1.9'
  cpus = 1

  input:
    file sam_file from sam_for_merfin_ch
    file genome from shhquis_mpileup_ch
  output:
    file 'out.mpileup' into bcf_for_merfin_ch  
    stdout samtools_mpileup_output
  when:
    params.polishtype == "merfin"
  """
    touch samtools.mpileup.flag.txt
    samtools view -S -b ${sam_file} > aln.bam
    samtools mpileup -E -uf reference.fa file.bam > out.mpileup
    echo "finished mpileup"
    exit 0;
  """
}

process bcftools_refmt {
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
    echo "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">" > merfin_header.vcf
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> merfin_header.vcf
    echo "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" >> merfin_header.vcf
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
    exit 0;
  """
}

process merfin {
  publishDir "${params.outdir}/merfin", mode: 'rellink'
  container = 'dmolik/merfin'
  cpus = params.threads

  input:
    file genome from shhquis_merfin_ch
    file kcov_file from kcov_ch
    file lookup_table from lookup_table_ch
    file ec_reads from meryl_ec_ch
    file vcf_file from vcf_for_merfin_ch 
  output:
    file 'merfin.polish.vcf' into merfin_vcf_ch
    stdout merfin_output
  when:
    params.polishtype == "merfin"

  script:
    kcov="""
      cat ${kcov_file}
    """.stripIndent()
    """
      echo "warning merfin is experimental"
      touch merfin.flag.txt
      meryl count k=21 ${ec_reads} output reads.meryl
      meryl greater-than 1 reads.meryl output reads.gt1.meryl
      merfin -polish -threads ${task.cpus} -sequence ${genome} -peak ${kcov} -prob ${lookup_table} -readmers reads.gt1.meryl -vcf ${vcf_file} -output merfin 
      echo "finished merfin"
      exit 0;
    """
}

process samtools_merge_for_deep_variant {
  container = 'mgibio/samtools:1.9'  
  cpus = params.threads

  input:
    file bam_reads from bam_dv_ch.collect()
  output:
    file 'merged.bam' into bam_dv_merged_ch
  when:
    params.polishtype == "dv"
  """
    touch samtools.merge.flag.txt
    samtools merge --threads ${task.cpus} -X merged.bam ${bam_reads}
    echo "finished merging"
    exit 0;
  """
}

process deep_variant {
  publishDir "${params.outdir}/deepvariant", mode: 'rellink'
  container = 'google/deepvariant'
  cpus = params.threads 

  input:
    file genome from shhquis_dv_ch 
    file bam_read from bam_dv_merged_ch 
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
    exit 0;
  """
}

process dv_bcftools {
  publishDir "${params.outdir}/genome", mode: 'rellink' 
  container = 'mgibio/bcftools:1.9'
  cpus = params.threads

  input:
    file genome from shhquis_bcftools_dv_ch
    file vcf from dv_vcf_ch
  output:
    file "${params.assembly}.vcf_polished_assembly.fasta" into dv_vcf_polished_genome_ch, dv_vcf_res_ch, dv_vcf_polished_busco_genome_ch
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
    exit 0;
  """
}

process merfin_bcftools {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'mgibio/bcftools:1.9'
  cpus = params.threads

  input:
    file genome from shhquis_bcftools_merfin_ch
    file vcf from merfin_vcf_ch
  output:
    file "${params.assembly}.vcf_polished_assembly.fasta" into merfin_vcf_polished_genome_ch, merfin_vcf_res_ch, merfin_vcf_polished_busco_genome_ch
    file "${params.assembly}.vcf_polished_assembly.fasta"
    stdout merfin_bcftools_output
  when:
    params.polishtype == "merfin"
  """
    touch merfin.bcftools.flag.txt
    bcftools view --threads ${task.cpus} -Oz ${vcf} > ${vcf}.gz
    bcftools index --threads ${task.cpus} ${vcf}.gz
    bcftools consensus ${vcf}.gz -f ${genome} -H 1 > ${params.assembly}.vcf_polished_assembly.fasta
    echo "finished bcftools from merfin"
    exit 0;
  """
}

process ragtag_dot_py_hap_simple_polish {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from simple_fasta_hap_polish_ch
    file fasta_genome from simple_polished_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/patched*"
    file "${params.assembly}_ragtag_scaffold/patched" into simple_hap_patch_res_ch
    stdout simple_ragtag_dot_py_hap_output
  when:
    params.polishtype == "simple"
  """
    touch ragtag.hap.flag.txt
    ragtag.py scaffold --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    exit 0;
  """
}

process ragtag_dot_py_hap_merfin_polish {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from merfin_fasta_hap_polish_ch
    file fasta_genome from merfin_vcf_polished_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/patched*"
    file "${params.assembly}_ragtag_scaffold/patched" into merfin_hap_patch_res_ch
    stdout merfin_ragtag_dot_py_hap_output
  when:
    params.polishtype == "merfin"
  """
    ragtag.py scaffold --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
    exit 0;
  """
}

process ragtag_dot_py_hap_deep_variant_polish {
  publishDir "${params.outdir}/genome", mode: 'rellink'
  container = 'dmolik/ragtag'
  cpus = params.threads

  input:
    file fasta_hap from dv_fasta_hap_polish_ch
    file fasta_genome from dv_vcf_polished_genome_ch
  output:
    file "${params.assembly}_ragtag_scaffold/patched*"
    file "${params.assembly}_ragtag_scaffold/patched" into dv_hap_patch_res_ch
    stdout dv_ragtag_dot_py_hap_output
  when:
    params.polishtype == "dv"
  """
    touch ragtag.hap.flag.txt
    ragtag.py scaffold --aligner unimap -t ${task.cpus} -o ./${params.assembly}_ragtag_scaffold ${fasta_genome} ${fasta_hap}
    mv ${params.assembly}_ragtag_scaffold/ragtag.scaffold.fasta ${params.assembly}_ragtag_scaffold/polished.${fasta_hap}
    echo "finished patching"
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

process merfin_busco_fasta {
  publishDir "${params.outdir}/busco_polish", mode: 'rellink'
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

process dv_busco_fasta {
  publishDir "${params.outdir}/busco_polish", mode: 'rellink'
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

process merfin_vcf_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from merfin_vcf_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch vcf_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process dv_vcf_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from dv_vcf_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch vcf_stats.flag.txt
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


process merfin_hap_patch_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from merfin_hap_patch_res_ch.flatten()
  output:
    file '*.stats'
  """
    touch hap_patch_stats.flag.txt
    stats.sh -Xmx4g ${fasta} > ${fasta}.stats
    echo "finished stats"
    exit 0;
  """
}

process dv_hap_patch_stats_dot_sh {
  publishDir "${params.outdir}/genome", mode: 'copy'
  container = 'bryce911/bbtools'
  cpus 1

  input:
    file fasta from dv_hap_patch_res_ch.flatten()
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

process deepvariant_Version {
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

process merfin_Version {
  cpus 1

  output:
    stdout merfin_version
  when:
    params.polishtype == 'merfin'

  """
    touch merfin_version.flag.txt
    echo "merfin  - - - - - - beta?"
    exit 0;
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
  when:
    params.polish

  """
    touch shhquis_version.flag.txt
    echo "Shhquis.jl - - - - - 0.1.0"
    exit 0;
  """
}

process HiFiAdapterFilt_Version {
  cpus 1

  output:
    stdout pbadapterfilt_version  

  """
    touch hifiadapterfilt_version.flag.txt
    echo "HiFiAdapterFilt  - - v1.0.0"
    exit 0;
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

simple_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

merfin_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

dv_ragtag_dot_py_hap_output
   .collectFile(name:'ragtag_hap.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

busco_gfa_output
   .collectFile(name:'busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_no_polish")

simple_busco_fasta_output
   .collectFile(name:'simple.busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_polish" )

merfin_busco_fasta_output
   .collectFile(name:'merfin.busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_polish" )

dv_busco_fasta_output
   .collectFile(name:'dv.busco.log.txt', newLine: true, storeDir:"${params.outdir}/busco_polish" )

jellyfish_output
   .collectFile(name:'jellyfish.log.txt', newLine: true, storeDir:"${params.outdir}/genomescope" )

genomescope2_output
   .collectFile(name:'genomescope2.log.txt', newLine: true, storeDir:"${params.outdir}/genomescope" )

bbmap_dot_sh_output
   .collectFile(name:'bbmap_dot_sh.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log" )

samtools_mpileup_output
   .collectFile(name:'mpileup.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log" )

bcftools_refmt_output
   .collectFile(name:'bcftools_refmt.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log" )

merfin_output
   .collectFile(name:'merfin.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

deep_variant_output
   .collectFile(name:'deepvariant.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

merfin_bcftools_output
   .collectFile(name:'bcftools.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

dv_bcftools_output
   .collectFile(name:'bcftools.log.txt', newLine: true, storeDir:"${params.outdir}/genome/log")

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
