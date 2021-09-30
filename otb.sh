#!/bin/bash
source scr/io.sh

help(){
  display_header
  describe "otb: Only The Best (Genome Assemblies)
  [NOTE] otb must be run from the directory it is in, only ./otb.sh will work
  [NOTE] run.nf must be in the same directory as otb.sh
  [NOTE] ./config holds grid config files
  [NOTE] ./scr holds scripts on which otb relies 

  utilities 
    -h or --help
       display this message and exit

    -v or --version
       display otb version and exit

    -s or --supress
       supress stop and checks, and submit without sending to background

    -c or --check
       perform checks to insure smoother operation
 
    -b or --background
       send nextflow run to the background, ensures that nextflow will continue to run if disconnected

  required:
    -f or --forward
       a fastq or fastq.gz file for the pipline, order does not matter

    -r or --reverse
       another fastq or fastq.gz file for the pipeline, order does not matter

    -b or --bam
       path to bam file, may include a wildcard for multiple ban files


  suggested:
    -m or --mode
       mode to use, must be one of \"phasing\",\"homozygous\",\"heterozygous\",\"trio\", default: homozygous

    -t or --threads
       number of threads to use, clusters sometimes use this as number of cores, default: 20

    -n or --name
       a name for the assembly 

    grid computing:
       select one of the following, defaults to local which is highly not-recomended
    --sge
    --slurm
    --slurm-usda

    --polish-type
       turn on polishing one of:
          \"simple\": ragtag patching and shhquis rearangement
          \"merfin\": merfin variant calls ontop of ragtag, bcftools consensus
          \"dv\": deep variant calls ontop of ragtag, bcftools consensus

    busco:
       busco options, requires a lineage option
    --busco
       busco flag turns on busco
    select one of the following:
       --auto-lineage
          try to use auto lineage finder from busco
       --auto-lineage-prok
          try to use auto lineage finder from busco, but limit to prokaryotes
       --auto-lineage-euk
          try to use auto linage finder from busco, but limit to eukaryotes
       -l or --lineage
          use a specific lineage with busco (recomended)"
  exit 0;
}

while [ $# -gt 0 ] ; do
  case $1 in
    -h | --help) help ;;
    -v | --version) version ;;
    -s | --supress) SUPRESS="true";;
    -c | --check) TEST="true";;
    -b | --background) BG="-bg";;
    -f | --forward) R1="$2" ;;
    -r | --reverse) R2="$2" ;;
    -b | --bam) BAM="$2" ;;
    -m | --mode) MODE="$2";;
    -t | --threads) THREADS="$2";;
    -n | --name) NAME="$2";;
    --sge) RUNNER="sge";;
    --slurm) RUNNER="slurm";;
    --slurm-usda) RUNNER="slurm_usda";;
    --busco) BUSCO="--busco ";;
    --polish-type) POLISHTYPE="$2";;
    --auto-lineage) LINEAGE="auto-lineage";;
    --auto-lineage-prok) LINEAGE="auto-lineage-prok";;
    --auto-lineage-euk) LINEAGE="auto-lineage-euk";;
    -l | --lineage) LINEAGE="$2";;
  esac
  shift
done

display_header

state "checking for run.nf"
[ -f "./run.nf" ] || error "run.nf not found aborting"
state "using ./work/ for work directory"
state "checking for nextflow"
command -v nextflow >/dev/null 2>&1 || error "nextflow could not be found, aborting"
state "using $(which nextflow) for nextflow"
state "checking for singularity"
command -v singularity >/dev/null 2>&1 || error "singularity could not be found, aborting" 
state "using $(which singularity) for singularity"
state "output user environment"
bash scr/check_env.sh 

state "checking runner"
if [ -n "$RUNNER" ]; then
  case $RUNNER in
    "sge") state "$RUNNER being used";;
    "slurm") state "$RUNNER being used";;
    "slurm_usda") state "$RUNNER being used";;
    *) 
      [ -f "config/${RUNNER}" ] && state "using custom config: $RUNNER" || error "runner type ${RUNNER} not found";;
  esac
fi

state "checking polishing"
if [ -n "$POLISHTYPE" ]; then
  case $POLISHTYPE in
    "merfin") state "merfin polishing";;
    "simple") state "simple/ragtag polishing";;
    "dv") state "deep variant polishing";;
    *) error "polishing type ${POLISHTYPE} not found, exiting";;
  esac
fi

RUN="nextflow run run.nf "
[ -n "$RUNNER" ] && RUN+="-c config/${RUNNER}.cfg " || warn "no grid computing environment set, using local. this is not recomended."
if [ -n "$MODE" ]; then
  case $MODE in
    phasing) RUN+="--mode=\"phasing\" ";;
    homozygous) RUN+="--mode=\"homozygous\" ";;
    heterozygous) RUN+="--mode=\"heterozygous\" ";;
    trio) RUN+="--mode=\"trio\" ";;
    *) error "mode set to $MODE, not an actual mode";;
  esac
else
  warn "mode not set, assuming homozygous run"
  RUN+="--mode=\"homozygous\" "
fi
[ -n "$THREADS" ] && RUN+="--threads=\"$THREADS\" " || warn "threads not set, setting to 20 maximum threads" 
[ -z "$THREADS" ] && RUN+="--threads=\"20\" "
[ -f "$R1" ] && RUN+="--readf=\"$R1\" " || error "read pair file one not found, exiting"
[ -f "$R2" ] && RUN+="--readr=\"$R2\" " || error "read pair file two not found, exiting"
[ -n "$BUSCO" ] && RUN+="$BUSCO " || state "not running busco"
[ -n "$POLISHTYPE" ] && RUN+="--polish --polishtype=\"$POLISHTYPE\" " || warn "not polishing, it is recomended that you polish"
[ -n "$BUSCO" ] && [ -z "$LINEAGE" ] && error "you want to run BUSCO, but busco lineage not set, exiting"
[ -n "$LINEAGE" ] && RUN+="--linreage=\"$LINEAGE\" "
[ -n "$BAM" ] && RUN+="--readbam=\"$BAM\" " || error "bam file(s) not given, exiting"
[ -z "$NAME" ] && NAME="$(date +%s)" && state "name not given, setting name to: $NAME"
RUN+="--assembly=\"$NAME\" "
[ -n "$BG" ] && RUN+="-bg"

pizzaz "$RUN"
[ -z "$SUPRESS" ] && stop_check "check that the command is expected, continue"

state "Prefetching singularity containers"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] && "Nextflow Singularity cache directory set: $NXF_SINGULARITY_CACHEDIR, will use for singularity images" || warn "NXF_SINGULARITY_CACHEDIR not set, using ./work/singularity instead"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] && ./scr/prefetch_containers.sh || ( mkdir -p "./work/singularity"; ./scr/prefetch_containers.sh "./work/singularity" )

if [ -n "$TEST" ]; then 
  state "Testing that containers work"
  if [ -n "$NXF_SINGULARITY_CACHEDIR" ]; then
    container_location="$NXF_SINGULARITY_CACHEDIR"
  else
    container_location="work/singularity"
  fi
  singularity exec "$container_location/bryce911-bbtools.img" echo "   ...hello from bbtools container" || error "bbtools container broken, exiting"
  singularity exec "$container_location/dmolik-genomescope2.img" echo "   ...hello from genomescope2 container" || error "genomescope2 container broken, exiting"
  singularity exec "$container_location/dmolik-hifiasm.img"  echo "   ...hello from hifiasm container" || error "hifiasm container broken, exiting"
  singularity exec "$container_location/dmolik-jellyfish.img" echo "   ...hello from jellyfish container" || error "jellyfish container broken, exiting"
  singularity exec "$container_location/dmolik-pbadapterfilt.img" echo "   ...hello from pbadapterfilt container" || error "pbadapterfilt container broken, exiting"
  singularity exec "$container_location/mgibio-samtools-1.9.img" echo "   ...hello from samtools container" || error "samtools container broken, exiting"
  singularity exec "$container_location/koszullab-hicstuff.img" echo "   ...hello from hicstuff container" || error "hicstuff container broken, exiting"
  singularity exec "$container_location/pvstodghill-any2fasta.img" echo "   ...hello from any2fasta container" || error "any2fasta container borken, exiting"
  if [ -n "$BUSCO" ]; then 
    state "busco container will be required for this run, testing busco"
    singularity exec "$container_location/ezlabgva-busco-v5.2.2_cv1.img" echo "   ...hello from busco container" || error "busco container broken, exiting"
  fi 
  case $POLISHTYPE in
    "simple")
      state "this will be a simple polish, checking ragtag and shhquis"
      singularity exec "$container_location/dmolik-ragtag.img" echo "   ...hello from ragtag container" || error "ragtag container broken, exiting"
      singularity exec "$container_location/dmolik-shhquis.img" echo "    ...hello from shhquis container" || error "shhquis container broken, exiting"
    ;;
    "merfin") 
      state "this will be a merfin polish, checking ragtag, shhquis, bcftools, and merfin"
      singularity exec "$container_location/dmolik-ragtag.img" echo "   ...hello from ragtag container" || error "ragtag container broken, exiting"
      singularity exec "$container_location/dmolik-shhquis.img" echo "    ...hello from shhquis container" || error "shhquis container broken, exiting"
      singularity exec "$container_location/mgibio-bcftools-1.9.img" echo "   ...hello from bcftools container" || error "bcftools container broken, exiting"
      singularity exec "$container_location/dmolik-merfin.img" echo "    ...hello from merfin container" || error "merfin container broken, exiting"
    ;;
    "dv")
      state "this will be a deep variant polish, checking ragtag, shhquis, bcftools, and merfin"
      singularity exec "$container_location/dmolik-ragtag.img" echo "   ...hello from ragtag container" || error "ragtag container broken, exiting"
      singularity exec "$container_location/dmolik-shhquis.img" echo "    ...hello from shhquis container" || error "shhquis container broken, exiting"
      singularity exec "$container_location/mgibio-bcftools-1.9.img" echo "   ...hello from bcftools container" || error "bcftools container broken, exiting"
      singularity exec "$container_location/google-deepvariant.img" echo "   ...hello from deepvariant container" || error "deepvariant container broken, exiting"
    ;;
  state "all required containers checked with an intial pass"
fi 

[ -z "$SUPRESS" ] && stop_check "proceed with run"
eval $RUN
