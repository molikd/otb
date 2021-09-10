#!/bin/bash
source scr/io.sh

help(){
  describe "otb: Only The Best (Genome Assemblies)
  help"
  exit 0;
}

while [ $# -gt 0 ] ; do
  case $1 in
    -h | --help) help ;;
    -v | --version) version ;;
    -f | --forward) R1="$2" ;;
    -r | --reverse) R2="$2" ;;
    -b | --bam) BAM="$2" ;;
    -m | --mode) MODE="$2";;
    -t | --threads) THREADS="$2";;
    -n | --name) NAME="$2";;
    --sge) RUNNER="sge";;
    --slurm) RUNNER="slurm";;
    --slurm-usda) RUNNER="slurm_usda";;
    --busco) BUSCO="--busco=\"true\"";;
    --polish) POLISH="--polish=\"true\"";;
    --auto-lineage) LINEAGE="auto-lineage";;
    --auto-lineage-prok) LINEAGE="auto-lineage-prok";;
    --auto-lineage-euk) LINEAGE="auto-lineage-euk";;
    -l | --lineage) LINEAGE="$2";;
  esac
  shift
done

display_header

state "checking for nextflow"
command -v nextflow >/dev/null 2>&1 || error "nextflow could not be found, aborting"
state "using $(which nextflow) for nextflow"

state "checking for singularity"
command -v singularity >/dev/null 2>&1 || error "singularity could not be found, aborting" 
state "using $(which singularity) for singularity"

state "output user environment"
bash scr/check_env.sh 

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
[ -n "$POLISH" ] && RUN+="$POLISH " || warn "not polishing, this is not recomended"
[ -n "$BUSCO" ] && [ -z "$LINEAGE" ] && error "you want to run BUSCO, but busco lineage not set, exiting"
[ -n "$LINEAGE" ] && RUN+="--linreage=\"$LINEAGE\" "
[ -n "$BAM" ] && RUN+="--readbam=\"$BAM\" " || error "bam file(s) not given, exiting"
[ -z "$NAME" ] && NAME="$(date +%s)" && state "name not given, setting name to: $NAME"
RUN+="--assembly=\"$NAME\" "
RUN+="-bg"

pizzaz "$RUN"
stop_check "check that the command is expected, continue"

state "Prefetching singularity containers"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] && "Nextflow Singularity cache directory set: $NXF_SINGULARITY_CACHEDIR, will use for singularity images" || warn "NXF_SINGULARITY_CACHEDIR not set, using ./work/singularity instead"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] && ./scr/prefetch_containers.sh || ( mkdir -p "./work/singularity"; ./scr/prefetch_containers.sh "./work/singularity" )

#TODO, check that all the containers work

stop_check "proceed with run"
eval $RUN
