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

  required:
    -f or --forward
       a fastq or fastq.gz file for the pipline, order does not matter

    -r or --reverse
       another fastq or fastq.gz file for the pipeline, order does not matter

    -in or --reads
       path to reads (generally from pacbio), may include a wildcard for multiple files, can be fastq or bam


  suggested:
    -m or --mode
       mode to use, must be one of \"phasing\",\"homozygous\",\"heterozygous\",\"trio\", default: homozygous

    -t or --threads
       number of threads to use, clusters sometimes use this as number of cores, default: 20

    -n or --name
       a name for the assembly

    -y or --yahs
       run yahs as well

    grid computing:
       select one of the following, defaults to local which is highly not-recomended
    -u or --runner
       runner choices are:
          \"sge\": generic sge
          \"slurm\": generic slurm
          \"slurm_usda\": USDA Ceres Cluster
          \"slurm_usda_atlas\": USDA Atlas Cluster, unstable
          \"slurm_usd_mem\": USDA Ceres Large Memory on Ceres
          \"none\": Local running
        you can also overload runner by giving a custom config in the /config directory without the .cfg.
        An exmple: there exists a config in the config directory called torque.cfg the following would
        be used:
          --runner torque
        the /config would not be included, otb.sh only looks in the config directory, and the .cfg is appended

    --polish-type
       turn on polishing one of:
          \"simple\": ragtag patching and shhquis rearangement
          \"merfin\": merfin variant calls ontop of ragtag, bcftools consensus
          \"dv\": deep variant calls ontop of ragtag, bcftools consensus

    busco:
       busco options, requires a lineage option
    --busco
       busco flag turns on busco
    select one of the following (all imply --busco):
       --auto-lineage
          try to use auto lineage finder from busco
       --auto-lineage-prok
          try to use auto lineage finder from busco, but limit to prokaryotes
       --auto-lineage-euk
          try to use auto linage finder from busco, but limit to eukaryotes
       -l or --lineage
          use a specific lineage with busco (recomended)
       -p or --busco-path
          run busco in offline mode, with path to database, or with database name to try and download

    K-mer counting
        options for k-mer counting
    -k or --kmer
        one of:
          \"kmc\": use kmc for k-mer counting
          \"jellyfish\": use jellyfish for k-mer counting

    Shhquis.jl:
       Shhquis options, optional
       --hclust-linkage
          linkage is one of single, average, complete, ward, ward_presquared. Deafults to average.
    "
  exit 0;
}

while [ $# -gt 0 ] ; do
  case $1 in
    -h | --help) help ;;
    -v | --version) version ;;
    -s | --supress) SUPRESS="true";;
    -c | --check) TEST="true";;
    -f | --forward) R1="$2" ;;
    -r | --reverse) R2="$2" ;;
    -in | --reads) READS="$2" ;;
    -m | --mode) MODE="$2";;
    -t | --threads) THREADS="$2";;
    -n | --name) NAME="$2";;
    -y | --yahs) YAHS="true";;
    -u | --runner) RUNNER="$2";;
    -k | --kmer) KMER="$2";;
    --busco) BUSCO="--busco ";;
    --polish-type) POLISHTYPE="$2";;
    --auto-lineage) LINEAGE="auto-lineage";;
    --auto-lineage-prok) LINEAGE="auto-lineage-prok";;
    --auto-lineage-euk) LINEAGE="auto-lineage-euk";;
    -l | --lineage) LINEAGE="$2";;
    -p | --busco-path) BUSCOPATH="$2";;
    --hclust-linkage) SHHQUISHCLST="$2";;
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
    "slurm_usda_atlas") state "$RUNNER being used";;
    "slurm_usda_mem") state "$RUNNER being used";;
    "none") state "$RUNNER being used";;
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
else
  state "   ...not polishing"
fi

state "building run parameters"

#RUN="nextflow run run.nf -with-report ./reports/nextflow-${NAME}.report.html -with-trace ./reports/nextflow-${NAME}.trace.txt -with-timeline ./reports/nextflow-${NAME}.timeline.html -with-dag ./reports/nextflow-${NAME}.dag.png "
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
  warn "mode not set, assuming heterozygous run"
  RUN+="--mode=\"heterozygous\" "
fi

[ -n "$SHHQUISHCLST" ] || SHHQUISHCLST="average"
case $SHHQUISHCLST in
  single) RUN+="--hclust-linkage=\"single\" ";;
  average) RUN+="--hclust-linkage=\"average\" ";;
  complete) RUN+="--hclust-linkage=\"complete\" ";;
  ward) RUN+="--hclust-linkage=\"ward\" ";;
  ward_presquared) RUN+="--hclust-linkage=\"ward_presquared\" ";;
  *) error "shhquis linkage type set to $SHHQUISHCLST, not a vaild hclust linkage type";;
esac

[ -n "$KMER" ] || KMER="kmc"
case $KMER in
  kmc) RUN+="--kmer=\"kmc\" ";;
  jellyfish) RUN+="--kmer=\"jellyfish\" ";;
  *) error "K-mer counting tool $KMER does not valid";;
esac

RUN+="--outfasta=\"${NAME}.genome.out\" "
[ -n "$THREADS" ] && RUN+="--threads=\"$THREADS\" " || warn "threads not set, setting to 20 maximum threads"
[ -z "$THREADS" ] && RUN+="--threads=\"20\" "
[ -f "$R1" ] && RUN+="--readf=\"$R1\" " || error "read pair file one not found, exiting"
[ -f "$R2" ] && RUN+="--readr=\"$R2\" " || error "read pair file two not found, exiting"
[ -n "$YAHS" ] && RUN+="--yahs "
[ -n "$BUSCO" ] && RUN+="$BUSCO " || state "not running busco"
[ -n "$POLISHTYPE" ] && RUN+="--polish --polishtype=\"$POLISHTYPE\" " || warn "not polishing, it is recomended that you polish"
[ -n "$BUSCOSTRING" ] && RUN+="$BUSCOSTRING"
[ -n "$READS" ] && RUN+="--readin=\"$READS\" " || error "reads file(s) not given, exiting"
[ -z "$NAME" ] && NAME="$(date +%s)" && state "name not given, setting name to: $NAME"
RUN+="--assembly=\"$NAME\" "

pizzaz "$RUN"
[ -z "$SUPRESS" ] && stop_check "check that the command is expected, continue"

state "Prefetching singularity containers"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] && "Nextflow Singularity cache directory set: $NXF_SINGULARITY_CACHEDIR, will use for singularity images" || warn "NXF_SINGULARITY_CACHEDIR not set, using ./work/singularity instead"

prefetch_container="./scr/prefetch_containers.sh"
[ -n "$YAHS" ] && prefetch_container+=" -y"
[ -n "$BUSCO" ] && prefetch_container+=" -b"
[ -n "$POLISHTYPE" ] && prefetch_container+=" -p $POLISHTYPE"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] || ( mkdir -p "./work/singularity"; prefetch_container+=" -l ./work/singularity" )
eval $prefetch_container

if [ -n "$TEST" ]; then
  check_container="./scr/check_containers.sh"
  [ -n "$BUSCO" ] && check_container+=" -b"
  [ -n "$POLISHTYPE" ] && check_container+=" -p $POLISHTYPE"
  [ -n "$NXF_SINGULARITY_CACHEDIR" ] || check_container+=" -l ./work/singularity"
  eval $check_container
fi

state "checking for running busco"
if [ -n "$BUSCO" -o -n "$LINEAGE" -o -n "$BUSCOPATH" ]; then
state "running busco, checking busco things"
  [ -z "$LINEAGE" -a -z "$BUSCOPATH" ] && error "trying to setup busco, but no lineage/path/file given"
  if [ -n "$LINEAGE" ]; then
    state "   ...busco lineage described: ${LINEAGE}"
    BUSCOSTRING="--busco --linreage=\"${LINEAGE} \""
  elif [ -n "$BUSCOPATH" ]; then
    state "   ...attempting offline busco run"
    if [ -f "$BUSCOPATH" ]; then
      state "you have handed busco a file"
      BUSCODWNLOC=$( dirname $BUSCOPATH )
      BUSCODBFILE=$( basename $BUSCOPATH )
      BUSCOSTRING="--busco --linreage=\"${BUSCODBFILE}\" --buscooffline --buscodb=\"${BUSCODWNLOC}\" "
    else
      warn "trying to download busco dataset, this is not recomended"
      mkdir -p work/busco
      if [ -f "work/singularity/ezlabgva-busco-v5.2.2_cv1.img" ]; then
        state "busco found at work/singularity/ezlabgva-busco-*"
        singularity exec work/singularity/ezlabgva-busco-v5.2.2_cv1.img busco --download_path work/busco --download "${BUSCOPATH}" || error "unable to download busco dataset when asked too, exiting"
      elif [ -f "${SINGULARITY_LOCALCADHEDIR}/ezlabgva-busco-v5.2.2_cv1.img" ]; then
        state "busco found at ${SINGULARITY_LOCALCADHEDIR}/ezlabgva-busco-*"
        singularity exec ${SINGULARITY_LOCALCADHEDIR}/ezlabgva-busco-v5.2.2_cv1.img busco --download_path work/busco --download "${BUSCOPATH}" || error "unable to download busco dataset when asked too, exiting"
      elif [ -f "${SINGULARITY_CADHEDIR}/ezlabgva-busco-v5.2.2_cv1.img" ]; then
        state "busco found at ${SINGULARITY_CADHEDIR}/ezlabgva-busco-*"
        singularity exec ${SINGULARITY_CADHEDIR}/ezlabgva-busco-v5.2.2_cv1.img busco --download_path work/busco --download "${BUSCOPATH}" || error "unable to download busco dataset when asked too, exiting"
      elif [ $( command -v busco >/dev/null ) ]; then
        state "busco command found"
        warn "using non-otb version of busco, this may result in problems"
        busco --download_path work/busco --download "${BUSCOPATH}" || error "unable to download busco dataset when asked too, exiting"
      else
        error "no busco command found, I can not dowload any datasets, exiting"
      fi
      BUSCOSTRING="--busco --linreage=\"${BUSCOPATH}\" --buscoffline --buscodb=\"work/busco\" "
    fi
  fi
else
  state "   ...not running busco"
fi

[ -z "$SUPRESS" ] && RUN+="-bg"

[ -z "$SUPRESS" ] && stop_check "proceed with run"
state "making /reports dir"
mkdir -p reports
pizzaz "running only the best"
echo $RUN > "./reports/${NAME}.nextflow.command.txt"
echo $RUN > "./reports/nextflow-${NAME}.stdout.txt"
eval $RUN &> "./reports/nextflow-${NAME}.stdout.txt"
