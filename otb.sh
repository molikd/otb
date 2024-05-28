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
       supress \"stop and gather user input\", and submit without sending to background

    -c or --check
       perform checks to insure smoother operation

    --lite
       don't perform anything with HiC

  required:
    -in or --reads
       path to HiFi reads (generally from pacbio), may include a wildcard for multiple files, can be fastq or bam

  required if otb-lite isn't being used:
    -f or --forward
       a fastq or fastq.gz file for the pipline, in HiC runs this is one of the sequencing files, in a trio run, this is either the maternal or paternal sequences, order does not matter

    -r or --reverse
       another fastq or fastq.gz file for the pipeline, in HiC runs this is one of the sequencing files, in a trio run, this is either the maternal or paternal sequences, order does not matter

  suggested:
    --matf and --matr
       optional maternal sequences forward and reverse

    --patf and --patr
       optional paternal sequences forward and reverse

    -m or --mode
       mode to use, must be one of \"phasing\",\"default\",\"trio\",\"primary\". default: \"default\"

    -t or --threads
       number of threads to use, clusters sometimes use this as number of cores, default: 20

    -n or --name
       a name for the assembly

    -y or --yahs
       run yahs as well

    --purge-dups
       [0-3], the amount of dups purging requested in HiFiASM, by default no purging is done. default: 0

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

    --patch
       perform error corrected read patching as part of polishing

    --hapscaffold
       turn on ragtag scaffolding for the haplotypes to algin and add N's to the haplotypes from the primary/parental genome

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
       --evalue
          change default busco evalue, default is 0.001

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
    --lite) LITE="true";;
    --matf) MATR1="$2";;
    --marr) MATR2="$2";;
    --patf) PATR1="$2";;
    --patr) PATR2="$2";;
    --purge-dups) PURGE_DUPS="$2";;
    --busco) BUSCO="--busco ";;
    --hapscaffold) HAPSCAFFOLD="true";;
    --patch) PATCH="true";;
    --evalue) BUSCOEVALUE="$2";;
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

if [ -n "$PATCH" ]; then
  [ -n $POLISHTYPE ] || error "can't patch without polishing first, select a polish type";
  RUN+="--patch ";
fi

state "building run parameters"

#RUN="nextflow run run.nf -with-report ./reports/nextflow-${NAME}.report.html -with-trace ./reports/nextflow-${NAME}.trace.txt -with-timeline ./reports/nextflow-${NAME}.timeline.html -with-dag ./reports/nextflow-${NAME}.dag.png "
RUN="nextflow run run.nf "
[ -n "$RUNNER" ] && RUN+="-c config/${RUNNER}.cfg " || warn "no grid computing environment set, using local. this is not recomended."
if [ -n "$MODE" ]; then
  case $MODE in
    phasing) RUN+="--mode=\"phasing\" ";;
    default) RUN+="--mode=\"default\" ";;
    trio) RUN+="--mode=\"trio\" ";;
    primary) RUN+="--mode=\"primary\" ";;
    *) error "mode set to $MODE, not an actual mode";;
  esac
else
  warn "mode not set, assuming default run"
  RUN+="--mode=\"default\" "
fi

state "checking required parameters for selected mode: $MODE"

case $MODE in
  phasing)
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then error "phasing mode can not be run without HiC data, not found"; fi
    if [ -f "$R1" ] && [ -f "$R2" ]; then warn "not using HiC data in assembly step, only as HiC polishing"; fi
    if [ -f "$MATR" ] || [ -f "$MATF" ] || [ -f "$PATR" ] || [ -f "$PATF" ]; then warn "some trio data given, but HiFiASM mode $MODE, does not utilize trio data, data will be ignored"; fi
  ;;
  default)
     if [ -z "$LITE" ]; then if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then error "trying to run HiC polishing without HiC data in non-lite mode, some files not found"; fi; fi
     if [ -f "$R1" ] && [ -f "$R2" ]; then warn "not using HiC data in assembly step, only as HiC polishing"; fi
     if [ -f "$MATR" ] || [ -f "$MATF" ] || [ -f "$PATR" ] || [ -f "$PATF" ]; then warn "some trio data given, but HiFiASM mode $MODE, does not utilize trio data, data will be ignored"; fi
  ;;
  trio)
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then error "trying to run trio mode without HiC data"; fi
    if [ ! -f "$MATR" ] || [ ! -f "$MATF" ] || [ ! -f "$PATR" ] || [ ! -f "$PATF" ]; then error "trying to run tro mode without Mat/Pat data"; fi
  ;;
  primary)
    if [ -z "$LITE" ]; then if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then error "trying to run HiC polishing without HiC data in non-lite mode, some files not found"; fi; fi
    if [ -f "$R1" ] && [ -f "$R2" ]; then warn "not using HiC data in assembly step, only as HiC polishing"; fi
    if [ -f "$MATR" ] || [ -f "$MATF" ] || [ -f "$PATR" ] || [ -f "$PATF" ]; then warn "some trio data given, but HiFiASM mode $MODE, does not utilize trio data, data will be ignored"; fi
  ;;
  *)
  error "Mode not found"
  ;;
esac

[[ -f "$R1" ]] && R1=$( readlink -f $R1 )
[[ -f "$R2" ]] && R2=$( readlink -f $R2 )
[[ -f "$MATR" ]] && MATR=$( readlink -f $MATR )
[[ -f "$PATR" ]] && PATR=$( readlink -f $PATR )
[[ -f "$MATF" ]] && MATF=$( readlink -f $MATF )
[[ -f "$PATF" ]] && PATF=$( readlink -f $PATF )

[ -n "$SHHQUISHCLST" ] || SHHQUISHCLST="average"
case $SHHQUISHCLST in
  single) RUN+="--hclustlinkage=\"single\" ";;
  average) RUN+="--hclustlinkage=\"average\" ";;
  complete) RUN+="--hclustlinkage=\"complete\" ";;
  ward) RUN+="--hclustlinkage=\"ward\" ";;
  ward_presquared) RUN+="--hclustlinkage=\"ward_presquared\" ";;
  *) error "shhquis linkage type set to $SHHQUISHCLST, not a vaild hclust linkage type";;
esac

[ -n "$PURGE_DUPS" ] || PURGE_DUPS="0"
case $PURGE_DUPS in
  0) RUN+="--l=\"0\" ";;
  1) RUN+="--l=\"1\" ";;
  2) RUN+="--l=\"2\" ";;
  3) RUN+="--l=\"3\" ";;
  *) error "purge dups level can not be set to $PURGE_DUPS. not a valid level";;
esac

[ -n "$KMER" ] || KMER="kmc"
case $KMER in
  kmc) RUN+="--kmer=\"kmc\" ";;
  jellyfish) RUN+="--kmer=\"jellyfish\" ";;
  *) error "K-mer counting tool $KMER does not valid";;
esac

[ -n "$HAPSCAFFOLD" ] && [ -z $POLISHTYPE ] && warn "running hapscaffold without any polishing results in no change, not going to run ragtag"

state "checking for running busco"
if [ -n "$BUSCO" -o -n "$LINEAGE" -o -n "$BUSCOPATH" ]; then
state "running busco, checking busco things"
  [ -z "$LINEAGE" -a -z "$BUSCOPATH" ] && error "trying to setup busco, but no lineage/path/file given"
  if [ -n "$LINEAGE" ]; then
    state "   ...busco lineage described: ${LINEAGE}"
    BUSCOSTRING="--busco --linreage=\"${LINEAGE}\" "
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

RUN+="--outfasta=\"${NAME}.genome.out\" "
[ -n "$THREADS" ] && RUN+="--threads=\"$THREADS\" " || warn "threads not set, setting to 20 maximum threads"
[ -z "$THREADS" ] && RUN+="--threads=\"20\" "
[ -f "$R1" ] && RUN+="--hicreadf=\"$( readlink -f $R1 )\" "
[ -f "$R2" ] && RUN+="--hicreadr=\"$( readlink -f $R2 )\" "
[ -f "$MATF" ] && RUN+="--matreadf=\"$( readlink -f $MATF )\" "
[ -f "$MATR" ] && RUN+="--matreadr=\"$( readlink -f $MATR )\" "
[ -f "$PATF" ] && RUN+="--patreadf=\"$( readlink -f $PATF )\" "
[ -f "$PATR" ] && RUN+="--patreadr=\"$( readlink -f $PATR )\" "
[ -n "$YAHS" ] && RUN+="--scaffold --yahs "
[ -n "$HAPSCAFFOLD" ] && RUN+="--hapscaffold "
[ -n "$BUSCO" ] && RUN+="$BUSCO " || state "not running busco"
[ -n "$BUSCO" ] && [ -n "$BUSCOEVALUE" ] && RUN+=" --buscoevalue=\"BUSCOEVALUE\" "
[ -n "$POLISHTYPE" ] && RUN+="--polish --polishtype=\"$POLISHTYPE\" " || warn "not polishing, it is recomended that you polish"
[ -n "$BUSCOSTRING" ] && RUN+="$BUSCOSTRING"
[ -n "$READS" ] && RUN+="--readin=\"${READS}\" " || error "reads file(s) not given, exiting"
[ -z "$NAME" ] && NAME="$(date +%s)" && state "name not given, setting name to: $NAME"
[ -n "$LITE" ] && RUN+="--lite "
RUN+="--assembly=\"$NAME\" "
[ -z "$SUPRESS" ] && RUN+="-bg"

pizzaz "$RUN"
[ -z "$SUPRESS" ] && stop_check "check that the command is expected, continue"

state "Prefetching singularity containers"
[ -n "$NXF_SINGULARITY_LIBRARYDIR" ] && "Nextflow Singularity Library directory set: $NXF_SINGULARITY_LIBRARYDIR, will use for singularity images" || warn "NXF_SINGULARITY_LIBRARYDIR not set"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] && "Nextflow Singularity cache directory set: $NXF_SINGULARITY_CACHEDIR" || warn "NXF_SINGULARITY_CACHEDIR not set"

prefetch_container="./scr/prefetch_containers.sh"
[ -n "$YAHS" ] && prefetch_container+=" -y"
[ -n "$BUSCO" ] && prefetch_container+=" -b"
[ -n "$POLISHTYPE" ] && prefetch_container+=" -p $POLISHTYPE"
[ -n "$NXF_SINGULARITY_CACHEDIR" ] || ( mkdir -p "./work/singularity"; prefetch_container+=" -l ./work/singularity" )
eval $prefetch_container

if [ -n "$TEST" ]; then
  check_container="./scr/check_containers.sh"
  [ -n "$YAHS" ] && check_container+=" -y"
  [ -n "$BUSCO" ] && check_container+=" -b"
  [ -n "$POLISHTYPE" ] && check_container+=" -p $POLISHTYPE"
  [ -n "$NXF_SINGULARITY_CACHEDIR" ] || check_container+=" -l ./work/singularity"
  eval $check_container
fi

[ -z "$SUPRESS" ] && stop_check "proceed with run"
state "making /reports dir"
mkdir -p reports
pizzaz "running only the best"
echo "${RUN}" > "./reports/${NAME}.nextflow.command.txt"
echo "${RUN}" > "./reports/nextflow-${NAME}.stdout.txt"
eval "${RUN}" &> "./reports/nextflow-${NAME}.stdout.txt"
