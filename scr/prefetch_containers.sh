#!/bin/bash

if [ command -v io.sh &> /dev/null ] || [ -f "io.sh" ]; then
  source io.sh
elif [ -f "scr/io.sh" ]; then
  source scr/io.sh
else
  echo >&2 "[$(date)] I require io.sh to run, but io.sh has not been found, not downloading, aborting"; exit 1;
fi

state "PREFETCH CONTAINERS"

while [ $# -gt 0 ] ; do
  case $1 in
      -v | --version) version ;;
      -p | --polish) POLISHTYPE="$2" ;;
      -b | --busco) BUSCO="true";;
      -l | --location) LOCATION="$2" ;;
      -y | --yahs) YAHS="true";;
  esac
  shift
done

location=$( pwd )
describe "fetch location is:"
if [ -n "$NXF_SINGULARITY_LIBRARYDIR" ]; then
 pizzaz $NXF_SINGULARITY_LIBRARYDIR
 cd $NXF_SINGULARITY_LIBRARYDIR
elif [ -n "$LOCATION" ]; then
 pizzaz $LOCATION
 cd $LOCATION
else
 error "..not set, please set NXF_SINGULARITY_LIBRARYDIR or give me a location"
fi

[ ! -f dmolik-gfastats.img ] && singularity pull dmolik-gfastats.img docker://dmolik/gfastats || state "   ...gfastats container found, not downloading"
[ ! -f dmolik-genomescope2.img ] && singularity pull dmolik-genomescope2.img docker://dmolik/genomescope2 || state "   ...genomescope2 container found, not downloading"
[ ! -f dmolik-hifiasm.img ] && singularity pull dmolik-hifiasm.img docker://dmolik/hifiasm || state "   ...hifiasm container found, not downloading"
[ ! -f dmolik-k-mer-counting-tools.img ] && singularity pull dmolik-k-mer-counting-tools.img docker://dmolik/k-mer-counting-tools || state "   ...k-mer counting container found, not downloading"
[ ! -f dmolik-pbadapterfilt.img ] && singularity pull dmolik-pbadapterfilt.img docker://dmolik/pbadapterfilt || state "   ...HiFi Filter container found, not downloading"
[ ! -f mgibio-samtools-1.9.img ] && singularity pull mgibio-samtools-1.9.img docker://mgibio/samtools:1.9 || state "   ...samtools container found, not downloading"
[ ! -f koszullab-hicstuff.img ] && singularity pull koszullab-hicstuff.img docker://koszullab/hicstuff || state "   ...hicstuff container found, not downloading"
[ ! -f pvstodghill-any2fasta.img ] && singularity pull pvstodghill-any2fasta.img docker://pvstodghill/any2fasta || state "   ...any2fasta container found, not downloading"

if [ -n "$BUSCO" ]; then
  state "need busco container for this run"
  [ ! -f ezlabgva-busco-v5.2.2_cv1.img ] && singularity pull ezlabgva-busco-v5.2.2_cv1.img docker://ezlabgva/busco:v5.2.2_cv1 || state "   ...busco container found, not downloading"
fi

if [ -n "$YAHS" ]; then
  state "need yahs container for this run"
  [ ! -f dceoy-bwa-mem2.img ] && singularity pull dceoy-bwa-mem2.img docker://dceoy/bwa-mem2 || state "bwa mem container found, not downloading"
  [ ! -f dmolik-yahs.img ] && singularity pull dmolik-yahs.img docker://dmolik/yahs || state "yahs container found, not downloading"
  [ ! -f dmolik-juicer-tools.img ] && singularity pull dmolik-juicer-tools.img docker://dmolik/juicer-tools || state "juicer tools container found, not downloading"
fi

case $POLISHTYPE in
  "simple")
    state "need minimum polishing containers for this run"
    [ ! -f dmolik-ragtag.img ] && singularity pull dmolik-ragtag.img docker://dmolik/ragtag || state "   ...ragtag.py container found, not downloading"
    [ ! -f dmolik-shhquis.img ] && singularity pull dmolik-shhquis.img docker://dmolik/shhquis || state "   ...shhquis container found, not downloading"
  ;;
  "merfin")
    state "need merfin containers for this run"
    [ ! -f dmolik-ragtag.img ] && singularity pull dmolik-ragtag.img docker://dmolik/ragtag || state "   ...ragtag.py container found, not downloading"
    [ ! -f dmolik-shhquis.img ] && singularity pull dmolik-shhquis.img docker://dmolik/shhquis || state "   ...shhquis container found, not downloading"
    [ ! -f mgibio-bcftools-1.9.img ] && singularity pull mgibio-bcftools-1.9.img docker://mgibio/bcftools:1.9 || state "   ...bcftools container found, not downloading"
    [ ! -f dmolik-merfin.img ] && singularity pull dmolik-merfin.img docker://dmolik/merfin || state "   ...merfin container found, not downloading"
  ;;
  "dv")
    state "need google deep variant containers for this run"
    [ ! -f dmolik-ragtag.img ] && singularity pull dmolik-ragtag.img docker://dmolik/ragtag || state "   ...ragtag.py container found, not downloading"
    [ ! -f dmolik-shhquis.img ] && singularity pull dmolik-shhquis.img docker://dmolik/shhquis || state "   ...shhquis container found, not downloading"
    [ ! -f mgibio-bcftools-1.9.img ] && singularity pull mgibio-bcftools-1.9.img docker://mgibio/bcftools:1.9 || state "   ...bcftools container found, not downloading"
    [ ! -f google-deepvariant.img ] && singularity pull google-deepvariant.img docker://google/deepvariant || state "   ...deepvariant container found, not downloading"
  ;;
esac

cd $location
