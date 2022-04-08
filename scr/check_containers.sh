#!/bin/bash

if [ command -v io.sh &> /dev/null ] || [ -f "io.sh" ]; then
  source io.sh
elif [ -f "scr/io.sh" ]; then
  source scr/io.sh
else
  echo >&2 "[$(date)] I require io.sh to run, but io.sh has not been found, not downloading, aborting"; exit 1;
fi

state "CHECKING CONTAINERS"

while [ $# -gt 0 ] ; do
  case $1 in
      -v | --version) version ;;
      -p | --polish) POLISHTYPE="$2" ;;
      -b | --busco) BUSCO="true";;
      -l | --location) LOCATION="$2" ;;
  esac
  shift
done

location=$( pwd )
describe "fetch location is:"
if [ -n "$NXF_SINGULARITY_CACHEDIR" ]; then
 pizzaz $NXF_SINGULARITY_CACHEDIR
 cd $NXF_SINGULARITY_CACHEDIR
elif [ -n "$LOCATION" ]; then
 pizzaz $LOCATION
 cd $LOCATION
else
 error "..not set, please set NXF_SINGULARITY_CACHEDIR or give me a location"
fi

singularity exec "bryce911-bbtools.img" echo "   ...hello from bbtools container" || error "bbtools container broken, exiting"
singularity exec "dmolik-genomescope2.img" echo "   ...hello from genomescope2 container" || error "genomescope2 container broken, exiting"
singularity exec "dmolik-hifiasm.img"  echo "   ...hello from hifiasm container" || error "hifiasm container broken, exiting"
singularity exec "dmolik-k-mer-counting-tools.img" echo "   ...hello from jellyfish container" || error "jellyfish container broken, exiting"
singularity exec "dmolik-pbadapterfilt.img" echo "   ...hello from pbadapterfilt container" || error "pbadapterfilt container broken, exiting"
singularity exec "mgibio-samtools-1.9.img" echo "   ...hello from samtools container" || error "samtools container broken, exiting"
singularity exec "koszullab-hicstuff.img" echo "   ...hello from hicstuff container" || error "hicstuff container broken, exiting"
singularity exec "pvstodghill-any2fasta.img" echo "   ...hello from any2fasta container" || error "any2fasta container borken, exiting"

if [ -n "$BUSCO" ]; then
  state "busco container will be required for this run, testing busco"
  singularity exec "ezlabgva-busco-v5.2.2_cv1.img" echo "   ...hello from busco container" || error "busco container broken, exiting"
fi

case $POLISHTYPE in
  "simple")
    state "this will be a simple polish, checking ragtag and shhquis"
    singularity exec "dmolik-ragtag.img" echo "   ...hello from ragtag container" || error "ragtag container broken, exiting"
    singularity exec "dmolik-shhquis.img" echo "    ...hello from shhquis container" || error "shhquis container broken, exiting"
  ;;
  "merfin")
    state "this will be a merfin polish, checking ragtag, shhquis, bcftools, and merfin"
    singularity exec "dmolik-ragtag.img" echo "   ...hello from ragtag container" || error "ragtag container broken, exiting"
    singularity exec "dmolik-shhquis.img" echo "    ...hello from shhquis container" || error "shhquis container broken, exiting"
    singularity exec "mgibio-bcftools-1.9.img" echo "   ...hello from bcftools container" || error "bcftools container broken, exiting"
    singularity exec "dmolik-merfin.img" echo "    ...hello from merfin container" || error "merfin container broken, exiting"
  ;;
  "dv")
    state "this will be a deep variant polish, checking ragtag, shhquis, bcftools, and merfin"
    singularity exec "dmolik-ragtag.img" echo "   ...hello from ragtag container" || error "ragtag container broken, exiting"
    singularity exec "dmolik-shhquis.img" echo "    ...hello from shhquis container" || error "shhquis container broken, exiting"
    singularity exec "mgibio-bcftools-1.9.img" echo "   ...hello from bcftools container" || error "bcftools container broken, exiting"
    singularity exec "google-deepvariant.img" echo "   ...hello from deepvariant container" || error "deepvariant container broken, exiting"
  ;;
esac

state "all required containers checked with an intial pass"

