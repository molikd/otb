#!/bin/bash

if [ command -v io.sh &> /dev/null ] || [ -f "io.sh" ]; then
  source io.sh
elif [ -f "scr/io.sh" ]; then
  source scr/io.sh
else
  echo >&2 "[$(date)] I require io.sh to run, but io.sh has not been found, aborting"; exit 1;
fi

state "PREFETCH CONTAINERS"

location=$( pwd )
describe "fetch location is:"
if [ -n "$NXF_SINGULARITY_CACHEDIR" ]; then
 pizzaz $NXF_SINGULARITY_CACHEDIR
 cd $NXF_SINGULARITY_CACHEDIR
elif [ -n $1 ]; then
 pizzaz $1
 cd $1
else
 error "..not set, please set NXF_SINGULARITY_CACHEDIR or give me a location"
fi

[ ! -f bryce911-bbtools.img ] && singularity pull bryce911-bbtools.img docker://bryce911/bbtools || state "   ...bbtools container found"
[ ! -f dmolik-genomescope2.img ] && singularity pull dmolik-genomescope2.img docker://dmolik/genomescope2 || state "   ...genomescope2 container found"
[ ! -f dmolik-hifiasm.img ] && singularity pull dmolik-hifiasm.img docker://dmolik/hifiasm || state "   ...hifiasm container found"
[ ! -f dmolik-jellyfish.img ] && singularity pull dmolik-jellyfish.img docker://dmolik/jellyfish || state "   ...jellyfish container found"
[ ! -f dmolik-pbadapterfilt.img ] && singularity pull dmolik-pbadapterfilt.img docker://dmolik/pbadapterfilt || state "   ...HiFi Filter container found"
[ ! -f dmolik-ragtag.img ] && singularity pull dmolik-ragtag.img docker://dmolik/ragtag || state "   ...ragtag.py container found"
[ ! -f mgibio-samtools-1.9.img ] && singularity pull mgibio-samtools:1.9.img docker://mgibio/samtools:1.9 || state "   ...samtools container found"
[ ! -f dmolik-shhquis.img ] && singularity pull dmolik-shhquis.img docker://dmolik/shhquis || state "   ...shhquis container found"
[ ! -f ezlabgva-busco-v5.2.2_cv1.img ] && singularity pull ezlabgva-busco:v5.2.2_cv1.img docker://ezlabgva/busco:v5.2.2_cv1 || state "   ...busco container found"
[ ! -f koszullab-hicstuff.img ] && singularity pull koszullab-hicstuff.img docker://koszullab/hicstuff || state "   ...hicstuff container found"
[ ! -f pvstodghill-any2fasta.img ] && singularity pull pvstodghill-any2fasta.img docker://pvstodghill/any2fasta || state "   ...any2fasta container found"
[ ! -f dmolik-merfin.img ] && singularity pull dmolik-merfin.img docker://dmolik/merfin || state "   ...merfin container found"
[ ! -f mgibio-bcftools-1.9.img ] && singularity pull mgibio-bcftools-1.9.img docker://mgibio/bcftools:1.9 || state "   ...bcftools container found"
[ ! -f google-deepvariant.img ] && singularity pull google-deepvariant.img docker://google/deepvariant || state "   ...deepvariant container found"

cd $location
