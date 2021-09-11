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

[ ! -f bryce911-bbtools.img ] && singularity pull bryce911-bbtools.img docker://bryce911/bbtools
[ ! -f dmolik-genomescope2.img ] && singularity pull dmolik-genomescope2.img docker://dmolik/genomescope2
[ ! -f dmolik-hifiasm.img ] && singularity pull dmolik-hifiasm.img docker://dmolik/hifiasm
[ ! -f dmolik-jellyfish.img ] && singularity pull dmolik-jellyfish.img docker://dmolik/jellyfish
[ ! -f dmolik-pbadapterfilt.img ] && singularity pull dmolik-pbadapterfilt.img docker://dmolik/pbadapterfilt
[ ! -f dmolik-ragtag.img ] && singularity pull dmolik-ragtag.img docker://dmolik/ragtag
[ ! -f mgibio-samtools:1.9.img ] && singularity pull mgibio-samtools:1.9.img docker://mgibio/samtools:1.9
[ ! -f dmolik-shhquis.img ] && singularity pull dmolik-shhquis.img docker://dmolik/shhquis
[ ! -f ezlabgva-busco:v5.2.2_cv1.img ] && singularity pull ezlabgva-busco:v5.2.2_cv1.img docker://ezlabgva/busco:v5.2.2_cv1
[ ! -f koszullab-hicstuff.img ] && singularity pull koszullab-hicstuff.img docker://koszullab/hicstuff
[ ! -f pvstodghill-any2fasta.img ] && singularity pull pvstodghill-any2fasta.img docker://pvstodghill/any2fasta

cd $location
