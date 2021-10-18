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

singularity pull --disable-cache --force bryce911-bbtools.img docker://bryce911/bbtools
singularity pull --disable-cache --force dmolik-genomescope2.img docker://dmolik/genomescope2
singularity pull --disable-cache --force dmolik-hifiasm.img docker://dmolik/hifiasm
singularity pull --disable-cache --force dmolik-jellyfish.img docker://dmolik/jellyfish
singularity pull --disable-cache --force dmolik-pbadapterfilt.img docker://dmolik/pbadapterfilt
singularity pull --disable-cache --force dmolik-ragtag.img docker://dmolik/ragtag
singularity pull --disable-cache --force mgibio-samtools:1.9.img docker://mgibio/samtools:1.9
singularity pull --disable-cache --force dmolik-shhquis.img docker://dmolik/shhquis
singularity pull --disable-cache --force ezlabgva-busco:v5.2.2_cv1.img docker://ezlabgva/busco:v5.2.2_cv1
singularity pull --disable-cache --force koszullab-hicstuff.img docker://koszullab/hicstuff
singularity pull --disable-cache --force pvstodghill-any2fasta.img docker://pvstodghill/any2fasta
singularity pull --disable-cache --force dmolik-merfin.img docker://dmolik/merfin
singularity pull --disable-cache --force mgibio-bcftools-1.9.img docker://mgibio/bcftools:1.9
singularity pull --disable-cache --force google-deepvariant.img docker://google/deepvariant

cd $location
