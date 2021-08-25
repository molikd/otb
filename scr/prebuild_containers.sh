#!/bin/bash


location=$( pwd )
echo "fetch location is:"
echo $NXF_SINGULARITY_CACHEDIR 

cd $NXF_SINGULARITY_CACHEDIR

singularity build --force --remote bryce911-bbtools.img docker://bryce911/bbtools
singularity build --force --remote dmolik-genomescope2.img docker://dmolik/genomescope2
singularity build --force --remote dmolik-hifiasm.img docker://dmolik/hifiasm
singularity build --force --remote dmolik-jellyfish.img docker://dmolik/jellyfish
singularity build --force --remote dmolik-pbadapterfilt.img docker://dmolik/pbadapterfilt
singularity build --force --remote dmolik-ragtag.img docker://dmolik/ragtag
singularity build --force --remote dmolik-samtools.img docker://dmolik/samtools
singularity build --force --remote dmolik-shhquis.img docker://dmolik/shhquis
singularity build --force --remote ezlabgva-busco:v5.2.2_cv1.img docker://ezlabgva/busco:v5.2.2_cv1
singularity build --force --remote koszullab-hicstuff.img docker://koszullab/hicstuff
singularity build --force --remote pvstodghill-any2fasta.img docker://pvstodghill/any2fasta

cd $location
