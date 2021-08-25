#!/bin/bash


location=$( pwd )
echo "fetch location is:"
echo $NXF_SINGULARITY_CACHEDIR 

cd $NXF_SINGULARITY_CACHEDIR

singularity pull --library https://library.sylabs.io bryce911-bbtools.img library://molikd/remote-builds/rb-6126cc70d78bf1ce9c212511:latest
singularity pull --library https://library.sylabs.io dmolik-genomescope2.img library://molikd/remote-builds/rb-6126cc70d78bf1ce9c212512:latest
singularity pull --library https://library.sylabs.io dmolik-hifiasm.img library://molikd/remote-builds/rb-6126cc70a267b3f6a3832aca:latest
singularity pull --library https://library.sylabs.io dmolik-jellyfish.img library://molikd/remote-builds/rb-6126cc71a267b3f6a3832acb:latest
singularity pull --library https://library.sylabs.io dmolik-pbadapterfilt.img library://molikd/remote-builds/rb-6126cc71596e300e6bbddba2:latest
singularity pull --library https://library.sylabs.io dmolik-ragtag.img library://molikd/remote-builds/rb-6126cc71d78bf1ce9c212513:latest
singularity pull --library https://library.sylabs.io mgibio-samtools:1.9.img library://molikd/remote-builds/rb-6126cc72a267b3f6a3832acc:latest
singularity pull --library https://library.sylabs.io dmolik-shhquis.img library://molikd/remote-builds/rb-6126cc72596e300e6bbddba3:latest
singularity pull --library https://library.sylabs.io ezlabgva-busco:v5.2.2_cv1.img library://molikd/remote-builds/rb-6126cc72d78bf1ce9c212514:latest
singularity pull --library https://library.sylabs.io koszullab-hicstuff.img library://molikd/remote-builds/rb-6126cc73a267b3f6a3832acd:latest
singularity pull --library https://library.sylabs.io pvstodghill-any2fasta.img library://molikd/remote-builds/rb-6126cc73d78bf1ce9c212515:latest

cd $location
