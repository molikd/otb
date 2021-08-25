#!/bin/bash


location=$( pwd )
echo "fetch location is:"
echo $NXF_SINGULARITY_CACHEDIR 

cd $NXF_SINGULARITY_CACHEDIR

singularity pull bryce911-bbtools.img library://molikd/remote-builds/rb-6126cc70d78bf1ce9c212511:latest
singularity pull dmolik-genomescope2.img library://molikd/remote-builds/rb-6126cc70d78bf1ce9c212512:latest
singularity pull dmolik-hifiasm.img library://molikd/remote-builds/rb-6126cc70a267b3f6a3832aca:latest
singularity pull dmolik-jellyfish.img library://molikd/remote-builds/rb-6126cc71a267b3f6a3832acb:latest
singularity pull dmolik-pbadapterfilt.img library://molikd/remote-builds/rb-6126cc71596e300e6bbddba2:latest
singularity pull dmolik-ragtag.img library://molikd/remote-builds/rb-6126cc71d78bf1ce9c212513:latest
singularity pull mgibio-samtools:1.9.img library://molikd/remote-builds/rb-6126cc72a267b3f6a3832acc:latest
singularity pull dmolik-shhquis.img library://molikd/remote-builds/rb-6126cc72596e300e6bbddba3:latest
singularity pull ezlabgva-busco:v5.2.2_cv1.img library://molikd/remote-builds/rb-6126cc72d78bf1ce9c212514:latest
singularity pull koszullab-hicstuff.img library://molikd/remote-builds/rb-6126cc73a267b3f6a3832acd:latest
singularity pull pvstodghill-any2fasta.img library://molikd/remote-builds/rb-6126cc73d78bf1ce9c212515:latest

cd $location
