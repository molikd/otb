#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p long
#SBATCH -o "otb.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "otb.stderr.%j.%N" 		#optional but it prints our standard error
#SBATCH --mail-user=<username>
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load nextflow/22.04.3
module load singularityCE

#sometimes HiC data comes out as bz2, otb likes gz so here is how you would convert:
#bzcat RawHiC/JNBN_HiC_NA_NA_GTGAAA_Bombus_huntii-Bombus_huntii_HiC_I1143_L4_R1.fastq.bz2 | gzip -c > RawHiC/JNBN_HiC_NA_NA_GTGAAA_Bombus_huntii-Bombus_huntii_HiC_I1143_L4_R1.fastq.gz
#bzcat RawHiC/JNBN_HiC_NA_NA_GTGAAA_Bombus_huntii-Bombus_huntii_HiC_I1143_L4_R2.fastq.bz2 | gzip -c > RawHiC/JNBN_HiC_NA_NA_GTGAAA_Bombus_huntii-Bombus_huntii_HiC_I1143_L4_R2.fastq.gz

# We need an assembly name, generally this is just the name of the organism"
# Assembly_Name="Bombus_huntii"
Assembly_Name=$(basename "$( pwd )")

# Foward and reverse reads
CCS='RawData/*.fastq'

#Comment/Uncommment for busco
Busco="--busco" #Busco will be run
Busco_Location="-l /project/ag100pest/software/OTB_test/busco"
Busco_DB="-p insecta_odb10"

#Busco="" #Busco will not be run

#Comment/Uncomment for Type (only select one of)
#HiFi_Type="phasing"
HiFi_Type="default"
#HiFi_Type="trio"

#Comment/Uncomment for Runner (only select one of)
#Runner="slurm_usda"
Runner="slurm_usda"

Threads="32"

if [[ ! -z "$BUSCO" ]]; then
  ./otb.sh --lite -n ${Assembly_Name}  -in "$(echo "${CCS}")" -m ${HiFi_Type} -t ${Threads} --runner ${Runner} -c -s
else
  ./otb.sh --lite -n ${Assembly_Name}  -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --runner ${Runner} ${Busco} ${Busco_Location} ${Busco_DB} -c -s
fi
