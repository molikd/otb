```
                                                                         
       .----------------. .----------------. .----------------.          
      | .--------------. | .--------------. | .--------------. |         
      | |     ____     | | |  _________   | | |   ______     | |         
      | |   .'    `.   | | | |  _   _  |  | | |  |_   _ \    | |         
      | |  /  .--.  \  | | | |_/ | | \_|  | | |    | |_) |   | |         
      | |  | |    | |  | | |     | |      | | |    |  __'.   | |         
      | |  \  `--'  /  | | |    _| |_     | | |   _| |__) |  | |         
      | |   `.____.'   | | |   |_____|    | | |  |_______/   | |         
      | |              | | |              | | |              | |         
      | '--------------' | '--------------' | '--------------' |         
       '----------------' '----------------' '----------------'          
      ____       __       ________         ___         __    __          
     / __ \___  / /_ __  /_  __/ /  ___   / _ )___ ___/ /_  / /          
    / /_/ / _ \/ / // /   / / / _ \/ -_) / _  / -_|_-< __/ /_/           
   _\____/_//_/_/\_, /   /_/ /_//_/\__/ /____/\__/___|__/_(_) ___        
  / ___/__ ___  /___/__ _  ___   / _ | ___ ______ __ _  / /  / (_)__ ___ 
 / (_ / -_) _ \/ _ \/  ' \/ -_) / __ |(_-<(_-< -_)  ' \/ _ \/ / / -_|_-< 
 \___/\__/_//_/\___/_/_/_/\__/ /_/ |_/___/___|__/_/_/_/_.__/_/_/\__/___/ 
                                                                         
```

*o*nly *t*he *b*est (genome assemblies) is a Hi-C / HiFi pipeline specifically designed for phasing

In order to utilize this pipeline:
  - singularity \(https://sylabs.io/singularity/) must be installed 
  - nextflow \(https://www.nextflow.io/) must be installed 
Most HPC support groups will either already have these tools installed, or be willing to install them as [modules](http://modules.sourceforge.net/)

otb operates from a local directory, and must be operated as ./otb.sh since it sources some shell scripts

in order to use otb download this directory and use ./otb.sh, an example:

```bash
 ./otb.sh --sge --mode homozygous --threads 40 -f otb_test_file_R2.fastq -r otb_test_file_R1.fastq --polish-type simple --bam otb_test.bam
```

otb is fully featured and utilzes the following softwares:
- bbtools
- genomescope2
- hifiasm
- jellyfish
- pbadapterfilt
- ragtag
- samtools
- shhquis
- busco
- hicstuff
- any2fasta
- blobtools

otb is in the public domain in the United States per 17 U.S.C. § 105
