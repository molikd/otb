```
 ------------------------------------------------------------------------------------- 
                                              ,,                                       
                                       mm    *MM                                       
                                       MM     MM                                       
                            ,pW"Wq.  mmMMmm   MM,dMMb.                                 
                           6W'   `Wb   MM     MM    `Mb                                
                           8M     M8   MM     MM     M8                                
                           YA.   ,A9   MM     MM.   ,M9                                
                            `Ybmd9'    `Mbmo  P^YbmdP'                                 
                                                                                       
 ------------------------------------------------------------------------------------- 
               _   ._   |        _|_  |_    _     |_    _    _  _|_  |                 
              (_)  | |  |  \/     |_  | |  (/_    |_)  (/_  _>   |_  o                 
                           /                                                           
  /   _    _   ._    _   ._ _    _      _.   _   _   _   ._ _   |_   |  o   _    _  \  
 |   (_|  (/_  | |  (_)  | | |  (/_    (_|  _>  _>  (/_  | | |  |_)  |  |  (/_  _>   | 
  \   _|                                                                            /  
 ------------------------------------------------------------------------------------- 
                                                                         
```
*o*nly *t*he *b*est (genome assemblies) is a Hi-C / HiFi pipeline specifically designed for phasing

Checkout the [wiki](https://github.com/molikd/otb/wiki), and [tutorial](https://github.com/molikd/otb/wiki/Tutorial)

In order to utilize this pipeline:
  - singularity \(https://sylabs.io/singularity/) must be installed 
  - nextflow \(https://www.nextflow.io/) must be installed 
Most HPC support groups will either already have these tools installed, or be willing to install them as [modules](http://modules.sourceforge.net/)

otb operates from a local directory, and must be operated as ./otb.sh since it sources some shell scripts

in order to use otb download this directory and use ./otb.sh, an example:

```bash
 ./otb.sh --sge --mode homozygous --threads 40 -f otb_test_file_R2.fastq -r otb_test_file_R1.fastq --polish-type simple --bam otb_test.bam
```

otb runs in the following fashion:

[![](https://mermaid.ink/img/pako:eNqNVsFu2zAM_RXBlzlA2wE99rBhQDfsslN3SwpDselYqy25ktzUSPvvIyXZlhN7XdHEFPneE0VRik9JrgpI7pKD5m3Fft_vJMO_QknY0tfNo3cou992BjTLeV0bGt6Y6pFdX39hIF9S_AilZQPSMi4LVooaDMsryJ-g2HgJxDi8EfLQ1VwL25-4BibhGPtYrqTlQoI2TMNzJzQUX9-9RAwjqbceDH0yCqSFOspa8eIfipsVHanw36uspAKvwthAH2Z0XAmvtqzVcTsYTGDinQyFC7IzaIgMeAq5UmUlN5anzmZkP392ns0C4QWTEkqmBZhciz0wo0p7pHqGiNmsTrPnTZgELb9XS2AMZhi0qJf65wU-gjjKT_FDfHv4lYYnlYIaZuyBaJkO_wfqui-FqdLRCuULhMlN8ANI1eCCVQtpZM8pw9xEqERubFeW6WDMoTsZltGhzglDzmAHRTtmFcOyInjovlHDFYeQEXvqR2RkzpWi5YObC6BvuAlaKOwUy2JGxLldz429x7hRnM7uZWCWXsYP2NtTksyNQ6pnKEf3muEI8MpQUvT8r3oRcKLO0yGPy4OMzTlqKpXDRZUK8Ilwu5rSUCaHWqhS5I-zimoU5TaPOlbD-z1kylagM99xzuPSdF7nvGJH-IRHNO-0xquy7rHyoi4oVZe16oYr5kLvYgeEyVpV48mgNXuLdOKFE3xxM0bumdS8AJofLD-4xXvzpu2nmaJa-Kgj09BUFV7bfkeDHYGDZ0R7cqZkVvE2mowpydBTK9u3sCZwrLjNnoQsMlUO5SAfIx9TZUiXTv1Kk56tfWwKH13Q9zgjmraGN3-4PLYk2Av-cHBpzZZGbBg9fqDWgEb820zhA0oB0AboOW9I3qtmSmcEPvkhU5pF3KEMMwnf0TO2R8V3bpy4f7jN8-ZmKYNIdjU8W1dkO-1ovJmWuc9Lq1Tt-20Y0E-3AWk6M0vFX90BE151JtGF6Ki32K4LmLgjZgF3wyRXCSbScFHgK9eJULsEz3gDu-QOzVJpMHaX7OQ7Iru24Ba-F8IqndyVvDZwlfDOqode5qPDo-4Fxze4Jnjf_wJ6329s)](https://mermaid-js.github.io/mermaid-live-editor/edit/#pako:eNqNVsFu2zAM_RXBlzlA2wE99rBhQDfsslN3SwpDselYqy25ktzUSPvvIyXZlhN7XdHEFPneE0VRik9JrgpI7pKD5m3Fft_vJMO_QknY0tfNo3cou992BjTLeV0bGt6Y6pFdX39hIF9S_AilZQPSMi4LVooaDMsryJ-g2HgJxDi8EfLQ1VwL25-4BibhGPtYrqTlQoI2TMNzJzQUX9-9RAwjqbceDH0yCqSFOspa8eIfipsVHanw36uspAKvwthAH2Z0XAmvtqzVcTsYTGDinQyFC7IzaIgMeAq5UmUlN5anzmZkP392ns0C4QWTEkqmBZhciz0wo0p7pHqGiNmsTrPnTZgELb9XS2AMZhi0qJf65wU-gjjKT_FDfHv4lYYnlYIaZuyBaJkO_wfqui-FqdLRCuULhMlN8ANI1eCCVQtpZM8pw9xEqERubFeW6WDMoTsZltGhzglDzmAHRTtmFcOyInjovlHDFYeQEXvqR2RkzpWi5YObC6BvuAlaKOwUy2JGxLldz429x7hRnM7uZWCWXsYP2NtTksyNQ6pnKEf3muEI8MpQUvT8r3oRcKLO0yGPy4OMzTlqKpXDRZUK8Ilwu5rSUCaHWqhS5I-zimoU5TaPOlbD-z1kylagM99xzuPSdF7nvGJH-IRHNO-0xquy7rHyoi4oVZe16oYr5kLvYgeEyVpV48mgNXuLdOKFE3xxM0bumdS8AJofLD-4xXvzpu2nmaJa-Kgj09BUFV7bfkeDHYGDZ0R7cqZkVvE2mowpydBTK9u3sCZwrLjNnoQsMlUO5SAfIx9TZUiXTv1Kk56tfWwKH13Q9zgjmraGN3-4PLYk2Av-cHBpzZZGbBg9fqDWgEb820zhA0oB0AboOW9I3qtmSmcEPvkhU5pF3KEMMwnf0TO2R8V3bpy4f7jN8-ZmKYNIdjU8W1dkO-1ovJmWuc9Lq1Tt-20Y0E-3AWk6M0vFX90BE151JtGF6Ki32K4LmLgjZgF3wyRXCSbScFHgK9eJULsEz3gDu-QOzVJpMHaX7OQ7Iru24Ba-F8IqndyVvDZwlfDOqode5qPDo-4Fxze4Jnjf_wJ6329s)

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
- merfin
- deep variant
- yahs
- bcftools
```
otb is in the public domain in the United States per 17 U.S.C. § 105
```
