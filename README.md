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

[![](https://mermaid.ink/img/pako:eNqNVsFunDAQ_RWLE0hJKuWYQ6tKadVLT-ltN0JeGBY3YBPbZLPa7L93xjZgdqFtlA3jmfeen4eB7CkpVAnJQ7LXvKvZr8etZPhTKgkb-nP37BPK7ja9Ac0K3jSGlnemfma3t58ZyLcUP0Jp2YK0jMuSVaIBw4oaihcoMy-BGIc3Qu77hmthjyeugUk4xDlWKGm5kKAN0_DaCw3ll7OXiGEk9XEEQ5-cCmmpDrJRvPyLYraiIxX-epUVK_AujA30YUfHlfBuq0YdNkPABBrvZWhckJ1BQ2XAU8m1Kq-4sTx1MaP49ZPLZAuENzQllExLMIUWO2BGVfZA_QwVk61us-Nt2AQjf6-WwFjMsWhRL_XXK3wEcZQf4rv4-vQzpSkJMbWDluMcREd1nN_QNMdKmDodo9DCQJjSBN-DVC0eWnWQRvGcMuxNhFoUxvZVlQ7BHLqV4Sg96pyw5AK2V3TXrGLYWgQPEzhquAYRMmJPM4mM3KVSjHwxuwL6oZugpcJpsSxmRJz7dW_sHONGcXp-rwszeznf43xPJplbB6sXKEf3muEx4LUhU3T9r34RcKLO7VDG-aAgu0RNrXK4qFMBPhHuVy0NbXKohS5F-dhV1CMne9miCXPVIWHyTjU4ueTJR2QpNkbwxWaN3AupuUHN95bvnTkf3nXHaafIpa86Mi1NXeOr1Xc8xBE4ZEa0J-dK5jXvos2YkgwzjbLHDtYEDjW3-YuQZa6qoR2UY5Rjqgp26alcGaKLs483zVcX9D3OiLZr4MMPv8dWBHvDlzuX1mxoxYbV8z_UWtCI_5gp_INSAnQBeskbzHvVXOmcwCe_ZEqziDu0YSbhdpizPSp-J8bG_cXdPB9mSw4i2dXy7FxR7LSjdTYdc1dUVqnGz9uwoH-vBqTpzcyKf7UGTPg6MokuVEe9xXFdwMQTMSu4N0Byk6CRlosSvxadCLVNbA0tbJMHDCulwdhtspVnRPZdyS18K4VVOnmoeGPgJuG9VU9HWYwJj3oUHL9ltSF7_gNGp02z)](https://mermaid-js.github.io/mermaid-live-editor/edit/#pako:eNqNVsFunDAQ_RWLE0hJKuWYQ6tKadVLT-ltN0JeGBY3YBPbZLPa7L93xjZgdqFtlA3jmfeen4eB7CkpVAnJQ7LXvKvZr8etZPhTKgkb-nP37BPK7ja9Ac0K3jSGlnemfma3t58ZyLcUP0Jp2YK0jMuSVaIBw4oaihcoMy-BGIc3Qu77hmthjyeugUk4xDlWKGm5kKAN0_DaCw3ll7OXiGEk9XEEQ5-cCmmpDrJRvPyLYraiIxX-epUVK_AujA30YUfHlfBuq0YdNkPABBrvZWhckJ1BQ2XAU8m1Kq-4sTx1MaP49ZPLZAuENzQllExLMIUWO2BGVfZA_QwVk61us-Nt2AQjf6-WwFjMsWhRL_XXK3wEcZQf4rv4-vQzpSkJMbWDluMcREd1nN_QNMdKmDodo9DCQJjSBN-DVC0eWnWQRvGcMuxNhFoUxvZVlQ7BHLqV4Sg96pyw5AK2V3TXrGLYWgQPEzhquAYRMmJPM4mM3KVSjHwxuwL6oZugpcJpsSxmRJz7dW_sHONGcXp-rwszeznf43xPJplbB6sXKEf3muEx4LUhU3T9r34RcKLO7VDG-aAgu0RNrXK4qFMBPhHuVy0NbXKohS5F-dhV1CMne9miCXPVIWHyTjU4ueTJR2QpNkbwxWaN3AupuUHN95bvnTkf3nXHaafIpa86Mi1NXeOr1Xc8xBE4ZEa0J-dK5jXvos2YkgwzjbLHDtYEDjW3-YuQZa6qoR2UY5Rjqgp26alcGaKLs483zVcX9D3OiLZr4MMPv8dWBHvDlzuX1mxoxYbV8z_UWtCI_5gp_INSAnQBeskbzHvVXOmcwCe_ZEqziDu0YSbhdpizPSp-J8bG_cXdPB9mSw4i2dXy7FxR7LSjdTYdc1dUVqnGz9uwoH-vBqTpzcyKf7UGTPg6MokuVEe9xXFdwMQTMSu4N0Byk6CRlosSvxadCLVNbA0tbJMHDCulwdhtspVnRPZdyS18K4VVOnmoeGPgJuG9VU9HWYwJj3oUHL9ltSF7_gNGp02z)

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
