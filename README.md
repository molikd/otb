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

otb operates from a local directory, and must be ran as ./otb.sh since it sources some shell scripts in the scr directory. 

in order to use otb download this repository and use ./otb.sh, an example:

```bash
 ./otb.sh --sge --mode homozygous --threads 40 -f otb_test_file_R2.fastq -r otb_test_file_R1.fastq --polish-type simple --bam otb_test.bam
```

otb runs in the following fashion:

[![](https://mermaid.ink/img/pako:eNqFVsFu2zAM_RXBlzlA2wE79rBhQDfsslN3SwpDselYqy25ktzUSPPvIyXZlhN7LZqYIh-fSIpUfEpyVUBynxw0byv252EnGf4VSsKWvu6evELZ_bYzoFnO69rQ8s5UT-z29isD-ZriRygtG5CWcVmwUtRgWF5B_gzFxlMgxuGNkIeu5lrY_sQ1MAnHWMdyJS0XErRhGl46oaH4dvYUMYyo3nsw9MnIkBbqKGvFi_8wblZ4pMJ_z7ISCrwJY4P7sKPzlfBmy1odt4PABAbeyVC4QDuDBsuAJ5MrVVZyY3nqZEbyy2en2Sw4vGJQQsm0AJNrsQdmVGmPVM9gMZvVbfa8CZug5M9qCYzGDI0W-VL_vMJHEOfyS_wU3x9_p-FJpaCGGXsgStPh_0Jd96UwVTpKoXzBYVIT_ABSNZiwaiGN5LnLsDc5VCI3tivLdBDm0J0MaXTIc0KTE9hB0YlZxbCsCB66b-RwxSFk5D31I3pkTpWi5I2bK6BvuAlaKOwUy2KPyOfLemzsHONGcprda8MsvIwfsLenIJlbh1AvUM7dc4YR4JVxUZGwGpRDLcQU6YeQSBVFRMsolMnqvBre7yFTtgKd-fN1Glc-p3XKG3aETzgQeac1Xkx1j3mKuqBQXdSqGwb6iu8qX2GyVtXYh5Szl4gnTpzgi60y-l5QzQug-cHyg0vei3dtP-0U1cJbnTMtTVXhJWmcY5AjcNCMaO-cKZlVvI02Y0oy1NTK9i2sERwrbrNnIYtMlUM5SMdIx1QZwqUZWxmhi9yxKegk4qZa7amopa4757JnJtQ0aA4XzdlCi0XHHu7uBYu3LdTCb2hE09bw7sfOY0uCveJPCpfWbGnFhtXTB2wNaMS_zxg-cCkA2gC99BuC96yZ0hmBT37JlGaR71DwGYWfvpm3R8W3cRy4f7jj8eJmKYKIdtU8yyuSHXe03kxp7vPSKlX77hgW9KNuQJrOzELxl3rAhJegiXTBOvItjtYCJu6ImcHdhslNgoE0XBT4MnYi1C7B-6iBXXKPYqk0GLtLdvKMyK4tuIUfhbBKJ_clrw3cJLyz6rGX-ajwqAfB8d2uCdrzP9wyd-Y)](https://mermaid-js.github.io/mermaid-live-editor/edit/#pako:eNqFVsFu2zAM_RXBlzlA2wE79rBhQDfsslN3SwpDselYqy25ktzUSPPvIyXZlhN7LZqYIh-fSIpUfEpyVUBynxw0byv252EnGf4VSsKWvu6evELZ_bYzoFnO69rQ8s5UT-z29isD-ZriRygtG5CWcVmwUtRgWF5B_gzFxlMgxuGNkIeu5lrY_sQ1MAnHWMdyJS0XErRhGl46oaH4dvYUMYyo3nsw9MnIkBbqKGvFi_8wblZ4pMJ_z7ISCrwJY4P7sKPzlfBmy1odt4PABAbeyVC4QDuDBsuAJ5MrVVZyY3nqZEbyy2en2Sw4vGJQQsm0AJNrsQdmVGmPVM9gMZvVbfa8CZug5M9qCYzGDI0W-VL_vMJHEOfyS_wU3x9_p-FJpaCGGXsgStPh_0Jd96UwVTpKoXzBYVIT_ABSNZiwaiGN5LnLsDc5VCI3tivLdBDm0J0MaXTIc0KTE9hB0YlZxbCsCB66b-RwxSFk5D31I3pkTpWi5I2bK6BvuAlaKOwUy2KPyOfLemzsHONGcprda8MsvIwfsLenIJlbh1AvUM7dc4YR4JVxUZGwGpRDLcQU6YeQSBVFRMsolMnqvBre7yFTtgKd-fN1Glc-p3XKG3aETzgQeac1Xkx1j3mKuqBQXdSqGwb6iu8qX2GyVtXYh5Szl4gnTpzgi60y-l5QzQug-cHyg0vei3dtP-0U1cJbnTMtTVXhJWmcY5AjcNCMaO-cKZlVvI02Y0oy1NTK9i2sERwrbrNnIYtMlUM5SMdIx1QZwqUZWxmhi9yxKegk4qZa7amopa4757JnJtQ0aA4XzdlCi0XHHu7uBYu3LdTCb2hE09bw7sfOY0uCveJPCpfWbGnFhtXTB2wNaMS_zxg-cCkA2gC99BuC96yZ0hmBT37JlGaR71DwGYWfvpm3R8W3cRy4f7jj8eJmKYKIdtU8yyuSHXe03kxp7vPSKlX77hgW9KNuQJrOzELxl3rAhJegiXTBOvItjtYCJu6ImcHdhslNgoE0XBT4MnYi1C7B-6iBXXKPYqk0GLtLdvKMyK4tuIUfhbBKJ_clrw3cJLyz6rGX-ajwqAfB8d2uCdrzP9wyd-Y)

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
