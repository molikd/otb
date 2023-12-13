```
 +-------------------------------------------------------------------------------------+
 |                                             ,,                                      |
 |                                      mm    *MM                                      |
 |                                      MM     MM                                      |
 |                           ,pW"Wq.  mmMMmm   MM,dMMb.                                |
 |                          6W'   `Wb   MM     MM    `Mb                               |
 |                          8M     M8   MM     MM     M8                               |
 |                          YA.   ,A9   MM     MM.   ,M9                               |
 |                           `Ybmd9'    `Mbmo  P^YbmdP'                                |
 |                                                                                     |
 |-------------------------------------------------------------------------------------|
 |              _   ._   |        _|_  |_    _     |_    _    _  _|_  |                |
 |             (_)  | |  |  \/     |_  | |  (/_    |_)  (/_  _>   |_  o                |
 |                          /                                                          |
 | /   _    _   ._    _   ._ _    _      _.   _   _   _   ._ _   |_   |  o   _    _  \ |
 ||   (_|  (/_  | |  (_)  | | |  (/_    (_|  _>  _>  (/_  | | |  |_)  |  |  (/_  _>   ||
 | \   _|                                                                            / |
 +-------------------------------------------------------------------------------------+

```
*o*nly *t*he *b*est (genome assemblies) is a Hi-C / HiFi pipeline specifically designed for phasing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6689817.svg)](https://doi.org/10.5281/zenodo.6689817) [![License: Public Domain](https://img.shields.io/badge/License-Public%20Domain-green)](https://raw.githubusercontent.com/molikd/otb/1/LICENSE)
---

Checkout the [wiki](https://github.com/molikd/otb/wiki), and [tutorial](https://github.com/molikd/otb/wiki/Tutorial)

In order to utilize this pipeline:
  - singularity \(https://sylabs.io/singularity/) must be installed
  - nextflow \(https://www.nextflow.io/) must be installed
Most HPC support groups will either already have these tools installed, or be willing to install them as [modules](http://modules.sourceforge.net/)

otb operates from a local directory, and must be ran as ./otb.sh since it sources some shell scripts in the scr directory.

we recomend that you set `NXF_SINGULARITY_LIBRARYDIR` in your bashrc (or similar CLI) environment, otb containers will be stored at that location.

in order to use otb download this repository and use ./otb.sh, an example:

```bash
 ./otb.sh --runner sge --mode homozygous --threads 40 -f otb_test_file_R2.fastq -r otb_test_file_R1.fastq --polish-type simple --reads otb_test.bam
```

otb runs in the following fashion:

[![](https://mermaid.ink/img/pako:eNqlVktv2zAM_iuCTw7QdkCPPWwYsA277NTdksJQbDrWakuuJDc10v73US9bTqwW2Iom4eMjRVKfZZ2yUlSQ3WUHSfuG_P624wT_KsFha75uHpxB6P12UCBJSdtWGfVGNQ_k-vozAf6c44cJyTvgmlBekZq1oEjZQPkI1calQIzFK8YPQ0sl0-OJSiAcjrGNlIJryjhIRSQ8DUxC9eXNpYhhJtXrCMp8CuPIK3HkraDVOxk3iTxc4L_LkigFXpjSPjysaGNtk1bdWpHohupkQwNHVUoodTv60fqFz5I536xbN4cXXbfiuA0CYTZnyBSsc6qaKk1zV5iRnz5Zy2Yl4BkLZILnFahSsj0QJWp9NDvkPWqTXGZPO78ISm7318DoLNCpMV_ufi_wEcSG_GQ_2Nf7X7n_NQ0bCk6sitq0-D_QtmPNVJNPkh-SD5jNBn4ALjpsWPSQR_IyJKxtAhpWKj3UdR6EJXTHfRsD5jmhywrkIMwWakFwrAgOfJ5y2OEYZBQ9MxwjCmvKDX-stLkAOgrP0EogszSJI6KY23Rt5C3GTcnNaXDpWJRX0APyfC6SWJ0ITnrR4sShImLQ_RCeo7NAm9Et458z2ihbqBGSdVrUSpmRPVRpTFGRRo1Kmb02qqPjHgqhG5CF3aXkyFyfYVMv4lZqS2BWppmeqR3Ke_NMTJWpwm2H6cdJppu4JwNf5egUe5ZqWb6kB00PtmYn3vTjvFJUrPPaYKOqpsHzXtlAL0dgb5nQtKK9tlgrCenOgijC2ie8W6wQvGhoHxVnZomWVuixh9SCRzzUi0fGq0LUYXzGRoyNiNq3Zw6DxLN-Niukg9nAmOpJpkdEv-TzOZNn1IJxAb5ObUeVfyV4RLQPaB4jk2T_b8qvEN-hV3bR1aJY17fw6k42h60N7Bnf4JRrtTUaCdrDB9k6kIh_XWT4IKQC6D30PC4U77IWQhYGfHIqQdpHsWGLFincli-iHSp-4cWFux-7AU7crFUQpU26F31Fss0d6Zu5zX1ZayFax-ugmDuUAq6G6ZoQHPMdxD7tt8kDwc_Lvml9sL_rzmWseOPc60dJAhdz6cJp303ZVYZFdZRVeP8-GeQuw0elg112h2ItJCi9y3b8DZFDX1EN3yuGfWV3NW0VXGV00OJ-5OVkcKhvjOJ1vvPWt79vNCyB)](https://mermaid-js.github.io/mermaid-live-editor/edit/#pako:eNqlVktv2zAM_iuCTw7QdkCPPWwYsA277NTdksJQbDrWakuuJDc10v73US9bTqwW2Iom4eMjRVKfZZ2yUlSQ3WUHSfuG_P624wT_KsFha75uHpxB6P12UCBJSdtWGfVGNQ_k-vozAf6c44cJyTvgmlBekZq1oEjZQPkI1calQIzFK8YPQ0sl0-OJSiAcjrGNlIJryjhIRSQ8DUxC9eXNpYhhJtXrCMp8CuPIK3HkraDVOxk3iTxc4L_LkigFXpjSPjysaGNtk1bdWpHohupkQwNHVUoodTv60fqFz5I536xbN4cXXbfiuA0CYTZnyBSsc6qaKk1zV5iRnz5Zy2Yl4BkLZILnFahSsj0QJWp9NDvkPWqTXGZPO78ISm7318DoLNCpMV_ufi_wEcSG_GQ_2Nf7X7n_NQ0bCk6sitq0-D_QtmPNVJNPkh-SD5jNBn4ALjpsWPSQR_IyJKxtAhpWKj3UdR6EJXTHfRsD5jmhywrkIMwWakFwrAgOfJ5y2OEYZBQ9MxwjCmvKDX-stLkAOgrP0EogszSJI6KY23Rt5C3GTcnNaXDpWJRX0APyfC6SWJ0ITnrR4sShImLQ_RCeo7NAm9Et458z2ihbqBGSdVrUSpmRPVRpTFGRRo1Kmb02qqPjHgqhG5CF3aXkyFyfYVMv4lZqS2BWppmeqR3Ke_NMTJWpwm2H6cdJppu4JwNf5egUe5ZqWb6kB00PtmYn3vTjvFJUrPPaYKOqpsHzXtlAL0dgb5nQtKK9tlgrCenOgijC2ie8W6wQvGhoHxVnZomWVuixh9SCRzzUi0fGq0LUYXzGRoyNiNq3Zw6DxLN-Niukg9nAmOpJpkdEv-TzOZNn1IJxAb5ObUeVfyV4RLQPaB4jk2T_b8qvEN-hV3bR1aJY17fw6k42h60N7Bnf4JRrtTUaCdrDB9k6kIh_XWT4IKQC6D30PC4U77IWQhYGfHIqQdpHsWGLFincli-iHSp-4cWFux-7AU7crFUQpU26F31Fss0d6Zu5zX1ZayFax-ugmDuUAq6G6ZoQHPMdxD7tt8kDwc_Lvml9sL_rzmWseOPc60dJAhdz6cJp303ZVYZFdZRVeP8-GeQuw0elg112h2ItJCi9y3b8DZFDX1EN3yuGfWV3NW0VXGV00OJ-5OVkcKhvjOJ1vvPWt79vNCyB)

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
- bwa
- fcs-adaptor

otb now has a preprint. 
Preprint: https://doi.org/10.32942/X2T897


```
otb is in the public domain in the United States per 17 U.S.C. § 105
```
