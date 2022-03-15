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

[![](https://mermaid.ink/img/pako:eNqlVk1v2zAM_SuCTw7QdsCOOWwY0A277NTdksJQbDrWaouuJDcN0vz3iZJsy4ndAlvRJBT5SD2Rzx-nJMcCknWyV7yt2O_7rWT2r0AJG_q6e_QONLtNp0GxnNe1puWdrh7Z7e0XBvIltR-BSjYgDeOyYKWoQbO8gvwJipUvYTEOr4XcdzVXwhxPXAGTcIh9LEdpuJCgNFPw3AkFxdezLxHDqNTbETR9MgqkBR5kjbx4p-JqoY5E---rLFCBV6FNSO93dLkSXk1Z42HTG0xY4p0MjQtlJ9AQ6fEUcq3KSq4NT53NyH7-5DyrmYQXS0qgTAvQuRI7YBpLc6B-hoheLW6z403YxFp-VnNgG8xs0Nh6qf-9wkcQl_JT_BDfHn6l4ZdaQYIZNBAd0-H_QF0fS6GrdLBC-0LC6Cb4HiQ29sDYQhrZ05R-b0qoRK5NV5Zpb0yhWxmO0dk6JxtyBtsjTcwgs2214F59Qw3XHEJG2aMebUbmXKm1fHB1BfSCG6EFWqUYFmdEOZ-XubFzjBuK07V7HZjQy_jeanskydyaoWQt1rbjUDDsTNv1qr9IdBX9NuGq4JV2RMlY5OlQMzQjf8-SXBFJWkZUxqjLavhxBxmaClTmprTYMn_OfqhXeTPcFjAz3VzuqWvKe_1c6KrQmR8HncdbdJr4TASf1eiQe1FqSl_xveF7x9mbd-1x3Cki66MumZa6quzdWbvEYEfg4BnQPjlDmVW8jTaj3lhPjebYwlKBQ8VN9iRkkWHZt4N8jHwMy0CXLu6Fa_fi7Ha8NJBYuovKjYR7rc9LZY6oiYJ6-LxU_ej_VbCRcD6QbYxcFO9_S3hGyB49M0XPRYumreHN36k8tiTYi30Kc2n0hlasXz1-UK0BZfFvkwofpBQAbYBe5vXkfdUMVUbgk18yVCzK7Uc0KeFHPsn2qPgBFhP3P24A3lzNMYjKLoYn54psVztar8Zj7vLSINZe1_2C3oM0SN3pCRX_HAyY8N44Fp2JDvVmbwozmFgRk4B7WiQ3iSXScFHY99cTobaJFXsD22RtzRIVaLNNtvJskV1bcAPfC2FQJeuS1xpuEt4ZfDjKfHB41L3g9nW4Cd7zXynF6W8)](https://mermaid-js.github.io/mermaid-live-editor/edit/#pako:eNqlVk1v2zAM_SuCTw7QdsCOOWwY0A277NTdksJQbDrWaouuJDcN0vz3iZJsy4ndAlvRJBT5SD2Rzx-nJMcCknWyV7yt2O_7rWT2r0AJG_q6e_QONLtNp0GxnNe1puWdrh7Z7e0XBvIltR-BSjYgDeOyYKWoQbO8gvwJipUvYTEOr4XcdzVXwhxPXAGTcIh9LEdpuJCgNFPw3AkFxdezLxHDqNTbETR9MgqkBR5kjbx4p-JqoY5E---rLFCBV6FNSO93dLkSXk1Z42HTG0xY4p0MjQtlJ9AQ6fEUcq3KSq4NT53NyH7-5DyrmYQXS0qgTAvQuRI7YBpLc6B-hoheLW6z403YxFp-VnNgG8xs0Nh6qf-9wkcQl_JT_BDfHn6l4ZdaQYIZNBAd0-H_QF0fS6GrdLBC-0LC6Cb4HiQ29sDYQhrZ05R-b0qoRK5NV5Zpb0yhWxmO0dk6JxtyBtsjTcwgs2214F59Qw3XHEJG2aMebUbmXKm1fHB1BfSCG6EFWqUYFmdEOZ-XubFzjBuK07V7HZjQy_jeanskydyaoWQt1rbjUDDsTNv1qr9IdBX9NuGq4JV2RMlY5OlQMzQjf8-SXBFJWkZUxqjLavhxBxmaClTmprTYMn_OfqhXeTPcFjAz3VzuqWvKe_1c6KrQmR8HncdbdJr4TASf1eiQe1FqSl_xveF7x9mbd-1x3Cki66MumZa6quzdWbvEYEfg4BnQPjlDmVW8jTaj3lhPjebYwlKBQ8VN9iRkkWHZt4N8jHwMy0CXLu6Fa_fi7Ha8NJBYuovKjYR7rc9LZY6oiYJ6-LxU_ej_VbCRcD6QbYxcFO9_S3hGyB49M0XPRYumreHN36k8tiTYi30Kc2n0hlasXz1-UK0BZfFvkwofpBQAbYBe5vXkfdUMVUbgk18yVCzK7Uc0KeFHPsn2qPgBFhP3P24A3lzNMYjKLoYn54psVztar8Zj7vLSINZe1_2C3oM0SN3pCRX_HAyY8N44Fp2JDvVmbwozmFgRk4B7WiQ3iSXScFHY99cTobaJFXsD22RtzRIVaLNNtvJskV1bcAPfC2FQJeuS1xpuEt4ZfDjKfHB41L3g9nW4Cd7zXynF6W8)

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
