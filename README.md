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

[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggVERcbiAgICBkb25lW2RvbmUuXVxuICAgIG90Ylt1c2VyIGNhbGxzIG90Yi5zaF0gLS0-IGVudihlbnZpb3JubWVudCBhbmQgZmlsZXMgY2hlY2tlZClcbiAgICBlbnYgLS0-IHNpbmd1bGFyaXR5e2FyZSBuZXcgc2luZ3VsYXJpdHkgY29udGFpbmVycyByZXF1aXJlZD99XG4gICAgc2luZ3VsYXJpdHkgLS0-IHx5ZXN8eWVzX3NpbmcoZG93bmxvYWQgbmV3IHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMpXG4gICAgc2luZ3VsYXJpdHkgLS0-IHxub3xub19zaW5nKHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMgZXhpc3QpXG4gICAgeWVzX3NpbmcgLS0-IG5leHRmbG93W25leHRmbG93IGlzIHJ1bl1cbiAgICBub19zaW5nIC0tPiBuZXh0Zmxvd1xuICAgIG5leHRmbG93IC0tPiBjaGVja19mYXN0YShjaGVjayBmYXN0cS9mYXN0YSlcbiAgICBuZXh0ZmxvdyAtLT4gdmVyc2lvbihkZXNjcmliZSBzb2Z0d2FyZSB2ZXJzaW9ucylcbiAgICBuZXh0ZmxvdyAtLT4gY2hlY2tfYmFtKGNoZWNrIGJhbSBmaWxlcylcbiAgICBuZXh0ZmxvdyAtLT4gYmFtX2ZpbHRlcnMoZmlsdGVyIGJhbSBmaWxlcylcbiAgICBiYW1fZmlsdGVycyAtLT4gSGlGaUFTTShjYWxsSGlGaUFTTSBpcyBjYWxsZWQpXG4gICAgY2hlY2tfZmFzdGEgLS0-IGplbGx5ZmlzaChqZWxseWZpc2ggaXMgcnVuKVxuICAgIGplbGx5ZmlzaCAtLT4gZ2Vub21lc2NvcGUoZ2Vub21lc2NvcGUgaXMgcnVuKVxuICAgIEhpRmlBU00gLS0-IGhpY3N0dWZmKGhpY3N0dWZmIGlzIHJ1bilcbiAgICBcblxuICAgIGJ1c2Nve2lzIGJ1c2NvIGdvaW5nIHRvIGJlIHJ1bj99XG4gICAgaGljc3R1ZmYgLS0-IGJ1c2NvXG4gICAgYnVzY28gLS0-IHx5ZXN8cnVuX2J1c2NvKHJ1biBidXNjbylcbiAgICBidXNjbyAtLT4gfG5vfG5vX3J1bl9idXNjbyhkbyBub3QgcnVuIGJ1c2NvKVxuXG5cbiAgICBidXNjbzJ7aXMgYnVzY28gZ29pbmcgdG8gYmUgcnVuP31cbiAgICBidXNjbzIgLS0-IHxub3xkb25lXG4gICAgYnVzY28yIC0tPiB8eWVzfHJ1bl9idXNjb19hZ2FpbihydW4gYnVzY28gYWdhaW4pXG4gICAgcnVuX2J1c2NvX2FnYWluIC0tPiBkb25lXG5cblxuICAgIGlzX3BvbGlzaHtpcyBwb2xpc2hpbmcgZ29pbmcgdG8gYmUgZG9uZT99XG4gICAgaGljc3R1ZmYgLS0-IGlzX3BvbGlzaFxuICAgIGlzX3BvbGlzaCAtLT4gfHllc3xydW5fcmFndGFnKHJ1biByYWd0YWcucHkgcG9saXNoaW5nKVxuICAgIHJ1bl9yYWd0YWcgLS0-IHJ1bl9zaGhxdWlzKHJ1biBzaGhxdWlzKVxuICAgIHJ1bl9zaGhxdWlzIC0tPiBydW5fcmFndGFnX29uX2hhcChydW4gcmFndGFnIG9uIGhhcGxvdHlwZXMpXG4gICAgcnVuX3NoaHF1aXMgLS0-IHdoYXRfa2luZF9vZl9wb2xpc2h7d2hhdCBraW5kIG9mIHBvbGlzaCBpcyBnb2luZyB0byBiZSBydW4_fVxuICAgIGlzX3BvbGlzaCAtLT4gfG5vfGRvbmVcblxuICAgIHdoYXRfa2luZF9vZl9wb2xpc2ggLS0-IHxzaW1wbGV8YnVzY28yXG4gICAgZmluZF92YXJpYW50c1tmaW5kIHZhcmlhbnRzXVxuICAgIHdoYXRfa2luZF9vZl9wb2xpc2ggLS0-IHxtZXJmaW58ZmluZF92YXJpYW50c1xuICAgIHdoYXRfa2luZF9vZl9wb2xpc2ggLS0-IHxkZWVwdmFyaWFudHxmaW5kX3ZhcmlhbnRzXG5cblxuICAgIG1lcmZpbl9vcl9kZWVwe21lcmZpbiBvciBkZWVwdmFyaWFudD99XG4gICAgZmluZF92YXJpYW50cyAtLT4gbWVyZmluX29yX2RlZXBcbiAgICBnZW5vbWVzY29wZSAtLT4gfG1lcmZpbnxtZXJmaW4ocnVuIG1lcmZpbilcbiAgICBtZXJmaW5fb3JfZGVlcCAtLT4gbWVyZmluXG4gICAgbWVyZmluX29yX2RlZXAgLS0-IHxkZWVwdmFyaWFudHxkZWVwdmFyaWFudChydW4gZGVlcHZhcmlhbnQpXG5cblxuICAgIGJjZnRvb2xzKHJ1biBiY2Z0b29scyBjb25zZW5zdXMpXG4gICAgbWVyZmluIC0tPiBiY2Z0b29sc1xuICAgIGRlZXB2YXJpYW50IC0tPiBiY2Z0b29sc1xuICAgIGJjZnRvb2xzIC0tPiByYWd0YWcyKHJ1biByYWd0YWcgYWdhaW4pXG4gICAgc2hocXVpczIgLS0-IHJ1bl9yYWd0YWdfb25faGFwXG4gICAgcmFndGFnMiAtLT4gc2hocXVpczIocnVuIHNoaHF1aXMgYWdhaW4pXG4gICAgc2hocXVpczIgLS0-IGJ1c2NvMlxuXG5cbiIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0In0sInVwZGF0ZUVkaXRvciI6ZmFsc2UsImF1dG9TeW5jIjpmYWxzZSwidXBkYXRlRGlhZ3JhbSI6ZmFsc2V9)](https://mermaid-js.github.io/mermaid-live-editor/edit/#eyJjb2RlIjoiZ3JhcGggVERcbiAgICBkb25lW2RvbmUuXVxuICAgIG90Ylt1c2VyIGNhbGxzIG90Yi5zaF0gLS0-IGVudihlbnZpb3JubWVudCBhbmQgZmlsZXMgY2hlY2tlZClcbiAgICBlbnYgLS0-IHNpbmd1bGFyaXR5e2FyZSBuZXcgc2luZ3VsYXJpdHkgY29udGFpbmVycyByZXF1aXJlZD99XG4gICAgc2luZ3VsYXJpdHkgLS0-IHx5ZXN8eWVzX3NpbmcoZG93bmxvYWQgbmV3IHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMpXG4gICAgc2luZ3VsYXJpdHkgLS0-IHxub3xub19zaW5nKHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMgZXhpc3QpXG4gICAgeWVzX3NpbmcgLS0-IG5leHRmbG93W25leHRmbG93IGlzIHJ1bl1cbiAgICBub19zaW5nIC0tPiBuZXh0Zmxvd1xuICAgIG5leHRmbG93IC0tPiBjaGVja19mYXN0YShjaGVjayBmYXN0cS9mYXN0YSlcbiAgICBuZXh0ZmxvdyAtLT4gdmVyc2lvbihkZXNjcmliZSBzb2Z0d2FyZSB2ZXJzaW9ucylcbiAgICBuZXh0ZmxvdyAtLT4gY2hlY2tfYmFtKGNoZWNrIGJhbSBmaWxlcylcbiAgICBuZXh0ZmxvdyAtLT4gYmFtX2ZpbHRlcnMoZmlsdGVyIGJhbSBmaWxlcylcbiAgICBiYW1fZmlsdGVycyAtLT4gSGlGaUFTTShjYWxsSGlGaUFTTSBpcyBjYWxsZWQpXG4gICAgY2hlY2tfZmFzdGEgLS0-IGplbGx5ZmlzaChqZWxseWZpc2ggaXMgcnVuKVxuICAgIGplbGx5ZmlzaCAtLT4gZ2Vub21lc2NvcGUoZ2Vub21lc2NvcGUgaXMgcnVuKVxuICAgIEhpRmlBU00gLS0-IGhpY3N0dWZmKGhpY3N0dWZmIGlzIHJ1bilcbiAgICBcblxuICAgIGJ1c2Nve2lzIGJ1c2NvIGdvaW5nIHRvIGJlIHJ1bj99XG4gICAgaGljc3R1ZmYgLS0-IGJ1c2NvXG4gICAgYnVzY28gLS0-IHx5ZXN8cnVuX2J1c2NvKHJ1biBidXNjbylcbiAgICBidXNjbyAtLT4gfG5vfG5vX3J1bl9idXNjbyhkbyBub3QgcnVuIGJ1c2NvKVxuXG5cbiAgICBidXNjbzJ7aXMgYnVzY28gZ29pbmcgdG8gYmUgcnVuP31cbiAgICBidXNjbzIgLS0-IHxub3xkb25lXG4gICAgYnVzY28yIC0tPiB8eWVzfHJ1bl9idXNjb19hZ2FpbihydW4gYnVzY28gYWdhaW4pXG4gICAgcnVuX2J1c2NvX2FnYWluIC0tPiBkb25lXG5cblxuICAgIGlzX3BvbGlzaHtpcyBwb2xpc2hpbmcgZ29pbmcgdG8gYmUgZG9uZT99XG4gICAgaGljc3R1ZmYgLS0-IGlzX3BvbGlzaFxuICAgIGlzX3BvbGlzaCAtLT4gfHllc3xydW5fcmFndGFnKHJ1biByYWd0YWcucHkgcG9saXNoaW5nKVxuICAgIHJ1bl9yYWd0YWcgLS0-IHJ1bl9zaGhxdWlzKHJ1biBzaGhxdWlzKVxuICAgIHJ1bl9zaGhxdWlzIC0tPiBydW5fcmFndGFnX29uX2hhcChydW4gcmFndGFnIG9uIGhhcGxvdHlwZXMpXG4gICAgcnVuX3NoaHF1aXMgLS0-IHdoYXRfa2luZF9vZl9wb2xpc2h7d2hhdCBraW5kIG9mIHBvbGlzaCBpcyBnb2luZyB0byBiZSBydW4_fVxuICAgIGlzX3BvbGlzaCAtLT4gfG5vfGRvbmVcblxuICAgIHdoYXRfa2luZF9vZl9wb2xpc2ggLS0-IHxzaW1wbGV8YnVzY28yXG4gICAgZmluZF92YXJpYW50c1tmaW5kIHZhcmlhbnRzXVxuICAgIHdoYXRfa2luZF9vZl9wb2xpc2ggLS0-IHxtZXJmaW58ZmluZF92YXJpYW50c1xuICAgIHdoYXRfa2luZF9vZl9wb2xpc2ggLS0-IHxkZWVwdmFyaWFudHxmaW5kX3ZhcmlhbnRzXG5cblxuICAgIG1lcmZpbl9vcl9kZWVwe21lcmZpbiBvciBkZWVwdmFyaWFudD99XG4gICAgZmluZF92YXJpYW50cyAtLT4gbWVyZmluX29yX2RlZXBcbiAgICBnZW5vbWVzY29wZSAtLT4gfG1lcmZpbnxtZXJmaW4ocnVuIG1lcmZpbilcbiAgICBtZXJmaW5fb3JfZGVlcCAtLT4gbWVyZmluXG4gICAgbWVyZmluX29yX2RlZXAgLS0-IHxkZWVwdmFyaWFudHxkZWVwdmFyaWFudChydW4gZGVlcHZhcmlhbnQpXG5cblxuICAgIGJjZnRvb2xzKHJ1biBiY2Z0b29scyBjb25zZW5zdXMpXG4gICAgbWVyZmluIC0tPiBiY2Z0b29sc1xuICAgIGRlZXB2YXJpYW50IC0tPiBiY2Z0b29sc1xuICAgIGJjZnRvb2xzIC0tPiByYWd0YWcyKHJ1biByYWd0YWcgYWdhaW4pXG4gICAgc2hocXVpczIgLS0-IHJ1bl9yYWd0YWdfb25faGFwXG4gICAgcmFndGFnMiAtLT4gc2hocXVpczIocnVuIHNoaHF1aXMgYWdhaW4pXG4gICAgc2hocXVpczIgLS0-IGJ1c2NvMlxuXG5cbiIsIm1lcm1haWQiOiJ7XG4gIFwidGhlbWVcIjogXCJkZWZhdWx0XCJcbn0iLCJ1cGRhdGVFZGl0b3IiOmZhbHNlLCJhdXRvU3luYyI6ZmFsc2UsInVwZGF0ZURpYWdyYW0iOmZhbHNlfQ)

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

```
otb is in the public domain in the United States per 17 U.S.C. § 105
```
