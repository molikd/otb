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

[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggVERcbiAgICBkb25lW2RvbmUuXVxuICAgIG90Ylt1c2VyIGNhbGxzIG90Yi5zaF0gLS0-IGVudihlbnZpb3JubWVudCBhbmQgZmlsZXMgY2hlY2tlZClcbiAgICBlbnYgLS0-IHNpbmd1bGFyaXR5e2FyZSBuZXcgc2luZ3VsYXJpdHkgY29udGFpbmVycyByZXF1aXJlZD99XG4gICAgc2luZ3VsYXJpdHkgLS0-IHx5ZXN8eWVzX3NpbmcoZG93bmxvYWQgbmV3IHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMpXG4gICAgc2luZ3VsYXJpdHkgLS0-IHxub3xub19zaW5nKHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMgZXhpc3QpXG4gICAgeWVzX3NpbmcgLS0-IG5leHRmbG93W25leHRmbG93IGlzIHJ1bl1cbiAgICBub19zaW5nIC0tPiBuZXh0Zmxvd1xuICAgIG5leHRmbG93IC0tPiBjaGVja19mYXN0YShmYXN0cS9mYXN0YSBjaGVja2VkKVxuICAgIG5leHRmbG93IC0tPiB2ZXJzaW9uKHNvZnR3YXJlIHZlcnNpb25zIGRlc2NyaWJlZClcbiAgICBuZXh0ZmxvdyAtLT4gY2hlY2tfYmFtKGJhbSBmaWxlcyBjaGVja2VkKVxuICAgIG5leHRmbG93IC0tPiBiYW1fZmlsdGVycyhiYW0gZmlsZXMgZmlsdGVyZWQpXG4gICAgYmFtX2ZpbHRlcnMgLS0-IEhpRmlBU00oSGlGaUFTTSBpcyBjYWxsZWQpXG4gICAgY2hlY2tfZmFzdGEgLS0-IGplbGx5ZmlzaChqZWxseWZpc2ggaXMgcnVuKVxuICAgIGplbGx5ZmlzaCAtLT4gZ2Vub21lc2NvcGUoZ2Vub21lc2NvcGUgaXMgcnVuKVxuICAgIEhpRmlBU00gLS0-IGhpY3N0dWZmKGhpY3N0dWZmIGlzIHJ1bilcbiAgICBcblxuICAgIGJ1c2Nve2lzIGJ1c2NvIGdvaW5nIHRvIGJlIHJ1bj99XG4gICAgaGljc3R1ZmYgLS0-IGJ1c2NvXG4gICAgYnVzY28gLS0-IHx5ZXN8cnVuX2J1c2NvKHJ1biBidXNjbylcbiAgICBidXNjbyAtLT4gfG5vfG5vX3J1bl9idXNjbyhkbyBub3QgcnVuIGJ1c2NvKVxuXG5cbiAgICBidXNjbzJ7aXMgYnVzY28gZ29pbmcgdG8gYmUgcnVuP31cbiAgICBidXNjbzIgLS0-IHxub3xkb25lXG4gICAgYnVzY28yIC0tPiB8eWVzfHJ1bl9idXNjb19hZ2FpbihydW4gYnVzY28gYWdhaW4pXG4gICAgcnVuX2J1c2NvX2FnYWluIC0tPiBkb25lXG5cblxuICAgIGlzX3BvbGlzaHtpcyBwb2xpc2hpbmcgZ29pbmcgdG8gYmUgZG9uZX1cbiAgICBoaWNzdHVmZiAtLT4gaXNfcG9saXNoXG4gICAgaXNfcG9saXNoIC0tPiB8eWVzfHJ1bl9yYWd0YWcocnVuIHJhZ3RhZy5weSBwb2xpc2hpbmcpXG4gICAgcnVuX3JhZ3RhZyAtLT4gcnVuX3NoaHF1aXMocnVuIHNoaHF1aXMpXG4gICAgcnVuX3NoaHF1aXMgLS0-IHJ1bl9yYWd0YWdfb25faGFwKHJ1biByYWd0YWcgb24gaGFwbG90eXBlcylcbiAgICBydW5fc2hocXVpcyAtLT4gd2hhdF9raW5kX29mX3BvbGlzaHt3aGF0IGtpbmQgb2YgcG9saXNoIGlzIGdvaW5nIHRvIGJlIHJ1bj99XG4gICAgaXNfcG9saXNoIC0tPiB8bm98ZG9uZVxuXG4gICAgd2hhdF9raW5kX29mX3BvbGlzaCAtLT4gfHNpbXBsZXxidXNjbzJcbiAgICBmaW5kX3ZhcmlhbnRzW2ZpbmQgdmFyaWFudHNdXG4gICAgd2hhdF9raW5kX29mX3BvbGlzaCAtLT4gfG1lcmZpbnxmaW5kX3ZhcmlhbnRzXG4gICAgd2hhdF9raW5kX29mX3BvbGlzaCAtLT4gfGRlZXB2YXJpYW50fGZpbmRfdmFyaWFudHNcblxuXG4gICAgbWVyZmluX29yX2RlZXB7bWVyZmluIG9yIGRlZXB2YXJpYW50P31cbiAgICBmaW5kX3ZhcmlhbnRzIC0tPiBtZXJmaW5fb3JfZGVlcFxuICAgIGdlbm9tZXNjb3BlIC0tPiB8bWVyZmlufG1lcmZpbihydW4gbWVyZmluKVxuICAgIG1lcmZpbl9vcl9kZWVwIC0tPiBtZXJmaW5cbiAgICBtZXJmaW5fb3JfZGVlcCAtLT4gfGRlZXB2YXJpYW50fGRlZXB2YXJpYW50KHJ1biBkZWVwdmFyaWFudClcblxuXG4gICAgYmNmdG9vbHMocnVuIGJjZnRvb2xzIGNvbnNlbnN1cylcbiAgICBtZXJmaW4gLS0-IGJjZnRvb2xzXG4gICAgZGVlcHZhcmlhbnQgLS0-IGJjZnRvb2xzXG4gICAgYmNmdG9vbHMgLS0-IHJhZ3RhZzIocnVuIHJhZ3RhZyBhZ2FpbilcbiAgICBzaGhxdWlzMiAtLT4gcnVuX3JhZ3RhZ19vbl9oYXBcbiAgICByYWd0YWcyIC0tPiBzaGhxdWlzMihydW4gc2hocXVpcyBhZ2FpbilcbiAgICBzaGhxdWlzMiAtLT4gYnVzY28yXG5cblxuIiwibWVybWFpZCI6eyJ0aGVtZSI6ImRlZmF1bHQifSwidXBkYXRlRWRpdG9yIjpmYWxzZSwiYXV0b1N5bmMiOmZhbHNlLCJ1cGRhdGVEaWFncmFtIjpmYWxzZX0)](https://mermaid-js.github.io/mermaid-live-editor/edit/#eyJjb2RlIjoiZ3JhcGggVERcbiAgICBkb25lW2RvbmUuXVxuICAgIG90Ylt1c2VyIGNhbGxzIG90Yi5zaF0gLS0-IGVudihlbnZpb3JubWVudCBhbmQgZmlsZXMgY2hlY2tlZClcbiAgICBlbnYgLS0-IHNpbmd1bGFyaXR5e2FyZSBuZXcgc2luZ3VsYXJpdHkgY29udGFpbmVycyByZXF1aXJlZD99XG4gICAgc2luZ3VsYXJpdHkgLS0-IHx5ZXN8eWVzX3NpbmcoZG93bmxvYWQgbmV3IHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMpXG4gICAgc2luZ3VsYXJpdHkgLS0-IHxub3xub19zaW5nKHNpbmd1bGFyaXR5IGNvbnRhaW5lcnMgZXhpc3QpXG4gICAgeWVzX3NpbmcgLS0-IG5leHRmbG93W25leHRmbG93IGlzIHJ1bl1cbiAgICBub19zaW5nIC0tPiBuZXh0Zmxvd1xuICAgIG5leHRmbG93IC0tPiBjaGVja19mYXN0YShmYXN0cS9mYXN0YSBjaGVja2VkKVxuICAgIG5leHRmbG93IC0tPiB2ZXJzaW9uKHNvZnR3YXJlIHZlcnNpb25zIGRlc2NyaWJlZClcbiAgICBuZXh0ZmxvdyAtLT4gY2hlY2tfYmFtKGJhbSBmaWxlcyBjaGVja2VkKVxuICAgIG5leHRmbG93IC0tPiBiYW1fZmlsdGVycyhiYW0gZmlsZXMgZmlsdGVyZWQpXG4gICAgYmFtX2ZpbHRlcnMgLS0-IEhpRmlBU00oSGlGaUFTTSBpcyBjYWxsZWQpXG4gICAgY2hlY2tfZmFzdGEgLS0-IGplbGx5ZmlzaChqZWxseWZpc2ggaXMgcnVuKVxuICAgIGplbGx5ZmlzaCAtLT4gZ2Vub21lc2NvcGUoZ2Vub21lc2NvcGUgaXMgcnVuKVxuICAgIEhpRmlBU00gLS0-IGhpY3N0dWZmKGhpY3N0dWZmIGlzIHJ1bilcbiAgICBcblxuICAgIGJ1c2Nve2lzIGJ1c2NvIGdvaW5nIHRvIGJlIHJ1bj99XG4gICAgaGljc3R1ZmYgLS0-IGJ1c2NvXG4gICAgYnVzY28gLS0-IHx5ZXN8cnVuX2J1c2NvKHJ1biBidXNjbylcbiAgICBidXNjbyAtLT4gfG5vfG5vX3J1bl9idXNjbyhkbyBub3QgcnVuIGJ1c2NvKVxuXG5cbiAgICBidXNjbzJ7aXMgYnVzY28gZ29pbmcgdG8gYmUgcnVuP31cbiAgICBidXNjbzIgLS0-IHxub3xkb25lXG4gICAgYnVzY28yIC0tPiB8eWVzfHJ1bl9idXNjb19hZ2FpbihydW4gYnVzY28gYWdhaW4pXG4gICAgcnVuX2J1c2NvX2FnYWluIC0tPiBkb25lXG5cblxuICAgIGlzX3BvbGlzaHtpcyBwb2xpc2hpbmcgZ29pbmcgdG8gYmUgZG9uZX1cbiAgICBoaWNzdHVmZiAtLT4gaXNfcG9saXNoXG4gICAgaXNfcG9saXNoIC0tPiB8eWVzfHJ1bl9yYWd0YWcocnVuIHJhZ3RhZy5weSBwb2xpc2hpbmcpXG4gICAgcnVuX3JhZ3RhZyAtLT4gcnVuX3NoaHF1aXMocnVuIHNoaHF1aXMpXG4gICAgcnVuX3NoaHF1aXMgLS0-IHJ1bl9yYWd0YWdfb25faGFwKHJ1biByYWd0YWcgb24gaGFwbG90eXBlcylcbiAgICBydW5fc2hocXVpcyAtLT4gd2hhdF9raW5kX29mX3BvbGlzaHt3aGF0IGtpbmQgb2YgcG9saXNoIGlzIGdvaW5nIHRvIGJlIHJ1bj99XG4gICAgaXNfcG9saXNoIC0tPiB8bm98ZG9uZVxuXG4gICAgd2hhdF9raW5kX29mX3BvbGlzaCAtLT4gfHNpbXBsZXxidXNjbzJcbiAgICBmaW5kX3ZhcmlhbnRzW2ZpbmQgdmFyaWFudHNdXG4gICAgd2hhdF9raW5kX29mX3BvbGlzaCAtLT4gfG1lcmZpbnxmaW5kX3ZhcmlhbnRzXG4gICAgd2hhdF9raW5kX29mX3BvbGlzaCAtLT4gfGRlZXB2YXJpYW50fGZpbmRfdmFyaWFudHNcblxuXG4gICAgbWVyZmluX29yX2RlZXB7bWVyZmluIG9yIGRlZXB2YXJpYW50P31cbiAgICBmaW5kX3ZhcmlhbnRzIC0tPiBtZXJmaW5fb3JfZGVlcFxuICAgIGdlbm9tZXNjb3BlIC0tPiB8bWVyZmlufG1lcmZpbihydW4gbWVyZmluKVxuICAgIG1lcmZpbl9vcl9kZWVwIC0tPiBtZXJmaW5cbiAgICBtZXJmaW5fb3JfZGVlcCAtLT4gfGRlZXB2YXJpYW50fGRlZXB2YXJpYW50KHJ1biBkZWVwdmFyaWFudClcblxuXG4gICAgYmNmdG9vbHMocnVuIGJjZnRvb2xzIGNvbnNlbnN1cylcbiAgICBtZXJmaW4gLS0-IGJjZnRvb2xzXG4gICAgZGVlcHZhcmlhbnQgLS0-IGJjZnRvb2xzXG4gICAgYmNmdG9vbHMgLS0-IHJhZ3RhZzIocnVuIHJhZ3RhZyBhZ2FpbilcbiAgICBzaGhxdWlzMiAtLT4gcnVuX3JhZ3RhZ19vbl9oYXBcbiAgICByYWd0YWcyIC0tPiBzaGhxdWlzMihydW4gc2hocXVpcyBhZ2FpbilcbiAgICBzaGhxdWlzMiAtLT4gYnVzY28yXG5cblxuIiwibWVybWFpZCI6IntcbiAgXCJ0aGVtZVwiOiBcImRlZmF1bHRcIlxufSIsInVwZGF0ZUVkaXRvciI6ZmFsc2UsImF1dG9TeW5jIjpmYWxzZSwidXBkYXRlRGlhZ3JhbSI6ZmFsc2V9)

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
