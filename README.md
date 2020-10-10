# midseq.py

---

```
Get midseq form fastq reads

usage: midseq.py [-h] [-f F] [-t T] [-o OUTPUT] [-of OUTPUTFASTA]
                 [Fastq files [Fastq files ...]]

Get the internal sequence between 5' and 3' boundary sequence from fastq file.
E.g., get string "ATLAS" from "GGGATLASAAA" with -f "GGG" and -t "AAA".

positional arguments:
  Fastq files           HTS reads in fastq format (not in .gz format). Support
                        muti-fastq-files.

optional arguments:
  -h, --help            show this help message and exit
  -f F                  The 5' sequence. Default is: CCGGCGACGTTGGGTCAACT
  -t T                  The 3' sequence. Default is: TGTCCTCTTCCTCTTTAGCG
  -o OUTPUT, --output OUTPUT
                        The output file name.
  -of OUTPUTFASTA, --outputfasta OUTPUTFASTA
                        The output fasta name. The default is not to output
                        this file unless assigned its file name.
```
