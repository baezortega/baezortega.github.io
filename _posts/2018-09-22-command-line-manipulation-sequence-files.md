---
layout: post
title: Command-line manipulation of sequence files
author: Adrian Baez-Ortega
date: 2018/09/22
---


Here I present some command-line approaches, hoarded over the years, for manipulating and converting between different popular sequence data files formats, namely [BAM/SAM](https://en.wikipedia.org/wiki/SAM_(file_format\)), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) and [FASTA](https://en.wikipedia.org/wiki/FASTA_format). While they are by no means the only methods available for these tasks, they are perhaps the most simple.

### 1. Converting BAM to SAM

BAM is just a binary compressed version of the SAM format, and the easier way to convert between them is via `samtools view` command from the [SAMtools](http://www.htslib.org/) software suite. The `-h` option allows including the BAM file header in the output SAM, which is often necessary.

``` sh
samtools view -h INPUT.bam > OUTPUT.sam
```

### 2. Converting SAM to BAM

Here, the `-bS` options indicate that the input is SAM and the output is BAM.

``` sh
samtools view -bS INPUT.sam > OUTPUT.bam
```

### 3. Converting SAM/BAM to FASTQ

Converting from SAM to FASTQ (both of which are text-based formats) may be as simple as selecting specific columns of the SAM file and arranging them into lines of the FASTQ file. The initial `grep` command discards the SAM file header, which cannot go into the FASTQ file.

``` sh
grep -v ^@ INPUT.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > OUTPUT.fastq
```

Obviously, the SAM file itself can be produced from a BAM file using `samtools view` as shown above. There is, however, a problem with this `grep|awk` approach: each record in a SAM/BAM file represent a read *alignment* to a reference sequence, and it is possible for a given sequence read to have more than one alignment. In contrast, each record in a FASTQ should represent a *unique read*, which implies that we should only consider one alignment (the primary alignment) for each read in the SAM/BAM file. 

The easiest way (to my knowledge) of generating a FASTQ file from a BAM file, such that reads with multiple alignments are considered only once, is using the `bamtofastq` command from the [biobambam2](https://www.sanger.ac.uk/science/tools/biobambam) software suite, which automatically ignores secondary alignments.

By default, `bamtofastq` generates two FASTQ files, containing the forward and reverse reads in each read pair, respectively. The output file names are indicated via the `F` and `F2` arguments. Orphan forward and reverse reads lacking a paired mate are also output to two files, indicated by `O` and `O2`. For simplicity, the same file name can be used for `F` and `O`, and for `F2` and `O2` (meaning that we don't care whether the reads have a mate, but only whether they are forward or reverse).

``` sh
bamtofastq filename=INPUT.bam F=OUTPUT_1.fastq O=OUTPUT_1.fastq \
F2=OUTPUT_2.fastq O2=OUTPUT_2.fastq
```

Additionally, the `inputformat` argument allows using SAM input files, while `gz=1` can be used to generate compressed FASTQ files (`fq.gz` extension).

``` sh
bamtofastq inputformat=sam gz=1 filename=INPUT.sam F=OUTPUT_1.fq.gz O=OUTPUT_1.fq.gz \
F2=OUTPUT_2.fq.gz O2=OUTPUT_2.fq.gz
```

Finally, if we simply want all the reads to be in a single FASTQ file, we can just omit all the output name arguments and redirect the output to a file.

``` sh
bamtofastq filename=INPUT.bam > OUTPUT.fastq
```

For more info, run `bamtofastq -h`.

### 4. Converting FASTQ to FASTA

This involves simply selecting the first two lines of every four-line record in the FASTQ file, and replacing the `@` at the beginning of the read identifier with `>`. Note that, for this to work, each sequence must be on a single line (without line breaks), that is, every record in the FASTQ must be composed of exactly four lines.

``` sh
awk '{ if ('NR%4==1' || 'NR%4==2'){ gsub("@",">",$1); print } }' INPUT.fastq > OUTPUT.fasta
```

### 5. Sorting FASTQ files by read identifier

To sort a FASTQ file with four-line records (i.e. with each sequence on a single line) by read ID, we can transform it into a four-column table, sort by the value of the first column, and convert back to FASTQ format, all within a single line. The `-S` argument in the `sort` command allows us to define the maximum amount of memory to use (`4G` means a maximum of 4 GB).

``` sh
cat INPUT.fastq | paste - - - - | sort -k1,1 -S 4G | tr '\t' '\n' > OUTPUT.fastq
```

### 6. Merging FASTQ files with identical read names

When we have a pair of FASTQ files containing, respectively, the forward and reverse reads in each read pair, sometimes we find that the two reads in a same pair have identical identifiers, such that if we merge both FASTQ files into one we end up with duplicated read IDs. This can be solved by adding a suffix (in this case, `_1` or `_2`) to the IDs in each file while merging them (by concatenation).

``` sh
awk '{ if ('NR%4==1'){ $1=$1"_1" } print }' INPUT_1.fastq > OUTPUT.fastq
awk '{ if ('NR%4==1'){ $1=$1"_2" } print }' INPUT_2.fastq >> OUTPUT.fastq
```

### 7. Merging BAM files with different read groups

The `samtools merge` command from the [SAMtools](http://www.htslib.org/) suite allows us to easily merge two or more BAM files as follows.

``` sh
samtools merge OUTPUT.bam INPUT_1.bam INPUT_2.bam [INPUT_3.bam ...]
```

However, this command preserves only the header of the first input file (`INPUT_1.bam`) and discards the headers of all subsequent files. This means that, if the reads in the different BAM files belong to different read groups, we will lose their read group information, as this is contained in the `@RG` records within the BAM header. Therefore, the `@RG` records of the subsequent input files first need to be appended to a copy of the first file's header, before providing `samtools merge` with the resulting modified header (via the `-h` argument). After merging, the output BAM should be indexed with `samtools index`.

``` sh
samtools view -H INPUT_1.bam > header.txt
samtools view -H INPUT_2.bam | grep "@RG" >> header.txt
## Repeat the line above for INPUT_3.bam, etc.
samtools merge -h header.txt OUTPUT.bam INPUT_1.bam INPUT_2.bam [INPUT_3.bam ...]
samtools index OUTPUT.bam
```
