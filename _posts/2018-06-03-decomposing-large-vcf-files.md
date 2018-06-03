---
layout: post
title: Decomposing large VCF files with Python
author: Adrian Baez-Ortega
date: 2018/06/03
---


When looking at genetic variation in a set of sequencing data files, the process normally starts with a VCF file. [VCF (Variant Calling Format)](https://samtools.github.io/hts-specs/VCFv4.2.pdf) is a standard file format for genetic variation calls, and is used by most, if not all, variant calling software tools. A VCF file is a text file that looks something like this:

    ##fileformat=VCFv4.x
    ##...
    ##...file header
    ##...
    #CHROM  POS    ID  REF  ALT  QUAL  FILTER  INFO       FORMAT              Sample1        ...
    1       26530  .   T    C    2948  PASS    BRF=0;...  GT:GL:GOF:GQ:NR:NV  0/0:...:209:0  ...
    ...
    ...

The first 8 columns in each line (named `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`) contain information about the variant in question. I tend to call this the 'variant metadata', although it's probably not the most correct name. But I'm too old to change now. 

After these 8 columns follows the information about the genotype of the variant in each of the samples; these are numerical values separated by colons (`:`). The `FORMAT` column indicates what each of these colon-delimited values mean, which is something that changes depending on the variant calling software. In the example above, the `FORMAT` column tells us that the values in the following columns represent, respectively, the genotype (`GT`), the genotype likelihood (`GL`), the goodness of fit (`GOF`), the genotype quality (`GQ`), the number of reads at the variant locus (`NR`) and the number of reads supporting the variant at the variant locus (`NV`). This is normally explained somewhere in the header of the VCF file.

If you want to analyse these variants in, say, R, you would need to read this VCF and decompose it so that you can access the colon-delimited values in the genotype information. If your VCF is small (say, less than 1 GB) then you can read it directly into R with the `read.table` function. However, my VCFs are rarely small (my personal record is 50 GB for a single file), which means that directly reading the VCF from R is not an option. Even worse, most of the space in the VCF is taken up by the colon-separated values in the sample genotype columns, and we are not normally interested in all of these values (I only use two of them). 

So I thought, wouldn't it be nice to decompose the VCF _outside_ R (in a faster language) and import only the data that I'm interested in?

To do this, __I wrote a small Python script ([`decompose_vcf.py`](https://github.com/baezortega/misc/blob/master/decompose_vcf.py))__ which processes a VCF file (with `.vcf` extension), or a gzipped VCF file (with `.vcf.gz` extension) and extracts the variant metadata (first 8 columns), and any other values that we need from the genotype data (e.g. `GT`, `GL`, etc.), into separate text files. This works regardless of which variant caller was used to make the VCF, because the format of the metadata columns if fixed, whereas the format of the genotype data is specified in the `FORMAT` column. (I suspect this was cleverly made on purpose.)

I wrote the script for Python 2.7, but it shouldn't be too hard to adapt for Python 3. It requires the Python library `gzip` in order to open gzipped VCF files.

To know how to use the script, just call it with no arguments:

``` sh
./decompose_vcf.py
```

    ## decompose_vcf.py
    ## This script processes a VCF, or compressed (gzipped) VCF file, extracting
    ## the variant metadata fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO),
    ## and those values specified by the user (by default, all) from the genotype
    ## information of each variant in each sample. Each of these values is output
    ## as a table to a separate output text file.
    ##
    ## Input: Path to input VCF (.vcf) or gzipped VCF (.vcf.gz) file.
    ##        [Optional] List of value names of interest (as shown in the FORMAT field)
    ##
    ## Usage: decompose_vcf.py /path/to/file.vcf[.gz] [GT GL NR NV ...]

So, if you have a large VCF named, say, `toomanyvariants.vcf`, and you are interested only in the fields `GT`, `NR` and `NV`, then you can extract those via:

``` sh
./decompose_vcf.py toomanyvariants.vcf GT NR NV
```

    ## Input file:
    ##   toomanyvariants.vcf
    ## Extracting variant metadata and values:
    ##   GT, NR, NV
    ## Output files:
    ##   toomanyvariants_GT.txt
    ##   toomanyvariants_NR.txt
    ##   toomanyvariants_NV.txt
    ##   toomanyvariants_Meta.txt
    ## Processing VCF...
    ## Done!

If we wanted to extract all the fields in the genotype data, we would just specify no field names and provide only the VCF file name, i.e. `./decompose_vcf.py toomanyvariants.vcf`.

Let's compare the size of the input and output files.

``` sh
ls -lh toomanyvariants* | awk '{print $5"\t"$9}'
```

    ## 6.5G	toomanyvariants.vcf
    ## 802M	toomanyvariants_GT.txt
    ##  44M	toomanyvariants_Meta.txt
    ## 628M	toomanyvariants_NR.txt
    ## 429M	toomanyvariants_NV.txt

The four output files should be quite easy to read from R using `read.table` (with argument `header=TRUE`). And the best thing is that the data are already decomposed into separate tables, so we can forget about the VCF format from here on.


---
