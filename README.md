# Ariadne: Barcoded Read Deconvolution

De novo assemblies are critical for capturing the genetic composition of complex samples. Linked-read sequencing techniques such as 10x Genomics’ Linked-Reads, UST’s TELL-Seq, Loop Genomics’ LoopSeq, and BGI’s Long Fragment Read (LFR) combines 30 barcoding with standard short-read sequencing to expand the range of linkage resolution from hundreds to tens of thousands of base-pairs. The application of linked-read sequencing to genome assembly has demonstrated that barcoding-based technologies balance the tradeoffs between long-range linkage, per-base coverage, and costs.

Linked-reads come with their own challenges, chief among them the association of multiple long fragments with the same 3' barcode. The lack of a unique correspondence between a long fragment and a barcode, in conjunction with low sequencing depth, confounds the assignment of linkage between short-reads.

Here, we introduce Ariadne, a novel assembly graph-based algorithm, that can be used to deconvolve a large metagenomic linked-read dataset. As demonstrated substantial increases in the largest alignments, contigs, and genomic contiguity, *de novo* assemblies from deconvolved reads represent advancements in terms of completeness and accuracy. Ariadne is intuitive, computationally efficient, and scalable to other large-scale linked-read problems, such as human genome phasing. 

The Ariadne manuscript will be publicly available soon. 

## Installation

Currently, Ariadne is implemented as a module of an older version of the cloudSPAdes *de novo* assembly program (version 3.12.1). The following libraries need to be pre-installed:

* g++ (version 5.3.1 or higher)
* cmake (version 2.8.12 or higher)
* zlib
* libbz2

From source: 
```
git clone <url>   
cd ariadne/
./spades_compile.sh
```
The installation directory can be set by `PREFIX=<destination_dir>` in the compile step. 

In the future, Ariadne will be repackaged as a standalone program, along with a scaffold generator such that deconvolved reads can be directly used to generate *de novo* assemblies without a second cloudSPAdes run.

## Deconvolving Reads

Use the following command to run barcode deconvolution. The option `gemcode` must be set for the main cloudSPAdes modules to recognize barcoded reads. `<fastq_name>` should be separate or interleaved fastq file where reads have a `BX` tag designating the barcode (this is the default output of [longranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)). `<max_search_dist>` is the user-specified parameter for the maximum search distance, which should be smaller than the average length of a genomic fragment. `<min_cloud_size>` is the user-specified parameter for the minimum cloud size- in terms of number of reads- for which the deconvolution process will be run. By default, these parameters are set to 5 kbp and 6 reads respectively. Future versions of Ariadne will allow for the user to specify the barcode identifier. BayesHammer error correction is turned off because the introduced tags interfere with barcode recognition. BayesHammer may be run separately from the assembly procedure to generate error-corrected reads, as long as the barcode format described in **Input** below is followed as input for the actual cloudSPAdes command. 
```
spades.py [--meta] --only-assembler --gemcode1-1 <fastq_name>.R1.fastq --gemcode1-2 <fastq_name>.R2.fastq --search-distance <max_search_dist> --size_cutoff <min_cloud_size> -t <num_threads> -m <mem_in_gb> -o /path/to/output_dir
```

For more SPAdes options, refer to the [Spades manual](http://cab.spbu.ru/files/release3.13.1/manual.html) or the command-line options.
```
spades.py
```

### Input

The input for Ariadne are the paired-end FastQ files provided to cloudSPAdes. All barcoded reads must be suffixed with a `-1`, as the output module replaces this number with its own inferred cloud number. For example, `@D00547:847:HYHNTBCXX:1:1101:10000:10626 BX:Z:CCTTCCCTCCTTCAAT-1` will be appropriately processed, while `@D00547:847:HYHNTBCXX:1:1101:10000:10626 BX:Z:CCTTCCCTCCTTCAAT` will not. If your FastQ files do not carry the `-1` suffix, as with newer linked-read technologies such as TELL-Seq and LoopSeq, they can be modified as such:

```
sed '1~4s/$/-1/' <fastq_name>.R1.fastq > <fastq_suffixed>.R1.fastq
```

### Output

At the end of the deconvolution procedure, Ariadne outputs paired FastQ files with enhanced barcode assignments in the directory `/path/to/output_dir/K55/<max_search_dist>.RX.fastq`. The original barcode has been augmented with the enhanced grouping number. In this example, the first set of paired reads `@D00547:847:HYHNTBCXX:1:1101:10000:10626` have been assigned to the 13th group out of all reads with the barcode `BX:Z:CCTTCCCTCCTTCAAT`. This fastq can be directly provided as input to a *de novo* assembler or a read mapper optimized for linked-reads.

```
$ head /path/to/output/K55/<max_search_dist>.R1.fastq
@D00547:847:HYHNTBCXX:1:1101:10000:10626 BX:Z:CCTTCCCTCCTTCAAT-13
ATGCTGGGGTTTCCGCTGCAATTCTTTGTCCGGTTCTTTAAGAACCACGGCTTGCTGTCGATCAGCAACCGCCCACAGTGGTGCGTGATCGAAGGCGGCTCCAGCAGCTACATCGAGCCGCTGACCC
+
IIGGIIIGIIIIIGIIIGGGGGGIIIGIIIIIIGGGGIGGGIIIIIGIGGGGGGGIIGGGGGGGIIIIIIIIIIGGGIGIIIIGGIIGIIIGIIIAGGG.GGIIGGGGGGGGGGIIIGGIIIIGGGI
@D00547:847:HYHNTBCXX:1:1101:10000:11336 BX:Z:CAGGTATCAGCGTAAG-1
TTCAGCAAGCGCAGCTTGATGGCATCCCGCAATTGATACAGGTTCAGGCCGGTGGTGTTGATCTTGAGGTCGGCCAGATCGATGATCGGTCCCAGCAGCGAGGTCTCGTCCTCGATGGCTTCGGCCA
+
A<.AGGGGGGGAGAAGGGIGGAAGGGGGGGGGGGGIGGGGGGIGGGIG.GGGAGGGGGGGG.GGGGGAGGGGG.AGA.AGG<GG<AAGGGGGGGGGAGGG.GGGIGGAG<<<AGAAGAAG<GGGGG.
@D00547:847:HYHNTBCXX:1:1101:10000:11336 BX:Z:CAGGTATCAGCGTAAG-1
GTCCTGGATGAACAAAACGGGCAGTAATCATGCGTTTGATCATCGTCAGCGGCCGCTCCGGCTCGGGTAAAAGCACCGCCCTCAACGTCCTTGAAGACAACGGCTTTTATTGCATCGACAACCTTCCCGCCGGTTTGCTGCCGGAGTTGG
+
G<..GGGGGGGGA<GGA.<GGGGGAGGG.GGAGAGAGG.AGG<<<..<...AAA<AA<G....<GAGGIGIGAGAG.....<AGG..<AAAA<..<...<.<...<<AGG...GGG<GAAAAGG..<A<AGGAG.<<.G...<7<<7GGG
```

### Performance

A full exploration of performance metrics will be available with the manuscript. For example, for metagenomics-sized dataset such as the full MOCK5 dataset (97 million reads approx., see below), deconvolution at a maximum search distance of 5 kbp takes 6 hours 20 minutes, 62 GB RAM (in addition to the memory requirements of cloudSPAdes), and 20 CPUs. 

## Datasets

The MOCK5 and MOCK20 10x datasets used in the paper may be downloaded from [AWS](https://s3.us-east-2.amazonaws.com/readclouds/cloudspades_data.tar.gz). The MOCK5 LoopSeq and MOCK20 TELL-Seq datasets can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA728470. 

## Credits and Citations

Find the Ariadne preprint [here](https://www.biorxiv.org/content/10.1101/2021.05.09.443255v1). If you've found Ariadne useful, please cite it as:

Ariadne: Barcoded Linked-Read Deconvolution Using de Bruijn Graphs
Lauren Mak, Dmitry Meleshko, David C. Danko, Waris N. Barakzai, Natan Belchikov, Iman Hajirasouliha
bioRxiv 2021.05.09.443255; doi: https://doi.org/10.1101/2021.05.09.443255

This algorithm was developed by Waris Barakzai and myself, and tested by myself with help from Dmitrii Meleshko, David Danko, Natan Belchikov and Iman Hajirasouliha.
