Marincevic 2017
===============

This repo contains the code accompanying the manuscript Transcriptome sequencing in pediatric acute lymphoblastic leukemia identifies fusion genes associated with distinct DNA methylation profiles, by Marincevic et. al.

The script presented here uses a simple strategy for calling fusions genes based on looking for read pairs where the different reads in the pair map to different genes. While sensitive, this method also produces false positives, and is therefore best suited for a targeted approach.

Usage instructions
------------------

**Installation**

 - Install Conda following the instructions here: https://conda.io/docs/install/quick.html
 - Create a new conda environment by using `conda create --name=fusion_detection python=2`
 - Load the newly created environment: `source activate fusion_detection`
 - Install the requirements for the script by using `pip install -r requirements.txt`.

**Running the detect_fusions.py script**

Here is an example of how to run the detection script:

`python detect_fusions.py --bam_file <input bam> --fusion_list_file <fusion list file> --genes_file <genes file> --custom_genes_file <custom genes file> --unknown_genes_output <unknown genes in sample> --output <output>`


The parameters are as follows:
 
 - `<input bam>` a bam file containing paired-end reads from an RNA-seq experiment
 - `<fusion list file>` a file with one fusion per line on the format: `gene1>gene2`
 - `<genes file>` a file containing the coordinate for refseq genes fetched from: http://genome.ucsc.edu/cgi-bin/hgTables?command=start (Make sure to select output format: "all fields from selected table")
 - `<custom genes file>` a file containing genes not included in refseq on the following format:

```
gene    chrom   start   stop
MLLT4   chr6    168225671       168374700
IGH     chr14   105967564       107300892
```
 - `<unknown genes in sample>` is a file where any genes which are found in the `<fusion list file>` but are not found in either the `<genes file>` or the `<custom genes file>` are place
 - `<output>` a tab separated file which contains the support for each fusion. Each column is a fusion and each row is a sample. The values are the number of reads supporting the fusion, normalized by the total number of reads sequenced per sample.


