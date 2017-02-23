# ----------------------------------------------------
# This is small script to extract reads from bam files to find support for fusion genes from RNA-seq data.
# Please note that all files should use 0 based indexing for genomic positions, i.e. the first base of
# chromosome 1 should be denoted as 1:0.
# ----------------------------------------------------

import os
import argparse
import logging as log
import sys

import pysam

root = log.getLogger()
root.setLevel(log.DEBUG)

ch = log.StreamHandler(sys.stdout)
ch.setLevel(log.DEBUG)
formatter = log.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
root.addHandler(ch)

CHROMOSOMES_OF_INTEREST = map(str, range(1,23)) + ["X", "Y"]

class Gene(object):
    def __init__(self, chromosome, strand, start, stop, name):
        if chromosome.startswith("chr"):
            self.chromosome = chromosome.replace("chr", "")
        else:
            self.chromosome = chromosome
        self.strand = strand
        self.start = int(start)
        self.stop = int(stop)
        self.name = name

    def __str__(self):
        return "{}:{}-{}\t{}".format(self.chromosome, self.start, self.stop, self.name)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __hash__(self):
        return hash("".join(map(str, self.__dict__.values())))


class Fusion(object):
    def __init__(self, name, five_prime_name, three_prime_name):
        self.name = name
        self.five_prime_name = five_prime_name
        self.three_prime_name = three_prime_name

        self.five_prime_genes = set()
        self.three_prime_genes = set()

        self.supporting_reads = 0
        self.normalized_supporting_reads = None
        self.matching_reads = []

    def add_five_prime_genes(self, five_prime_genes):
        if self.five_prime_genes:
            raise Exception("Trying to add to existing fusion (5')")
        self.five_prime_genes = five_prime_genes

    def add_three_prime_genes(self, three_prime_genes):
        if self.three_prime_genes:
            raise Exception("Trying to add to existing fusion (3')")
        self.three_prime_genes = three_prime_genes

    def __str__(self):
        return "{}>{}".format(self.five_prime_name, self.three_prime_name)


class UnknownFusions(object):
    def __init__(self):
        self.unknown_fusions = set()

    def add(self, fusion):
        self.unknown_fusions.add(fusion)

    def write_unknown_fusions_to_file(self, file_name):
        with file(file_name, "w") as f:
            for fusion in self.unknown_fusions:
                f.write(str(fusion) + "\n")


def read_fusions_from_file(fusion_file):
    with file(fusion_file, "r") as f:
        for line in f:
            if line.strip():
                five_prime, three_prime = line.split(">")
                yield Fusion(line.strip(), five_prime.strip(), three_prime.strip())


def generate_custom_genes(custom_gene_list):
    with file(custom_gene_list, "r") as f:
        for line in f:
            if not line.startswith("gene"):
                split_string = line.split("\t")
                new_gene = Gene(chromosome=split_string[1].strip(),
                                strand='N/A',
                                start=split_string[2].strip(),
                                stop=split_string[3].strip(),
                                name=split_string[0].strip())
                # Only keep 1 to 22 + X and Y
                if new_gene.chromosome in CHROMOSOMES_OF_INTEREST:
                    yield new_gene


def generate_genes(gene_list):
    with file(gene_list, "r") as f:
        for line in f:
            if not line.startswith("#"):
                split_string = line.split("\t")
                new_gene = Gene(chromosome=split_string[2].strip(),
                                strand=split_string[3].strip(),
                                start=split_string[4].strip(),
                                stop=split_string[5].strip(),
                                name=split_string[12].strip())
                # Only keep 1 to 22 + X and Y
                if new_gene.chromosome in CHROMOSOMES_OF_INTEREST:
                    yield new_gene


def gene_names_to_genes_list(genes):
    gene_names_to_genes = {}
    for gene in genes:
        gene_list_for_name = gene_names_to_genes.get(gene.name)

        if gene_list_for_name:
            gene_list_for_name.append(gene)
        else:
            gene_names_to_genes[gene.name] = [gene]
    return gene_names_to_genes


def add_genes_for_fusion(fusion, gene_lookup, unknown_fusions):
    try:
        five_prime = gene_lookup.get(fusion.five_prime_name)
        three_prime = gene_lookup.get(fusion.three_prime_name)

        if not five_prime or not three_prime:
            unknown_fusions.add(fusion)

        if five_prime:
            fusion.add_five_prime_genes(set(five_prime))
        if three_prime:
            fusion.add_three_prime_genes(set(three_prime))

    except Exception as e:
        log.error("five prime partner: {}, three prime partner: {}".format(fusion.five_prime_name,
                                                                           fusion.three_prime_name))
        raise e


def generate_fusions(fusions_list_file, names2genes, unknown_fusions):
    for fusion in read_fusions_from_file(fusions_list_file):
        add_genes_for_fusion(fusion, names2genes, unknown_fusions)
        yield fusion


def is_fused(read, three_prime_fusions):
    for tpf in three_prime_fusions:
        if str(read.next_reference_name) == str(tpf.chromosome) and tpf.start <= read.next_reference_start < tpf.stop:
            return True
    return False


def count_unique_reads(reads):
    s = set()
    for read in reads:
        s.add((read.reference_name, read.reference_start, read.next_reference_name, read.next_reference_start))
    return len(s)


def parse_supporting_reads_from_alignments(bam_file, fusions):
    with pysam.AlignmentFile(bam_file, "rb") as alignments:

        log.info("Starting to count total number of alignments in bam file")
        total_nbr_of_reads = alignments.count()
        log.info("Found a total of {} alignments in file.".format(total_nbr_of_reads))

        total_nbr_of_fusions = len(fusions)
        fusions_counter = 0

        for fusion in fusions:

            fusions_counter += 1
            if fusions_counter % 100 == 0:
                log.info("Parsed reads for {0:.0f}% of the input fusions".format(
                    100*(float(fusions_counter)/total_nbr_of_fusions)))

            matching_reads = []
            for fp in fusion.five_prime_genes:
                try:
                    for read in alignments.fetch(fp.chromosome, fp.start, fp.stop, multiple_iterators=True):
                        if not read.mate_is_unmapped and read.is_read1 and is_fused(read, fusion.three_prime_genes):
                            matching_reads.append(read)
                except Exception as e:
                    log.error("Trouble parsing read: {}".format(read.tostring(alignments)))
                    raise e
            fusion.matching_reads = matching_reads
            fusion.supporting_reads = count_unique_reads(matching_reads)
            fusion.normalized_supporting_reads = float(fusion.supporting_reads) / float(total_nbr_of_reads)

        log.info("Finished parsing, parsed a total of {} fusions".format(fusions_counter))

# ---------------------------------------------------
#  Actual script starts here
# ---------------------------------------------------


if __name__ == "__main__":

    # Parse input arguments

    parser = argparse.ArgumentParser(description="Script used to look for paired-end read supporting fusions in a "
                                                 "bam file. "
                                                 "Example of execution: \n"
                                                 "python detect_fusions.py --bam_file 00_252.bam "
                                                 "--fusion_list_file FINAL_COMPLETE_looping_file_list.txt "
                                                 "--genes_file refseq_hg19.txt "
                                                 "--unknown_genes_output 00_252_unknown_fusions.txt "
                                                 "--output 00_252_fusions_support.tsv")

    parser.add_argument("--bam_file", help="The bam file to for look fusions in.", required=True)
    parser.add_argument("--fusion_list_file", help="A list of fusions to look for. One fusion per row, with format: "
                                                   "[five prime partner]>[three prime partner].", required=True)
    parser.add_argument("--genes_file", help="File with gene annotation with format from ucsc genome browser, with "
                                             "output format 'all fields from selected table'.", required=True)
    parser.add_argument("--custom_genes_file", help="Add extra annotations for genes not found in refseq",
                        required=True)
    parser.add_argument("--unknown_genes_output", help="The output file where all fusions which "
                                                       "involve genes which are not found in the 'genes_file' are "
                                                       "placed.", required=True)
    parser.add_argument("--output", help="File where output is placed. The format used is one row per sample, with "
                                         "support for each fusion in the columns, measured as supporting reads per "
                                         "total reads sequenced. (Tab separated)", required=True)

    args = parser.parse_args()

    # Input data
    bam_file = args.bam_file
    fusions_list_file = args.fusion_list_file
    genes_file = args.genes_file
    custom_genes_file = args.custom_genes_file

    log.info("Start processing: {}".format(bam_file))

    # Output files
    output_file = args.output
    unknown_genes_output = args.unknown_genes_output

    # Unknown fusions will be gathered here
    unknown_fusions = UnknownFusions()

    # Get all genes from gene list file
    genes = generate_genes(genes_file)
    custom_genes = generate_custom_genes(custom_genes_file)
    all_genes = list(genes) + list(custom_genes)
    names2genes = gene_names_to_genes_list(all_genes)

    fusions_as_list = list(generate_fusions(fusions_list_file, names2genes, unknown_fusions))

    log.info("Parsing support for fusions from bam file")

    # Annotate the fusions list with how many supporting reads we have
    # for each particular fusion
    parse_supporting_reads_from_alignments(bam_file, fusions_as_list)

    log.info("Write output files, {} and {}".format(output_file, unknown_genes_output))

    # Write output files
    def write_row(f, row):
        f.write(row + "\n")

    header = ["name"] + map(lambda x: x.name, fusions_as_list)
    with file(output_file, "w") as f:
        write_row(f, "\t".join(header))
        sample_name = os.path.basename(bam_file).replace(".bam", "")
        norm_support_reads = map(lambda x: str(x.normalized_supporting_reads), fusions_as_list)
        write_row(f, "\t".join([sample_name] + norm_support_reads))

    unknown_fusions.write_unknown_fusions_to_file(unknown_genes_output)

    log.info("Finished, have a nice day...")

