"""
The goal of this module is to take a BED file containing regions of the genome we wish to exclude.
This module will then examine those regions and build a kmer mutability model based on those regions.
"""


# from cyvcf2 import VCF
# from pyfaidx import Fasta
# import multiprocessing as mp
from eskedit import VCFRegion


def get_bed_regions(bed_path, invert_selection=False, header=False):
    regions = []
    if invert_selection:
        # TODO: implement
        pass
    else:
        with open(bed_path, 'r') as bedfile:
            for line in bedfile.readlines():
                fields = line.split('\t')
                regions.append(VCFRegion(fields[0], fields[1], fields[2]))
    return regions


def train_model(bed_path, vcf_path, fasta_path, kmer_size, nprocs=None, invert_selection=False):
    regions = get_bed_regions(bed_path, invert_selection=False)
    pass
