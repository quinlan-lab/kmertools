import os
import random
from collections import defaultdict
from cyvcf2 import VCF, Writer
import numpy as np
from eskedit.constants import get_autosome_names_grch38, get_autosome_lengths_grch38
from eskedit.kclass import GRegion


def get_split_chrom_vcf(vcf_path, chrom, nbins):
    vcf = VCF(vcf_path)
    nbins = int(nbins)
    names = get_autosome_names_grch38()
    if chrom not in names:
        print('Chromosome key, %s, not found in VCF, please check for typos.' % str(chrom))
    chrom_len = 0
    for kv in zip(vcf.seqnames, vcf.seqlens):
        if kv[0] == chrom:
            chrom_len = kv[1]
    chunk_size = chrom_len // nbins
    regions = []
    start = 1
    stop = chunk_size
    for i in range(nbins):
        if i < nbins - 1:
            regions.append(GRegion(chrom, start, stop))
            start = stop
            stop += chunk_size
        else:
            regions.append(GRegion(chrom, start, chrom_len))
    return regions


def get_split_vcf_regions(vcf_path, nprocs):
    vcf = VCF(vcf_path)
    nprocs = int(nprocs)
    num_entries = np.sum(vcf.seqlens)
    chunk_size = int(num_entries / nprocs) + 1
    names = get_autosome_names_grch38()
    num_chunk = 0
    regions = []
    gen_pos = 0
    current_chromosome = 0
    chrom_pos = 0
    while num_chunk < nprocs and current_chromosome < len(names):
        current_chunk = 0
        region = []
        while current_chunk < chunk_size:
            if current_chromosome >= len(names):
                current_chunk = chunk_size
                continue
            remaining_chunk = chunk_size - current_chunk
            if remaining_chunk <= (vcf.seqlens[current_chromosome] - chrom_pos):
                new_region = GRegion(vcf.seqnames[current_chromosome], chrom_pos, chrom_pos + remaining_chunk)
                region.append(new_region)
                chrom_pos += remaining_chunk
                current_chunk += new_region.size()
                gen_pos += new_region.size()
                continue
            else:  # remaining chunk can fit remainder of chromosome and then some
                new_region = GRegion(vcf.seqnames[current_chromosome], chrom_pos, vcf.seqlens[current_chromosome])
                region.append(new_region)
                current_chunk += new_region.size()
                gen_pos += new_region.size()
                chrom_pos = 0
                current_chromosome += 1
        regions.append(region)
    # print(regions)
    return regions


def split_seq(sequence, nprocs, overlap=None):
    """
    :param sequence: string sequence to split
    :param nprocs: number of CPUs to use
    :param overlap: for use when counting k-mers. If considering 3-mers, overlap=3
    :return: a list of 'nprocs' roughly equal chunks to pass as an argument to python multiprocessing
    """
    chunk_size = int(len(sequence) / nprocs) + 1
    args = []
    start = 0
    end = chunk_size
    for proc in range(nprocs):
        if overlap is not None:
            args.append(sequence[start:(end + overlap - 1)])
        else:
            args.append(sequence[start:end])
        start = end
        end += chunk_size
        if end > len(sequence):
            end = len(sequence)
    return args


def get_bed_windows(window, overlap=0, outfile=None, bins=False):
    chrom_info = get_autosome_lengths_grch38()
    regions = []
    for chrom, length in chrom_info.items():
        if bins:
            window_size = length // window + 1
        else:
            window_size = window
        if window_size > length:
            print("Window size is larger than chromosome. Please reduce window size")
            exit(1)
        if overlap >= window_size:
            print("Overlap value (%d) is too large for window size (%d). Using overlap=0" % (overlap, window_size))
            overlap = 0
        current_position = 1
        while current_position < length:
            regions.append((chrom, max(1, current_position - overlap), min(length, current_position + window_size - 1)))
            current_position += window_size
    if outfile is not None:
        with open(outfile, 'w') as bedfile:
            for region in regions:
                bedfile.write('%s\t%s\t%s\n' % (str(region[0]), str(region[1]), str(region[2])))
    return regions


def append_variants_to_vcf(filename, chrom, start, stop, outfile='samp_build38.vcf'):
    return "tabix " + filename + " " + str(chrom) + ":" + str(
        start) + "-" + str(
        stop) + " >> " + outfile


def generate_sample_vcf(vcf_path, outfile='samp_build38.vcf'):
    """Takes a large VCF file and takes random samples from each chromosome to make a smaller VCF for testing"""
    vcf = VCF(vcf_path)
    write = Writer(outfile, vcf)
    write.write_header()
    key_values = zip(vcf.seqnames, vcf.seqlens)
    chrom_keys = defaultdict(int)
    chroms = get_autosome_names_grch38()
    for kv in list(key_values):
        if kv[0] in chroms:
            chrom_keys[kv[0]] = kv[1]
    for chrom_num, chrom_len in chrom_keys.items():
        begin = random.randint(100000, chrom_len - 100000)
        os.system(append_variants_to_vcf(vcf_path, chrom_num, begin, begin + 10_000, outfile=outfile))
    write.close()
