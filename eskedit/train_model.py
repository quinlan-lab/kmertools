"""
The goal of this module is to take a BED file containing regions of the genome we wish to exclude.
This module will then examine those regions and build a kmer mutability model based on those regions.
"""
from collections import defaultdict, Counter
import time
from cyvcf2 import VCF
from pyfaidx import Fasta
import eskedit as ek
import pandas as pd
import multiprocessing as mp
from eskedit import GRegion, get_autosome_lengths_grch38, get_grch38_chroms, RegionContainer, Variant
import sys

"""
This currently works with 3, 5, and 7mer contexts. The expected input is a tab separated
file with the following fields in this order:
CHROM   POS REF ALT`
    1. CHROM - chromosome should be as reported by grch38
    2. POS - position on chromosome aligned to hg38
    3. REF - reference allele
    4. ALT - alternate allele from VCF
    * Additional fields after this will be ignored
"""


def get_bed_regions(bed_path, invert_selection=True, header=False, clean_bed=False):
    """
    Returns an iterable of GRegions specified by the filepath in bed format.
    :param bed_path:            Path to bed file
    :param invert_selection:    True (default) will return GRegions not in the file
    :param header:              True if file has a header line. False (default)
    :param clean_bed:           False (default) means bed file may have overlapping regions which will be merged. True means each line is added independently of position
    :return:                    Iterable of GRegions
    """
    regions = RegionContainer()
    with open(bed_path, 'r') as bedfile:
        if header:
            bedfile.readline()
        for line in bedfile.readlines():
            fields = line.split('\t')
            if clean_bed:
                regions.add_distinct_region(GRegion(fields[0], fields[1], fields[2]))
            else:
                regions.add_region(GRegion(fields[0], fields[1], fields[2]))
    if invert_selection:
        return regions.get_inverse()
    else:
        return regions.get_regions()


def model_region(master_ref_counts, transitions_list, vcf_path, fasta_path, kmer_size, region):
    start = time.time()
    fasta = Fasta(fasta_path)
    vcf = VCF(vcf_path)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    sequence = fasta.get_seq(region[0], region[1], region[2]).seq.upper()
    region_ref_counts = ek.kmer_search(sequence, kmer_size)  # nprocs=1 due to short region
    r_string = str(region[0]) + ':' + str(region[1]) + '-' + str(region[2])
    transitions = defaultdict(Counter)
    for variant in vcf(r_string):
        if ek.is_quality_singleton(variant):
            new_var = Variant(variant=variant, fields=['vep'])
            # TODO: with less memory overhead variant_positions[new_var.INDEX] = new_var
            # take 7mer around variant. pyfaidx excludes start index and includes end index
            adj_seq = fasta[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
            if adj_seq[kmer_mid_idx] != variant.REF:
                print('WARNING: Reference mismatch\tFasta REF: %s\tVCF REF: %s' % (adj_seq[kmer_mid_idx], variant.REF),
                      file=sys.stderr)
            if ek.complete_sequence(adj_seq):
                transitions[adj_seq.upper()][new_var.ALT[0]] += 1
            # del new_var
    lock.acquire()
    try:
        for k, v in region_ref_counts.items():
            master_ref_counts[k] += v
        for k, v in transitions.items():
            for alt, count in v.items():
                transitions_list[k][alt] += count
    finally:
        lock.release()
    print('Finished region %s in %s' % (str(region), str(time.time() - start)))
    return  # region_ref_counts, transitions


def init_lock(l):
    # source: https://stackoverflow.com/questions/25557686/python-sharing-a-lock-between-processes
    global lock
    lock = l


def train_kmer_model(bed_path, vcf_path, fasta_path, kmer_size, nprocs=None, invert_selection=True, clean_bed=False,
                     header=False):
    """
    Builds the counts tables required for the k-mer model. Returned as 2 dictionaries.
    @param bed_path:            path to bed file
    @param vcf_path:            path to vcf file
    @param fasta_path:          path to reference fasta
    @param kmer_size:
    @param nprocs:              number of processors to use
    @param invert_selection:    True (default) process regions NOT specified by bed file
    @param clean_bed:           False (default) if the bed needs to be merged. True processes regions as is
    @param header:              False (default) if the bed file does not have a header
    @return:
    """
    l = mp.Lock()
    master_ref_counts = Counter()
    master_transitions = defaultdict(Counter)

    regions = get_bed_regions(bed_path, invert_selection=invert_selection, header=header, clean_bed=clean_bed)
    # Bundle arguments to pass to 'model_region' function
    arguments = [(master_ref_counts, master_transitions, vcf_path, fasta_path, kmer_size, region.flist) for region
                 in regions]
    pool = mp.Pool(nprocs, initializer=init_lock, initargs=(l,))
    # Distribute workload
    pool.starmap(model_region, arguments)
    pool.close()
    pool.join()
    # Data structure for counting k-mer occurrences in reference sequences
    # master_ref_counts = Counter()
    # # Contains data for which k-mers mutated to which nucleotide
    # transitions_list = []
    # # merge results from different threads
    # for result in results.get():
    #     for k, v in result[0].items():
    #         master_ref_counts[k] += v
    #
    #     transitions_list.append(result[1])
    # merged_transitions = ek.merge_transitions_ddc(transitions_list, outfile=None)
    return master_ref_counts, master_transitions


def generate_frequency_table(reference_counts, transition_counts, filepath=False, save_file=None):
    if filepath:
        counts = pd.read_csv(reference_counts, index_col=0).sort_index()
        transitions = pd.read_csv(transition_counts, index_col=0).sort_index()
    else:
        counts = pd.DataFrame.from_dict(reference_counts, orient='index').sort_index()
        transitions = pd.DataFrame.from_dict(transition_counts, orient='index').sort_index()
    freq_table = pd.DataFrame()
    if counts.shape[0] != transitions.shape[0]:
        raise ValueError(
            'The reference counts (read %d rows) and transition counts (read %d rows) must have the same number of rows' % (
                counts.shape[0], transitions.shape[0]))
    freq_table['frequency'] = transitions.sum(axis=1)
    freq_table['frequency'] = freq_table['frequency'] / counts.iloc[:, 0]
    if save_file is not None:
        freq_table.to_csv(save_file)
    return freq_table

# if __name__ == "__main__":
#     # test_path = '/Users/simonelongo/too_big_for_icloud/exons_grch38.bed'
#     # test_path = '/Users/simonelongo/too_big_for_icloud/exons_chr1_grch38.bed'
#     test_path = '/Users/simonelongo/too_big_for_icloud/chr1_test.bed'
#     # test_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/input_data/test_data/test.bed'
#     test_reg = get_bed_regions(test_path, invert_selection=True)
#     for tregion in test_reg:
#         print(tregion)
