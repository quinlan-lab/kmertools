"""
The goal of this module is to take a BED file containing regions of the genome we wish to exclude.
This module will then examine those regions and build a kmer mutability model based on those regions.
"""
from collections import defaultdict, Counter
import time
from cyvcf2 import VCF
from pyfaidx import Fasta, FetchError
import eskedit as ek
import pandas as pd
import multiprocessing as mp
from eskedit import GRegion, get_autosome_lengths_grch38, get_grch38_chroms, RegionContainer, Variant, DataContainer
import sys
import array

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


def get_bed_regions(bed_path, invert_selection=True, header=False, clean_bed=False, strand_col=None,
                    bed_names_col=None):
    """
    Returns an iterable of GRegions specified by the filepath in bed format.
    :param bed_path:            Path to bed file
    :param invert_selection:    True (default) will return GRegions not in the file
    :param header:              True if file has a header line. False (default)
    :param clean_bed:           False (default) means bed file may have overlapping regions which will be merged. True means each line is added independently of position
    :param strand_col:          Zero-based column index containing strand information
    :param bed_names_col:       Zero-based column index containing name information
    :return:                    Iterable of GRegions

    """
    additional_fields = defaultdict(int)
    if strand_col is not None:
        try:
            strand_col = int(strand_col)
            additional_fields['strand'] = strand_col
        except ValueError:
            strand_col = None
    if bed_names_col is not None:
        try:
            bed_names_col = int(bed_names_col)
            additional_fields['name'] = bed_names_col
        except ValueError:
            bed_names_col = None
    regions = RegionContainer()
    with open(bed_path, 'r') as bedfile:
        kwargs = defaultdict(str)
        if header:
            bedfile.readline()
        for line in bedfile.readlines():
            fields = line.split('\t')
            # add keyword arguments to pass to GRegion constructor
            for k, v in additional_fields.items():
                kwargs[k] = fields[v]

            if clean_bed:
                regions.add_distinct_region(GRegion(*[fields[0], fields[1], fields[2]], **kwargs))
            else:
                regions.add_region(GRegion(*[fields[0], fields[1], fields[2]], **kwargs))
    if invert_selection:
        return regions.get_inverse()
    else:
        return regions.get_regions()


def model_region_singletons(data_container, vcf_path, fasta_path, kmer_size, region):
    start = time.time()
    fasta = Fasta(fasta_path)
    vcf = VCF(vcf_path)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    try:
        if region.strand is not None:
            if ek.is_dash(region.strand):
                sequence = fasta.get_seq(region.chrom, region.start-kmer_mid_idx, region.stop+kmer_mid_idx).complement.seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start-kmer_mid_idx, region.stop+kmer_mid_idx).seq.upper()
        else:
            sequence = fasta.get_seq(region.chrom, region.start-kmer_mid_idx, region.stop+kmer_mid_idx).seq.upper()
    except (KeyError, FetchError):
        print('Region %s not found in fasta, continuing...' % str(region), file=sys.stderr, flush=True)
        return
    region_ref_counts = ek.kmer_search(sequence, kmer_size)  # nprocs=1 due to short region
    r_string = str(region.chrom) + ':' + str(region.start) + '-' + str(region.stop)
    transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    # Define indices for nucleotides
    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    idx_nuc = list('ACGT')
    for variant in vcf(r_string):
        if ek.is_singleton_snv(variant):
            new_var = Variant(variant=variant, fields=['vep'])
            # take 7mer around variant. pyfaidx excludes start index and includes end index
            adj_seq = fasta[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
            if str(adj_seq[kmer_mid_idx]).upper() != str(variant.REF).upper():
                print('WARNING: Reference mismatch\tFasta REF: %s\tVCF REF: %s' % (adj_seq[kmer_mid_idx], variant.REF),
                      file=sys.stderr, flush=True)
            if ek.complete_sequence(adj_seq):
                transitions[adj_seq.upper()][nuc_idx[new_var.ALT[0]]] += 1
    temp = data_container.get()
    temp.add_count(region_ref_counts)
    temp.add_transition(transitions)
    data_container.set(temp)
    print('Finished region %s in %s' % (str(region), str(time.time() - start)), flush=True)
    return


def model_region_nonsingletons(data_container, vcf_path, fasta_path, kmer_size, region, AC_cutoff):
    if AC_cutoff is not None:
        try:
            AC_cutoff = int(AC_cutoff)
        except ValueError:
            AC_cutoff = None
            print('AC cutoff must be a positive integer. Ignoring user value and using SNVs with any AC.', file=sys.stderr, flush=True)
    try:
        kmer_size = int(kmer_size)
        if kmer_size < 1: raise ValueError
    except ValueError:
        print('kmer_size must be a positive integer. Please check arguments.', file=sys.stderr, flush=True)
        exit(1)
    start = time.time()
    fasta = Fasta(fasta_path)
    vcf = VCF(vcf_path)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    try:
        if region.strand is not None:
            if ek.is_dash(region.strand):
                sequence = fasta.get_seq(region.chrom, region.start-kmer_mid_idx, region.stop+kmer_mid_idx).complement.seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start-kmer_mid_idx, region.stop+kmer_mid_idx).seq.upper()
        else:
            sequence = fasta.get_seq(region.chrom, region.start-kmer_mid_idx, region.stop+kmer_mid_idx).seq.upper()
    except (KeyError, FetchError):
        print('Region %s not found in fasta, continuing...' % str(region), file=sys.stderr, flush=True)
        return
    region_ref_counts = ek.kmer_search(sequence, kmer_size)  # nprocs=1 due to short region
    r_string = str(region.chrom) + ':' + str(region.start) + '-' + str(region.stop)
    ac_transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    an_transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    # Define indices for nucleotides
    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    idx_nuc = list('ACGT')
    for variant in vcf(r_string):
        if ek.is_quality_snv(variant, AC_cutoff=AC_cutoff):
            new_var = Variant(variant=variant)
            adj_seq = fasta[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
            if str(adj_seq[kmer_mid_idx]).upper() != str(variant.REF).upper():
                print('WARNING: Reference mismatch\tFasta REF: %s\tVCF REF: %s' % (adj_seq[kmer_mid_idx], variant.REF),
                      file=sys.stderr, flush=True)
            if ek.complete_sequence(adj_seq):
                ac_transitions[adj_seq.upper()][nuc_idx[new_var.ALT[0]]] += new_var.AC
                an_transitions[adj_seq.upper()][nuc_idx[new_var.ALT[0]]] += new_var.AN
        # if ek.is_singleton_snv(variant):
        #     new_var = Variant(variant=variant, fields=['vep'])
        #     # take 7mer around variant. pyfaidx excludes start index and includes end index
        #     adj_seq = fasta[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
        #     if str(adj_seq[kmer_mid_idx]).upper() != str(variant.REF).upper():
        #         print('WARNING: Reference mismatch\tFasta REF: %s\tVCF REF: %s' % (adj_seq[kmer_mid_idx], variant.REF), file=sys.stderr, flush=True)
        #     if ek.complete_sequence(adj_seq):
        #         transitions[adj_seq.upper()][nuc_idx[new_var.ALT[0]]] += 1
    temp = data_container.get()
    temp.add_count(region_ref_counts)
    temp.add_transition(ac_transitions)
    temp.add_transition2(an_transitions)
    data_container.set(temp)
    print('Finished region %s in %s' % (str(region), str(time.time() - start)), flush=True)
    return


def train_kmer_model(bed_path, vcf_path, fasta_path, kmer_size, nprocs=1, invert_selection=True,
                     clean_bed=False, singletons=False, nonsingletons=False,
                     header=False, strand_col=None, bed_names_col=None, AC_cutoff=None):
    """
    Builds the counts tables required for the k-mer model. Returned as 2 dictionaries.
    @param AC_cutoff:           Specify to filter out variants above a given AC (AC > cutoff will be filtered)
                                ** works only if keyword 'nonsingletons=True' **
    @param nonsingletons:       Set true to train model on all SNVs
    @param singletons:          Set true to train model based on singleton variants
    @param strand_col:          zero-based column index of strand information from bed file
    @param bed_names_col:       zero-based column index of name information from bed file
    @param bed_path:            path to bed file
    @param vcf_path:            path to vcf file
    @param fasta_path:          path to reference fasta
    @param kmer_size:           NOTE unpredictable behavior may occur if even numbers are used here.
    @param nprocs:              number of processors to use
    @param invert_selection:    True (default) process regions NOT specified by bed file
    @param clean_bed:           False (default) if the bed needs to be merged. True processes regions as is
    @param header:              False (default) if the bed file does not have a header
    @return:
    """
    try:
        nprocs = int(nprocs)
        if strand_col is not None:
            strand_col = int(strand_col)
        if bed_names_col is not None:
            bed_names_col = int(bed_names_col)
    except ValueError:
        print('ERROR: nprocs and column indices must be integers', file=sys.stderr, flush=True)
        exit(1)
    manager = mp.Manager()
    # set up so master data count stays in shared memory
    dc = manager.Value(DataContainer, DataContainer())
    regions = get_bed_regions(bed_path, invert_selection=invert_selection, header=header, clean_bed=clean_bed,
                              strand_col=strand_col, bed_names_col=bed_names_col)
    # Bundle arguments to pass to 'model_region' function

    pool = mp.Pool(nprocs)
    # Distribute workload
    results = defaultdict(tuple)
    if singletons:
        arguments = [(dc, vcf_path, fasta_path, kmer_size, region) for region in regions]
        pool.starmap_async(model_region_singletons, arguments)
        pool.close()
        pool.join()
        results['singletons'] = dc.value.get()
    if nonsingletons:
        arguments = [(dc, vcf_path, fasta_path, kmer_size, region, AC_cutoff) for region in regions]
        pool.starmap_async(model_region_nonsingletons, arguments)
        pool.close()
        pool.join()
        results['all_variants'] = dc.value.get()
    return results  # master_ref_counts, transitions_list


def generate_frequency_table(reference_counts, transition_counts, filepath=False, save_file=None):
    if filepath:
        counts = pd.read_csv(reference_counts, index_col=0).sort_index()
        transitions = pd.read_csv(transition_counts, index_col=0).sort_index()
    else:
        counts = pd.DataFrame.from_dict(reference_counts, orient='index').sort_index()
        transitions = pd.DataFrame.from_dict(transition_counts, orient='index').sort_index()
    freq_table = pd.DataFrame()
    counts.fillna(0, inplace=True)
    transitions.fillna(0, inplace=True)
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
