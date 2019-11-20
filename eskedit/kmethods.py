import random
import multiprocessing as mp
import time
from collections import Counter, defaultdict
import numpy as np
from cyvcf2 import VCF
import re
import pandas as pd
from pyfaidx import Fasta
from eskedit.kclass import Variant
from eskedit.ksplit import split_seq, get_split_vcf_regions


def ran_seq(seq_len):
    """
    Driver method for gen_random_sequence
    :param seq_len:
    :return: random lsequence of length seq_len
    """
    sequence = ""
    for i in range(seq_len):
        sequence += random.choice('ACTG')
    return sequence


def gen_random_sequence(length):
    """
    :param length: (integer) number of nucleotides in a random sequence
    :return: a random DNA sequence of length 'length'
    """
    if length > 5_000_000:
        pool = mp.Pool()
        nprocs = mp.cpu_count()
        chunk_size = int(length / nprocs) + 1
        args = []
        tot = chunk_size
        while tot <= length:
            args.append(chunk_size)
            tot += chunk_size
        new_tot = np.array(args).sum()
        if new_tot < length:
            args.append(length - new_tot)
        results = [funccall.get() for funccall in [pool.map_async(ran_seq, args)]]
        pool.close()
        random_seq = ""
        for seq in results[0]:
            random_seq += seq
    else:
        random_seq = ran_seq(length)
    return random_seq


def complement(c):
    """
    :param c: Nucleotide to get complement of
    :return: character representing the complement of 'c'
    """
    base_pairs = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N'
    }
    try:
        return base_pairs[c.upper()]
    except KeyError:
        raise ValueError(c + " is not a valid nucleotide.")


def get_complementary_sequence(sequence):
    """
    Returns a string of nucleotides complementary to the input string
    All letters in input sequence must A, C, T, or G, otherwise will raise a ValueError
    """
    comp_seq = []
    for c in sequence[::-1]:  # take the reverse complement
        comp_seq.append(complement(c))
    return "".join(comp_seq)


def is_quality_variant(var_to_test):
    """
    high quality variants will have FILTER == None
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def is_quality_nonsingleton(var_to_test):
    """
    high quality variants will have FILTER == None
    Additionally, variants shoud NOT be singletons ('AC' != 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') != 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def is_quality_singleton(var_to_test):
    """
    high quality variants will have FILTER == None
    Additionally, variants shoud be singletons ('AC' == 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') == 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def complete_sequence(seq):
    """
    Seq is normalized by capitalizing all of the characters
    :param seq: sequence to test
    :return: True if sequence only contains A, C, T, or G, False if contains N or any other character
    """
    allowed = set('ACTG')
    return set(seq.upper()).issubset(allowed)


def kmer_search(sequence, kmer_length):
    """
    Driver for get_kmer_count
    :param sequence:
    :param kmer_length:
    :return:
    """
    counts = Counter()
    for i in range(len(sequence) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = sequence[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq] += 1
    return counts


def get_vcf_info_fields(vcf_path):
    vcf = VCF(vcf_path)
    header = vcf.raw_header.split('\n')
    info = [x for x in header if '##INFO' in x]
    # NOTE: this is harcoded for gnomAD vcf
    return [re.split('[<>,=]', x)[3] for x in info]


def get_kmer_count(sequence, kmer_length, nprocs=None):
    """
    Counts the number of k-mers of a specified length in a given sequence
    :param sequence: string representation of a sequence
    :param kmer_length: desired length of k-mer
    :param nprocs: number of processors to use
    :return: collections.Counter (dictionary) object indexed by k-mer and mapped to number of occurences
    """
    start = time.time()
    if nprocs is None:
        nprocs = mp.cpu_count()
    args = split_seq(sequence, nprocs, overlap=kmer_length)
    args = [[seq, kmer_length] for seq in args]
    pool = mp.Pool(mp.cpu_count())
    results = [res.get() for res in [pool.starmap_async(kmer_search, args)]]
    pool.close()
    counts = Counter()
    for result in results[0]:
        for k, v in result.items():
            counts[k] += v
    print("Done in " + str(time.time() - start))
    return counts


def process_region(region, vcf_path, ref_fasta, kmer_size, singletons=True, nsingletons=False):
    """
    Process a VCFRegion to look for variants and sequence context
    :param region: VCFRegion as defined in kclass.py
    :param vcf_path: path to VCF file
    :param ref_fasta: path to reference genome in '.fa' format
    :param kmer_size: integer of total length of k-mer
    :param singletons: if True, will include singleton variants
    :param nsingletons: if True, will include nonsingleton variants
    :return: a dictionary containing the type of result mapped to the data structure
    """
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    # Open iterators
    vcf = VCF(vcf_path)
    ref = Fasta(ref_fasta)

    if singletons and not nsingletons:
        s_transitions = defaultdict(Counter)
        s_positions = defaultdict(Variant)
        for section in region:
            start = time.time()
            print('Processing ' + str(section))
            for variant in vcf(str(section)):
                if is_quality_singleton(variant):
                    new_var = Variant(variant=variant, fields=['vep'])
                    s_positions[new_var.INDEX] = new_var
                    # take 7mer around variant. pyfaidx excludes start index and includes end index
                    adj_seq = ref[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
                    if complete_sequence(adj_seq):
                        s_transitions[adj_seq][new_var.ALT[0]] += 1
            print('Time to complete section ' + str(section) + ': ' + str(time.time() - start))
        return {'singleton_transitions': s_transitions, 'singleton_positions': s_positions}

    elif nsingletons and not singletons:
        ns_transitions = defaultdict(Counter)
        ns_positions = defaultdict(Variant)
        for section in region:
            start = time.time()
            print('Processing ' + str(section))
            for variant in vcf(str(section)):
                if is_quality_nonsingleton(variant):
                    new_var = Variant(variant=variant, fields=['vep'])
                    ns_positions[new_var.INDEX] = new_var
                    # take 7mer around variant. pyfaidx excludes start index and includes end index
                    adj_seq = ref[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
                    if complete_sequence(adj_seq):
                        ns_transitions[adj_seq][new_var.ALT[0]] += 1
            print('Time to complete section ' + str(section) + ': ' + str(time.time() - start))
        return {'nonsingleton_transitions': ns_transitions, 'nonsingleton_positions': ns_positions}

    elif singletons and nsingletons:
        s_transitions = defaultdict(Counter)
        s_positions = defaultdict(Variant)
        ns_transitions = defaultdict(Counter)
        ns_positions = defaultdict(Variant)
        for section in region:
            start = time.time()
            print('Processing ' + str(section))
            for variant in vcf(str(section)):
                if is_quality_singleton(variant):
                    new_var = Variant(variant=variant,  fields=['vep'])
                    s_positions[new_var.INDEX] = new_var
                    # take 7mer around variant. pyfaidx excludes start index and includes end index
                    adj_seq = ref[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
                    if complete_sequence(adj_seq):
                        s_transitions[adj_seq][new_var.ALT[0]] += 1
                if is_quality_nonsingleton(variant):
                    new_var = Variant(variant)
                    ns_positions[new_var.INDEX] = new_var
                    # take 7mer around variant. pyfaidx excludes start index and includes end index
                    adj_seq = ref[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
                    if complete_sequence(adj_seq):
                        ns_transitions[adj_seq][new_var.ALT[0]] += 1
            print('Time to complete section ' + str(section) + ': ' + str(time.time() - start))
        return {'singleton_transitions': s_transitions, 'nonsingleton_transitions': ns_transitions,
                'singleton_positions': s_positions, 'nonsingleton_positions': ns_positions}

    else:
        print('No variants tested. Please have at least one keyword \'singletons\' or \'nsingletons\' set to \'True\'')
        return


def get_kmer_context(vcf_path, ref_fasta, kmer_length, nprocs=0, singletons=True, nsingletons=False):
    """
    Process VCF in parallel to obtain specified variants and their k-mer contexts
    :param vcf_path: Path to VCF file
    :param ref_fasta: Path to reference Fasta ('.fa')
    :param kmer_length: Integer for k-mer size (eg. 3 will give nucleotides 1 position to the left and right)
    :param nprocs: Number of CPUs to use
    :param singletons: True/False process singletons?
    :param nsingletons: True/False proces nonsingletons? (NOTE: These 2 can be used together)
    :return: A list containing dictionaries from each process. Maps title of information to the data structure. Post-processing will be a future feature.
    """
    if singletons and not nsingletons:
        print('Processing %s singletons' % vcf_path)
    elif nsingletons and not singletons:
        print('Processing %s nonsingletons' % vcf_path)
    elif singletons and nsingletons:
        print('Processing %s singletons and nonsingletons' % vcf_path)
    else:
        print('Not processing anything. Please have at least one keyword \'singletons\' or \'nsingletons\' set to \'True\'')

    if nprocs == 0:
        nprocs = mp.cpu_count()
    # intentionally create more chunks to even CPU distribution
    regions = get_split_vcf_regions(vcf_path, nprocs * 20)  # TODO: This number (20) can and should change to optimize
    arguments = [(region, vcf_path, ref_fasta, kmer_length, dict(singletons=singletons, nsingletons=nsingletons)) for
                 region in regions]
    pool = mp.Pool(nprocs)
    results = [result.get() for result in [pool.starmap_async(process_region, arguments)]]
    pool.close()
    return results[0]


def merge_transitions_ddc(dict_list, outfile=None):
    """
    Merge a list of defaultdict(Counter)'s into a single data structure
    :param outfile: name of file to write (leave empty or ignore to not write to disk)
    :param dict_list: type = iterable containing defaultdict(Counter)
    :return: a single defaultdict(Counter)
    """
    print('Writing kmer transitions to file')
    master_count = defaultdict(Counter)
    for counts in dict_list:
        for k, v in counts.items():
            for alt, count in v.items():
                master_count[k][alt] += count
    if outfile is not None and len(outfile) > 1:
        pd.DataFrame.from_dict(master_count, orient='index').to_csv(outfile)
    return master_count


def merge_positions_dd(positions_iter, outfile='positions.bed'):
    """
    Write variant positions to a bed-style file
    :param positions_iter: iterable containing defaultdict(Variant)'s
    :param outfile: name of merged file
    :return: not a damn thing (except your file should be on your filesystem)
    """
    print('Writing variant positions to file: %s' % outfile)
    output = open(outfile, "a+")
    output.write("CHROM\tPOS\tREF\tALT\tAC\n")  # header same for all
    for pos in positions_iter:
        for k, v in pos.items():
            output.write(v.print_variant(fields=['vep']))
    output.close()
    return
