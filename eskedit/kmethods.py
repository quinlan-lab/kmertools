import random
import multiprocessing as mp
import sys
import time
from collections import Counter, defaultdict
import numpy as np
from cyvcf2 import VCF
import re
import math
import pandas as pd
from pkg_resources import resource_filename
from pyfaidx import Fasta, FetchError

from eskedit.constants import get_autosome_names_grch38
from eskedit.kclass import Variant, Kmer, KmerWindow, GRegion
from eskedit.ksplit import split_seq, get_split_vcf_regions, get_split_chrom_vcf


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


def ref_genome_as_string(ref_fasta, keys=None):
    ref_genome = Fasta(ref_fasta)
    if keys is None:
        keys = ref_genome.keys()
    ref_seq = ""
    for chrom in keys:
        ref_seq += str(ref_genome[chrom])
    return ref_seq


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


def is_quality_snv(var_to_test, AC_cutoff=None):
    """
    high quality variants will have FILTER == None
    AND we are ignoring insertions and deltions here
    """
    if AC_cutoff is not None:
        try:
            AC_cutoff = int(AC_cutoff)
            return var_to_test.FILTER is None and var_to_test.INFO.get(
                'variant_type') == 'snv' and var_to_test.INFO.get('AC') <= AC_cutoff
        except ValueError:
            AC_cutoff = None
    return var_to_test.FILTER is None and var_to_test.INFO.get('variant_type') == 'snv'


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
    NOTE: This method is now deprecated, please use "is_singleton_snv"
    high quality variants will have FILTER == None
    Additionally, variants shoud be singletons ('AC' == 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') == 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def is_singleton_snv(var_to_test):
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') == 1 and var_to_test.INFO.get(
        'variant_type') == 'snv'


def complete_sequence(seq):
    """
    Seq is normalized by capitalizing all of the characters
    :param seq: sequence to test
    :return: True if sequence only contains A, C, T, or G, False if contains N or any other character
    """
    allowed = set('ACTG')
    return set(seq.upper()).issubset(allowed)


def kmer_search(sequence, kmer_length, count_gc=False, count_n=False):
    """
    Driver for get_kmer_count
    :param sequence:
    :param kmer_length:
    :param count_gc: Return GC content (default = False)
    :param count_n: Reurn count number of unreferenced nucleotides (default = False)
    :return:
    """
    counts = Counter()
    for i in range(len(sequence) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = sequence[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq.upper()] += 1
    if count_gc and count_n:
        nucleotides = Counter(sequence)
        gc_content = (nucleotides['G'] + nucleotides['C']) / max(
            (nucleotides['A'] + nucleotides['T'] + nucleotides['G'] + nucleotides['C']), 1)
        return counts, gc_content, nucleotides['N']
    elif count_gc and not count_n:
        nucleotides = Counter(sequence)
        gc_content = (nucleotides['G'] + nucleotides['C']) / max((
                nucleotides['A'] + nucleotides['T'] + nucleotides['G'] + nucleotides['C']), 1)
        return counts, gc_content
    elif not count_gc and count_n:
        nucleotides = Counter(sequence)
        return counts, nucleotides['N']
    else:
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
    pool = mp.Pool(nprocs)
    results = [res.get() for res in [pool.starmap_async(kmer_search, args)]]
    pool.close()
    counts = Counter()
    for result in results[0]:
        for k, v in result.items():
            counts[k] += v
    # print("Done in " + str(time.time() - start))
    return counts


def process_region(region, vcf_path, ref_fasta, kmer_size, singletons, nsingletons):
    """
    Process a GRegion to look for variants and sequence context
    :param region: GRegion as defined in kclass.py
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
        print('Processing only singleton variants')
        s_transitions = defaultdict(Counter)
        s_positions = defaultdict(Variant)
        for section in region:
            start = time.time()
            print('Processing ' + str(section))
            for variant in vcf(str(section)):
                if is_singleton_snv(variant):
                    new_var = Variant(variant=variant, fields=['vep'])
                    s_positions[new_var.INDEX] = new_var
                    # take 7mer around variant. pyfaidx excludes start index and includes end index
                    adj_seq = ref[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
                    if complete_sequence(adj_seq):
                        s_transitions[adj_seq.upper()][new_var.ALT[0]] += 1
            print('Time to complete section ' + str(section) + ': ' + str(time.time() - start))
        return {'singleton_transitions': s_transitions, 'singleton_positions': s_positions}

    elif nsingletons and not singletons:
        print('Procesing only non-singleton variants')
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
        print('Processing singleton and non-singleton variants')
        s_transitions = defaultdict(Counter)
        s_positions = defaultdict(Variant)
        ns_transitions = defaultdict(Counter)
        ns_positions = defaultdict(Variant)
        for section in region:
            start = time.time()
            print('Processing ' + str(section))
            for variant in vcf(str(section)):
                if is_singleton_snv(variant):
                    new_var = Variant(variant=variant, fields=['vep'])
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
            print('Time to complete section ' + str(section) + ': ' + str(time.time() - start), flush=True)
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
        print(
            'Not processing anything. Please have at least one keyword \'singletons\' or \'nsingletons\' set to \'True\'')

    if nprocs == 0:
        nprocs = mp.cpu_count()
    # intentionally create more chunks to even CPU distribution
    regions = get_split_vcf_regions(vcf_path, nprocs * 20)  # TODO: This number (20) can and should change to optimize
    arguments = [(region,
                  vcf_path,
                  ref_fasta,
                  kmer_length,
                  singletons,
                  nsingletons) for region in regions]
    pool = mp.Pool(nprocs)
    results = [result.get() for result in [pool.starmap_async(process_region, arguments)]]
    pool.close()
    pool.join()
    return results[0]


def merge_transitions_ddc(dict_list, outfile=None):
    """
    Merge a list of defaultdict(Counter)'s into a single data structure
    :param outfile: name of file to write (leave empty or ignore to not write to disk)
    :param dict_list: type = iterable containing defaultdict(Counter)
    :return: a single defaultdict(Counter)
    """
    master_count = defaultdict(Counter)
    for counts in dict_list:
        for k, v in counts.items():
            for alt, count in v.items():
                master_count[k][alt] += count
    if outfile is not None and len(outfile) > 1:
        print('Writing kmer transitions to file %s' % outfile)
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


def combine_df_kmer_indices(df):
    combined = Counter()
    df.fillna(0, inplace=True)
    for i, r in df.iterrows():
        try:
            k = Kmer(i)
        except (ValueError, KeyError):
            print("Row %s was skipped" % str(r.name), file=sys.stderr)
            continue
        combined[k] += r
    return pd.DataFrame.from_dict(combined, orient='index').sort_index()


def clean_counts_df(counts_fpath, get_freq=True, combine_complements=False, save_file=None):
    df = pd.read_csv(counts_fpath, index_col=0)
    df.fillna(0, inplace=True)
    cdict = defaultdict(Counter)
    for i, r in df.iterrows():
        if combine_complements:
            k = Kmer(i.upper())
        else:
            k = i.upper()
        for n in list('ACGT'):
            cdict[k][n] += r[n]
    df = pd.DataFrame.from_dict(cdict, orient='index')
    df.sort_index(inplace=True)
    if get_freq:
        df['total'] = df.sum(axis=1)
        df = df / df.sum().sum()
    if save_file is not None:
        df.to_csv(save_file)
    return df


def file_len(fname):
    """
    Taken from https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
    :param fname: path to file
    :return: number of lines in the file
    """
    with open(fname, 'r') as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def get_counts_dict(kmer_size, alt_path=None):
    if alt_path is None:
        filepath = resource_filename('eskedit',
                                     '../input_data/counts_data/{}mer_relative_mutation_freq_v3.csv'.format(kmer_size))
        # path = './input_data/counts_data/{}mer_relative_mutation_freq_v3.csv'.format(kmer_size)
    else:
        filepath = alt_path
    df = pd.read_csv(filepath, index_col=0)
    # df['probability'] = df.iloc[:, :4].apply
    return dict(zip(df.index, df.frequency))


def count_regional_variants(vcf_region):
    singletons = 0
    count = 0
    for v in vcf_region:
        if is_singleton_snv(v):
            singletons += 1
        if is_quality_snv(v):
            count += 1
    return count, singletons


def count_regional_AF(vcf_region):
    total = 0
    calc = 0
    for v in vcf_region:
        if is_quality_snv(v):
            total += v.INFO.get('AF')
            calc += v.INFO.get('AC') / v.INFO.get('AN')
    return calc, total


def query_region(vcf_path, fasta, chrom, kmer_size, bins=100, counts_path=None):
    # TODO: Add window size!!!
    vcf = VCF(vcf_path)
    fasta = Fasta(fasta)
    regions = get_split_chrom_vcf(vcf_path, chrom, bins)
    window = KmerWindow(kmer_size, counts_path=counts_path)
    print('region\tactual\texpected\tratio', flush=True)
    expected, actual = [], []
    for r in regions:
        try:
            exp = window.calculate_expected(str(fasta.get_seq(r.chrom, r.start, r.stop)))
        except (KeyError, FetchError):
            exp = 0
        all_vars, singletons = count_regional_variants(vcf(str(r)))
        expected.append(exp)
        actual.append(singletons)
        if exp == 0:
            ratio = 0
        else:
            ratio = singletons / exp
        print("%s\t%d\t%f\t%f" % (str(r), all_vars, exp, ratio), flush=True)


def query_bed_region(region, vcf_path, fasta, kmer_size, counts_path, count_frequency):
    """
    @param region:
    @param vcf_path:
    @param fasta:
    @param kmer_size:
    @param counts_path:         This field is critical for count_freuency to work. Needs table of expected AF
    @param count_frequency:
    @return:
    """
    # TODO: Add binning somehow (either keep equal size or equal number of bins
    start = time.time()
    vcf = VCF(vcf_path)
    fasta = Fasta(fasta)
    window = KmerWindow(kmer_size, counts_path=counts_path)
    # expected, actual = [], []
    try:
        if region.strand is not None:
            if is_dash(region.strand):
                sequence = fasta.get_seq(region.chrom, region.start, region.stop).complement.seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start, region.stop).seq.upper()
        else:
            sequence = fasta.get_seq(region.chrom, region.start, region.stop).seq.upper()
        exp = window.calculate_expected(sequence)  # this does account for strandedness
        if count_frequency:
            calc, total = count_regional_AF(vcf(str(region)))
            if not math.isclose(calc, total, rel_tol=1e-05):
                print('WARNING: Calculated AF and VCF AF are different!     Calculated AF: %f     VCF AF: %f' % (
                    calc, total), file=sys.stderr, flush=True)
            observed_variants = calc
            # Setting this so output shows difference between calculated and listed AF and won't throw an error
            all_vars = total
            pass
        else:
            # does not account for strandedness here
            all_vars, observed_variants = count_regional_variants(vcf(str(region)))
    except (KeyError, FetchError):
        exp = 0
        observed_variants = 0
        all_vars = 0
    if exp == 0:
        ratio = 0
    else:
        ratio = observed_variants / exp
    print('{0:<30} {1:>10} {2:>20} {3:>20} {4:>20}'.format((region.printstr(delim=' ')),
                                                           str(all_vars - observed_variants), str(observed_variants),
                                                           str(exp), str(ratio)), flush=True)
    # return "%s\t%s\t%s\t%d\t%f\t%f\n" % (str(region.chrom), str(region.start), str(region.stop), act, exp, ratio)
    return "%s\t%f\t%f\t%f\t%f\n" % (region.printstr(), all_vars - observed_variants, observed_variants, exp, ratio)


def check_bed_regions(bed_path, vcf_path, fasta_path, kmer_size, nprocs=4, counts_path=None, outfile=None,
                      strand_col=None, bed_names_col=None, singletons=False):
    # TODO: implement where bed regions are read and then analyzed for expected v actual
    try:
        kmer_size = int(kmer_size)
        nprocs = int(nprocs)
    except ValueError:
        print('ERROR: kmer_size and nprocs must be integers', file=sys.stderr, flush=True)
        exit(1)
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
    regions = []
    if outfile is None:
        outfile = '{}mer_expected_mutation_count.dat'.format(kmer_size)
    with open(outfile, 'w') as output:
        output.write('Region\tObserved\tExpected\tObsToExp\n')
    with open(bed_path, 'r') as bedfile:
        kwargs = defaultdict(str)
        for line in bedfile.readlines():
            fields = line.split('\t')
            for k, v in additional_fields.items():
                kwargs[k] = fields[v]
            regions.append(GRegion(*[fields[0], fields[1], fields[2]], **kwargs))
    print('{0:<30} {1:>10} {2:>20} {3:>20}'.format('Region', 'Observed', 'Expected', 'Obs/Exp'), flush=True)
    # Expected arguments query_bed_region(region, vcf_path, fasta, kmer_size, counts_path, count_frequency)

    arguments = [(region, vcf_path, fasta_path, kmer_size, counts_path, (not singletons)) for region in regions]
    pool = mp.Pool(nprocs)
    results = pool.starmap_async(query_bed_region, arguments)
    pool.close()
    pool.join()
    with open(outfile, 'w') as output:
        for result in results.get():
            output.write(result)
    pass


def is_dash(pdash):
    regex = '[\u002D\u058A\u05BE\u1400\u1806\u2010-\u2015\u2E17\u2E1A\u2E3A\u2E3B\u2E40\u301C\u3030\u30A0\uFE31\uFE32\uFE58\uFE63\uFF0D]'
    if re.match(regex, pdash) is not None:
        return True
    return False


def check_clinvar(vcf_path, fasta_path, kmer_size, left_context=0, right_context=0, counts_path=None):
    # VCF is 1-based, closed
    # BED is 0-based, half-open
    names = get_autosome_names_grch38()
    vcf = VCF(vcf_path)
    fasta = Fasta(fasta_path, read_ahead=10_000_000)
    window = KmerWindow(kmer_size, counts_path=counts_path)
    clinvar = []
    for variant in vcf:
        if variant.CHROM not in names: continue
        start = variant.POS - left_context
        stop = variant.POS + right_context
        # start = variant.POS - (kmer_size)
        # stop = variant.POS + (kmer_size)
        seq = fasta.get_seq(variant.CHROM, start, stop).seq.upper()
        clinvar.append((variant.CHROM, variant.POS, variant.INFO.get("CLNSIG"),
                        window.calculate_expected(seq), seq))
    return clinvar


if __name__ == "__main__":
    print("Testing from kmethods.")
    starttime = time.time()
    vp = '/Users/simonelongo/Downloads/chr_clinvar.hg38.bed.gz'
    # vp = '/Users/simonelongo/Downloads/clinvar_chr22.vcf.gz'
    fp = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
    cv = check_clinvar(vp, fp, 3, left_context=5, right_context=5)
    df = (pd.DataFrame(cv))
    df.columns = ['chrom', 'pos', 'clnsig', 'expected', 'kmer']
    df.to_csv('5left5right_clinvar_3mers.csv', index=False)

    print('Done in %s' % str(time.time() - starttime))
