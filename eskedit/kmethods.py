import math
import os
import random
import multiprocessing as mp
import sys
import time
from collections import Counter, defaultdict
import numpy as np
from cyvcf2 import VCF
import re
import pandas as pd
from pkg_resources import resource_filename
from pyfaidx import Fasta, FetchError
import datetime
from eskedit.constants import get_autosome_names_grch38
from eskedit.kclass import Kmer, KmerWindow, GRegion
from eskedit.ksplit import split_seq
import itertools
import matplotlib.pyplot as plt


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


def generate_kmers(k):
    """Generates a list of all possible DNA sequences of length k. E.g. generate_kmers(2) will return
    [ AA, AC, AT, AG, CA, CC, CT, CG, TA, TC, TT, TG, GA, GC, GT, GG ] """
    len_k = int(k)
    if len_k < 0:
        raise ValueError("Must be a positive integer")
    combos = list(itertools.product('ACTG', repeat=len_k))
    seqs = []
    for seq in combos:
        seqs.append(''.join(seq))
    return seqs


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


def get_counts_from_file(filepath):
    if filepath is None:
        print('ERROR: No input file detected.', file=sys.stderr, flush=True)
        exit(1)
    countsdf = pd.read_csv(filepath, index_col=0)
    # df['probability'] = df.iloc[:, :4].apply
    return dict(zip(countsdf.index, countsdf.iloc[:, 0]))


def get_counts_dict(kmer_size, name):
    try:
        filepath = resource_filename('eskedit',
                                     '../input_data/counts_data/{}mer_{}_table.csv'.format(kmer_size, name))
        countsdf = pd.read_csv(filepath, index_col=0)
        # df['probability'] = df.iloc[:, :4].apply
        return dict(zip(countsdf.index, countsdf.iloc[:, 0]))
    except FileNotFoundError:
        return None


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
    AF = 0.0
    AN = 0
    AC = 0
    for v in vcf_region:
        if is_quality_snv(v):
            try:
                AF += float(v.INFO.get('AF'))
                AC += int(v.INFO.get('AC'))
                AN += int(v.INFO.get('AN'))
            except ValueError:
                AF = 0
                AC = 0
                AN = 0
                print("ERROR: casting AC, AF, or AN in %s" % (str(vcf_region)), file=sys.stderr, flush=True)
    return AF, AC, AN


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
    # The first kmer actually begins centered around first nucleotide in sequence so
    # start position is shifted upstream by half the kmer length
    # end position is shifted downstream by the same
    shift = kmer_size // 2
    try:
        if region.strand is not None:
            if is_dash(region.strand):
                sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).complement.seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq.upper()
        else:
            sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq.upper()
        exp = window.calculate_expected(sequence)  # this does account for strandedness
        if count_frequency:
            AF, AC, AN = count_regional_AF(vcf(str(region)))
            # if not math.isclose(calc, total, rel_tol=1e-05):
            #     print('WARNING: Calculated AF and VCF AF are different!     Calculated AF: %f     VCF AF: %f' % (
            #         calc, total), file=sys.stderr, flush=True)
            field1 = AC
            field2 = AN
            field3 = AF
            field4 = exp
            # if exp == 0:
            #     field4 = exp
            # else:
            #     field4 = exp
        else:
            # does not account for strandedness here
            all_vars, observed_variants = count_regional_variants(vcf(str(region)))
            field1 = all_vars - observed_variants
            field2 = observed_variants
            field3 = exp
            if exp == 0:
                field4 = 0
            else:
                field4 = observed_variants / exp
    except (KeyError, FetchError):
        field1 = 0
        field2 = 0
        field3 = 0
        field4 = 0
    #     exp = 0
    #     observed_variants = 0
    #     all_vars = 0
    # if exp == 0:
    #     ratio = 0
    # else:
    #     ratio = observed_variants / exp
    print('{0:<30} {1:>10} {2:>20} {3:>20} {4:>20}'.format((region.printstr(delim=' ')),
                                                           str(field1), str(field2),
                                                           str(field3), str(field4)), flush=True)
    return "%s\t%s\t%s\t%s\t%s\n" % (region.printstr(), str(field1), str(field2), str(field3), str(field4))


def check_bed_regions(bed_path, vcf_path, fasta_path, kmer_size, nprocs=4, counts_path=None, outfile=None,
                      strand_col=None, bed_names_col=None, singletons=False):
    # TODO: implement where bed regions are read and then analyzed for expected v actual

    # field1 = AC
    # field2 = AN
    # field3 = AF
    # field1 = all_vars - observed_variants
    # field2 = observed_variants
    # field3 = exp
    if singletons:
        field1 = 'NotSingletons'
        field2 = 'Singletons'
        field3 = 'Expected'
        field4 = 'ObsExpRatio'
    else:
        field1 = 'AC'
        field2 = 'AN'
        field3 = 'AF'
        field4 = 'ExpectedAF'

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
        outfile = bed_path.split('/')[-1]
        outfile = "".join(outfile.split('.')[:-1])
        outfile = '%s_%s.%dmer.dat' % (outfile, str(datetime.date.today()), kmer_size)
    with open(outfile, 'w') as output:
        output.write('%s\t%s\t%s\t%s\t%s\n' % ('Region', field1, field2, field3, field4))
    with open(bed_path, 'r') as bedfile:
        kwargs = defaultdict(str)
        for line in bedfile.readlines():
            fields = line.split('\t')
            for k, v in additional_fields.items():
                kwargs[k] = fields[v]
            regions.append(GRegion(*[fields[0], fields[1], fields[2]], **kwargs))
    print('{0:<30} {1:>10} {2:>20} {3:>20} {4:>20}'.format('Region', field1, field2, field3, field4), flush=True)
    # Expected arguments query_bed_region(region, vcf_path, fasta, kmer_size, counts_path, count_frequency)

    arguments = [(region, vcf_path, fasta_path, kmer_size, counts_path, (not singletons)) for region in regions]
    pool = mp.Pool(nprocs)
    results = pool.starmap_async(query_bed_region, arguments)
    pool.close()
    pool.join()
    with open(outfile, 'w') as output:
        for result in results.get():
            output.write(result)
    return


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


def get_model_paths(modeldir):
    data = {'singleton': None, 'AF': None, 'AN': None, 'AC': None, 'reference': None}
    for file in os.listdir(modeldir):
        for k in data.keys():
            if k in file.split('_'):
                data[k] = (os.path.join(modeldir, file))
    return data


def _buildtable(path, counts):
    ts = pd.read_csv(path, index_col=0).sort_index()
    return pd.DataFrame(ts.mean(axis=1) / counts['0'])


def build_tables(modeldir):
    data = get_model_paths(modeldir)
    tables = {key: None for key in data.keys()}
    tables['reference'] = pd.read_csv(data['reference'], index_col=0).sort_index()
    kmersize = int(math.log(len(tables['reference'].index), 4))
    dirname = 'TABLES_{}'.format(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))
    dirname = os.path.join(modeldir, dirname)
    os.makedirs(dirname)
    for k, path in data.items():
        if tables[k] is None:
            df = _buildtable(path, tables['reference'])
            df.to_csv(os.path.join(dirname, '{}mer_{}_table.csv'.format(kmersize, k)))
            tables[k] = df
    return tables


def build_7mer_flanks(df, ref=None, alt=None, rcomp=False):
    rvcomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    nucs = list('ACGT')
    mindex = list(it.product(nucs, nucs, nucs))
    dindex = pd.MultiIndex.from_tuples(mindex)
    flanks = pd.DataFrame(index=dindex, columns=dindex)
    flanks = flanks.fillna(0)
    ct = []
    condition = [ref is None, alt is None]
    for i, r in df.iterrows():
        if ref is None:
            ref = i[3]
        if alt is None:
            alt = 0
        if i[3] == ref:
            ct.append((i, r[alt]))
        if rcomp:
            if i[3] == rvcomp[ref]:
                ct.append((i, r[rvcomp[alt]]))
    ct = np.array(ct)
    ctdf = pd.DataFrame.from_records(ct)
    ctdf.columns = ['Kmer', 'Count']
    ctdf.Count = ctdf.Count.astype('float64')
    for i, val in zip(list(ctdf.Kmer), list(ctdf.Count)):
        idx = tuple(list(i[:3]))
        col = tuple(list(i[4:]))
        flanks.loc[idx, col] = val
    return flanks.astype('float64')


def plot_7mer_heatmap(flanks):
    fontsize = 30
    spacing = fontsize + 5
    nucs = list('ACGT')

    plt.rcParams['figure.figsize'] = [30, 30]
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize

    level1 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 4)]
    midpt = (level1[1] - level1[0]) / 2
    level1 = [i + midpt for i in level1]
    l1_lab = nucs

    level2 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 16)]
    midpt = (level2[1] - level2[0]) / 2
    level2 = [i + midpt for i in level2]
    l2_lab = nucs * 4

    level3 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 64)]
    midpt = (level3[1] - level3[0]) / 2
    level3 = [i + midpt for i in level3]
    l3_lab = nucs * 16

    ax1 = plt.subplot(111)
    ax1.matshow(flanks.T)
    ax1.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax1.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_label_position('left')
    ax1.set_xticks(level3)
    ax1.set_xticklabels(l3_lab)
    ax1.set_yticks(level3)
    ax1.set_yticklabels(l3_lab)

    ax2 = ax1.twinx()  # .twinx()
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', spacing))
    ax2.set_xticks(level2)
    ax2.set_xticklabels(l2_lab)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_yticks(level2)
    ax2.set_yticklabels(l2_lab)
    ax2.spines['left'].set_position(('outward', spacing))
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.set_ylim(ax1.get_ylim())

    ax3 = ax2.twinx()  # .twinx()
    ax3.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax3.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax3.spines['bottom'].set_position(('outward', spacing * 2))
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')
    ax3.spines['left'].set_position(('outward', spacing * 2))
    ax3.set_xticks(level1)
    ax3.set_xticklabels(l1_lab)
    ax3.set_yticks(level1)
    ax3.set_yticklabels(l1_lab)
    ax3.set_ylim(ax1.get_ylim())


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
