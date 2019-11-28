from collections import Counter
import pandas as pd
import multiprocessing as mp
from pyfaidx import Fasta
from eskedit import split_seq
import eskedit as ek


def get_clean_seq(seq):
    """
    Collapse long stretches of N in a sequence
    :param seq: string or iterable representation of the sequence
    :return: tuple(a list containing the original sequence with 'N's collapsed and counted, total number of 'N's)
    """
    total_n = 0
    ncount = 0
    new_seq = []
    for i, n in enumerate(seq):
        if n.upper() == 'N':
            ncount += 1
            continue
        else:
            if ncount > 0:
                new_seq.append(str(ncount) + 'N')
                total_n += ncount
                ncount = 0
            new_seq.append(n)
    return new_seq, total_n


def section_seq_prep(seq, kmer_size):
    total_n = 0
    ncount = 0
    cg_count = 0
    new_seq = []
    for i, n in enumerate(seq):
        if n.upper() == 'N':
            ncount += 1
            continue
        else:
            if ncount > 0:
                new_seq.append(str(ncount) + 'N')
                total_n += ncount
                ncount = 0
            if n.upper() == 'G':
                if seq[i - 1].upper() == 'C':
                    cg_count += 1
            new_seq.append(n)
    return new_seq, cg_count, total_n


def dpgp_kmer_search(sequence, kmer_length):
    """
    Driver for get_kmer_count
    :param sequence:
    :param kmer_length:
    :return:
    """
    total_n = 0
    cg_count = 0
    previous_n_idx = 0
    counts = Counter()
    for i, n in enumerate(sequence):
        if i >= (len(sequence) - (kmer_length - 1)):
            break
        next_seq = sequence[i:(i + kmer_length)]
        if sequence[i] == 'G' and sequence[i - 1] == 'C':
            cg_count += 1
        nseq_str = "".join(next_seq).upper()
        if 'N' not in nseq_str:
            counts[nseq_str] += 1
        else:
            new_n_idx = 0
            for idx in range(i - (kmer_length // 2), i + (kmer_length // 2)):
                if 'N' in sequence[idx]:
                    new_n_idx = idx
                    if previous_n_idx < new_n_idx - (kmer_length // 2):
                        try:
                            total_n += int(sequence[new_n_idx].split('N')[0])
                        except ValueError:
                            total_n += 1  # these are simply N's
            previous_n_idx = new_n_idx

    # for i in range(len(sequence) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
    #     next_seq = sequence[i:(i + kmer_length)]
    #     if not ('N' in next_seq or 'n' in next_seq):
    #         counts[next_seq.upper()] += 1
    return counts, total_n, cg_count


def chrom_stats(seq, kmer_size, nbins, chrom=None):
    df = pd.DataFrame()
    genome_position = 0
    clean_seq, total_n = get_clean_seq(seq)
    chrom_split = split_seq(clean_seq, nbins)
    for idx, section in enumerate(chrom_split):
        counts, chrom_n, cg_count = dpgp_kmer_search(section, kmer_size)
        chrom_df = pd.DataFrame.from_dict(counts, orient='index', columns=[str(genome_position)])
        genome_position += len(section)
        chrom_df.loc['CG_count'] = cg_count
        chrom_df.loc['section_N_count'] = chrom_n
        chrom_df.loc['total_N_count'] = total_n
        if idx == 0:
            df = chrom_df
        else:
            df = pd.concat([df, chrom_df], axis=1, sort=False)
    return df


def named_chrom_stats(chromname, seq, kmer_size, nbins):
    df = pd.DataFrame()
    genome_position = 0
    clean_seq, total_n = get_clean_seq(seq)
    chrom_split = split_seq(clean_seq, nbins)
    for idx, section in enumerate(chrom_split):
        counts, chrom_n, cg_count = dpgp_kmer_search(section, kmer_size)
        chrom_df = pd.DataFrame.from_dict(counts, orient='index', columns=[str(chromname) + '_' + str(genome_position)])
        genome_position += len(section)
        chrom_df.loc['CG_count'] = cg_count
        chrom_df.loc['section_N_count'] = chrom_n
        chrom_df.loc['total_N_count'] = total_n
        if idx == 0:
            df = chrom_df
        else:
            df = pd.concat([df, chrom_df], axis=1, sort=False)
    return df


def kmers_per_chromosome_clustering(ref_fasta, kmer_size, chrom_bins, nprocs=None):
    if nprocs is None:
        nprocs = mp.cpu_count()
    fa = Fasta(ref_fasta)
    chrom_seqs = []
    for chrom in ek.get_autosome_names_grch38():
        chrom_seqs.append((chrom, str(fa[chrom])))
    args = [(*seq, kmer_size, chrom_bins) for seq in chrom_seqs]
    pool = mp.Pool(nprocs)
    results = [result.get() for result in [pool.starmap_async(named_chrom_stats, args)]]
    pool.close()
    df = pd.DataFrame()
    for idx, res in enumerate(results[0]):
        if idx == 0:
            df = res
        else:
            df = pd.concat([df, res], axis=1, sort=False)
    return df
