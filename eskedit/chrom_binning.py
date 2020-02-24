import array
from collections import defaultdict
from cyvcf2 import VCF
from pyfaidx import Fasta, FetchError
import eskedit as ek
import sys
import multiprocessing as mp
from eskedit import GRegion, Variant
import time
import pandas as pd


def bin_chrom_vcf(chrom, chrom_len, nbins):
    bin_size = chrom_len // nbins
    regions = []
    start = 1
    stop = bin_size
    for i in range(nbins):
        if i < nbins - 1:
            regions.append(GRegion(chrom, start, stop))
            start = stop
            stop += bin_size
        else:
            regions.append(GRegion(chrom, start, chrom_len))
    return regions


def process_chrom_bin(region, kmer_size, vcf_path, fasta_path):
    start = time.time()
    fasta = Fasta(fasta_path)
    vcf = VCF(vcf_path)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    try:
        sequence = fasta.get_seq(region.chrom, region.start, region.stop).seq.upper()
    except (KeyError, FetchError):
        print('Region %s not found in fasta, continuing...' % str(region), file=sys.stderr)
        return
    region_ref_counts, gc_content, n_count = ek.kmer_search(sequence, kmer_size, count_gc=True,
                                                            count_n=True)  # nprocs=1 due to short region
    r_string = str(region.chrom) + ':' + str(region.start) + '-' + str(region.stop)
    transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    # Define indices for nucleotides
    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    # count, singletons = ek.count_regional_variants(vcf(r_string))
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
    if len(transitions.keys()) > 0 and len(region_ref_counts.keys()) > 0:
        bin_trans = pd.DataFrame.from_dict(transitions, orient='index')
        bin_trans.sort_index(inplace=True)
        bin_trans['tot'] = bin_trans.sum(axis=1)
        bin_kcounts = pd.DataFrame.from_dict(region_ref_counts, orient='index')
        bin_kcounts.sort_index(inplace=True)
        bin_kcounts.columns = ['counts']
        kmer_freq = pd.concat([bin_trans.loc[:, 'tot'], bin_kcounts], join='outer', axis=1, sort=True)
        kmer_freq.fillna(0, inplace=True)
        kmer_freq['freq'] = kmer_freq.tot / kmer_freq.counts
        kmer_freq.loc['GC_content', 'freq'] = gc_content
        kmer_freq.loc['N_count', 'freq'] = n_count
        print('Finished region %s in %s' % (str(region), str(time.time() - start)), flush=True)
        return region, kmer_freq['freq'].to_dict()
    else:
        print('Finished region %s in %s' % (str(region), str(time.time() - start)), flush=True)
        return region, None


def chrom_bin_mutability(vcfpath, fastapath, kmer_size, nbins, chroms=None, nprocs=1):
    # prepare bins for each chromosome
    try:
        kmer_size = int(kmer_size)
        nbins = int(nbins)
        nprocs = int(nprocs)
    except ValueError:
        print("ERROR: kmer_size, nbins, and nprocs must be integers.")
        exit(1)
    default_idx = 'A' * kmer_size
    if chroms is None:
        chromdict = ek.get_autosome_lengths_grch38()  # dictionary which maps chrom name to length
    else:
        vcf = VCF(vcfpath)
        vcfchroms = dict(zip(vcf.seqnames, vcf.seqlens))
        chromdict = defaultdict(int)
        for c in chroms:
            try:
                chromdict[c] = vcfchroms[c]
            except KeyError:
                print("Chromosome %s not found in VCF. Ignoring..." % str(c), file=sys.stderr)
    chrom_bins = []
    for k, v in chromdict.items():
        chrom_bins.extend(bin_chrom_vcf(k, v, nbins))

    pool = mp.Pool(nprocs)
    arguments = [(region, kmer_size, vcfpath, fastapath) for region in chrom_bins]
    results = pool.starmap(process_chrom_bin, arguments)
    pool.close()
    pool.join()
    genome_df = pd.DataFrame()
    for res in results:
        region, freq = res
        if freq is not None:
            freq_df = pd.DataFrame.from_dict(freq, orient='index', columns=[str(region)])
        else:
            freq_df = pd.DataFrame({str(region): [0]}, index=[default_idx])
        genome_df = pd.concat([genome_df, freq_df], axis=1, sort=True)
    multiindex = pd.MultiIndex.from_tuples([(s.split(':')[0], s.split(':')[1]) for s in genome_df.columns.to_list()],
                                           names=['chrom', 'position'])
    genome_df.columns = multiindex
    genome_df.fillna(0, inplace=True)
    return genome_df


if __name__ == "__main__":
    testvcf, testfasta = ek.get_test_path()
    print(chrom_bin_mutability(testvcf, testfasta, 3, 20, chroms=['chr22'], nprocs=12))
