import pandas as pd
import pathos.multiprocessing as mp
import eskedit as ek
from pyfaidx import Fasta
from cyvcf2 import VCF

"""
This package will search a genome reference sliding window with variants aligned.

This currently works with 3, 5, and 7mer contexts. The expected input is a tab separated
file with the following fields in this order:
CHROM   POS REF ALT`
    1. CHROM - chromosome should be as reported by grch38
    2. POS - position on chromosome aligned to hg38
    3. REF - reference allele
    4. ALT - alternate allele from VCF
    * Additional fields after this will be ignored
"""
WINDOW_SIZE = 0


def process_chunk(df):
    ref = df.iloc[:, 0].to_numpy()
    alt = df.iloc[:, 1].to_numpy()


def sliding_window_from_csv(var_csv, window_size, nprocs=None):
    global WINDOW_SIZE
    WINDOW_SIZE = window_size
    if nprocs is None:
        nprocs = mp.cpu_count()
    chunk_size = ek.file_len(var_csv) // nprocs
    df_chunks = pd.read_csv(var_csv, sep='\t', chunksize=chunk_size)
    pool = mp.Pool(nprocs)
    results = [result.get() for result in [pool.starmap_async(process_chunk, df_chunks)]]
    pool.close()
    pass


def getcountsdict(kmer_size):
    path = './input_data/counts_data/{}mer_relative_mutation_freq_v3.csv'.format(kmer_size)
    df = pd.read_csv(path, index_col=0)
    # df['probability'] = df.iloc[:, :4].apply
    return dict(zip(df.index, df.frequency))


class KmerWindow:
    def __init__(self, kmer_size, gnomad_samples=71702, test_samples=71702):
        if kmer_size not in [3, 5, 7]:
            print('Only supports kmer sizes of 3, 5, and 7 right now.')
            exit(0)
        # self.window_size = window_size
        self.gnomad_chroms = gnomad_samples * 2
        self.test_chroms = test_samples * 2
        self.kmer_size = kmer_size
        self.counts_dict = getcountsdict(kmer_size)

    def calculate_expected(self, seq, raw_data=False):
        kmer_count = 0
        freq_sum = 0.0
        seq = seq.upper()
        # num_nucs = len(seq.replace('N', ''))
        for start in range(len(seq) - self.kmer_size + 1):
            next_k = seq[start:start + self.kmer_size]
            if 'N' not in next_k:
                kmer_count += 1
                freq_sum += self.counts_dict[next_k]
        if raw_data:
            return freq_sum, kmer_count, len(seq)  # num_nucs
        else:
            return freq_sum / self.gnomad_chroms * self.test_chroms  # / len(seq)


def count_singletons(vcf_region):
    count = 0
    for v in vcf_region:
        if ek.is_quality_singleton(v):
            count += 1
    return count


def query_region(vcf_path, fasta, chrom, kmer_size, bins=100):
    # TODO: Add window size!!!
    vcf = VCF(vcf_path)
    fasta = Fasta(fasta)
    regions = ek.get_split_chrom_vcf(vcf_path, chrom, bins)
    window = KmerWindow(kmer_size)

    # def count_region(r):
    #     expected = (window.calculate_expected(str(fasta.get_seq(r.chrom, r.start, r.stop))))
    #     actual = (count_singletons(vcf(str(r))))
    #     print("%s\t%d\t%f" % (str(r), actual, expected), flush=True)

    # pool = mp.ProcessPool()
    print('region\tactual\texpected\tratio')
    # pool.map(count_region, regions)
    # pool.close()
    expected, actual = [], []
    for r in regions:
        exp = window.calculate_expected(str(fasta.get_seq(r.chrom, r.start, r.stop)))
        act = count_singletons(vcf(str(r)))
        expected.append(exp)
        actual.append(act)
        if exp == 0:
            ratio = 0
        else:
            ratio = act / exp
        print("%s\t%d\t%f\t%f" % (str(r), act, exp, ratio))
    # for a, e in zip(actual, expected):
    #     print('Actual: %d\tExpected: %f' % (a, e))


if __name__ == "__main__":
    fp = ek.get_test_path()
    vc = VCF(fp[0])
    f = Fasta(fp[1])
    query_region(fp[0], fp[1], 'chr22', 3, bins=50)
    # testseq = str(f.get_seq('chr1', 1_000_000, 1_500_000)).upper()
    # # repeat_as = "A" * 500000
    # windows = [KmerWindow(i) for i in [3, 5, 7]]
    # # for w in windows:
    # #     print("%d-mer: %f" % (w.kmer_size, w.calculate_expected(testseq)))
    # num_actual = count_singletons(vcf('chr1:1000000-1500000'))
    # for w in windows:
    #     print("%d-mer: %f\t%d" % (w.kmer_size, w.calculate_expected(testseq), num_actual))
