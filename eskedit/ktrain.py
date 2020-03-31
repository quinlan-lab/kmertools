import array
import time
from collections import defaultdict
import multiprocessing as mp
import pandas as pd
import sys
from cyvcf2 import VCF
from pyfaidx import FetchError, Fasta
from eskedit import get_bed_regions, ModelContainer, is_dash, Variant, complete_sequence, is_quality_snv, kmer_search
from signal import signal, SIGINT


# def sigint_handler(signal_received, frame):
#     # print('SIGINT or CTRL-C detected. Exiting gracefully')
#     exit(0)


def model_region(datacontainer, vcf_path, fasta_path, kmer_size, region, AC_cutoff):
    if AC_cutoff is not None:
        try:
            AC_cutoff = int(AC_cutoff)
        except ValueError:
            AC_cutoff = None
            print('AC cutoff must be a positive integer. Ignoring user value and using SNVs with any AC.',
                  file=sys.stderr, flush=True)
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
            if is_dash(region.strand):
                sequence = fasta.get_seq(region.chrom, region.start - kmer_mid_idx,
                                         region.stop + kmer_mid_idx).complement.seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start - kmer_mid_idx,
                                         region.stop + kmer_mid_idx).seq.upper()
        else:
            sequence = fasta.get_seq(region.chrom, region.start - kmer_mid_idx, region.stop + kmer_mid_idx).seq.upper()
    except (KeyError, FetchError):
        print('Region %s not found in fasta, continuing...' % str(region), file=sys.stderr, flush=True)
        return
    region_ref_counts = kmer_search(sequence, kmer_size)  # nprocs=1 due to short region
    r_string = str(region.chrom) + ':' + str(region.start) + '-' + str(region.stop)
    singleton_transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    ac_transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    an_transitions = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    af_transitions = defaultdict(lambda: array.array('d', [0, 0, 0, 0]))
    # Define indices for nucleotides
    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    idx_nuc = list('ACGT')
    for variant in vcf(r_string):
        if is_quality_snv(variant, AC_cutoff=AC_cutoff):
            adj_seq = fasta[str(variant.CHROM)][(variant.POS - start_idx_offset):(variant.POS + kmer_mid_idx)].seq
            if str(adj_seq[kmer_mid_idx]).upper() != str(variant.REF).upper():
                print('WARNING: Reference mismatch\tFasta REF: %s\tVCF REF: %s' % (adj_seq[kmer_mid_idx], variant.REF),
                      file=sys.stderr, flush=True)
            if complete_sequence(adj_seq):
                ac_transitions[adj_seq.upper()][nuc_idx[variant.ALT[0]]] += variant.INFO.get('AC')
                an_transitions[adj_seq.upper()][nuc_idx[variant.ALT[0]]] += variant.INFO.get('AN')
                af_transitions[adj_seq.upper()][nuc_idx[variant.ALT[0]]] += variant.INFO.get('AF')
                if variant.INFO.get('AC') == 1:
                    singleton_transitions[adj_seq.upper()][nuc_idx[variant.ALT[0]]] += 1
    data = {'singleton': singleton_transitions,
            'AC': ac_transitions,
            'AN': an_transitions,
            'AF': af_transitions}
    temp = datacontainer.get()
    temp.add_kmer_counts(region_ref_counts)
    for k, v in data.items():
        temp.add_transition(v, k)

    datacontainer.set(temp)
    print('Finished region %s in %s' % (str(region), str(time.time() - start)), flush=True)
    return


def train_kmer_model(bed_path, vcf_path, fasta_path, kmer_size, nprocs=1,
                     header=False, strand_col=None, bed_names_col=None, AC_cutoff=None):
    """
    Builds the counts tables required for the k-mer model. Returned as 2 dictionaries.
    @param AC_cutoff:           Specify to filter out variants above a given AC (AC > cutoff will be filtered)
                                ** works only if keyword 'nonsingletons=True' **
    @param strand_col:          zero-based column index of strand information from bed file
    @param bed_names_col:       zero-based column index of name information from bed file
    @param bed_path:            path to bed file
    @param vcf_path:            path to vcf file
    @param fasta_path:          path to reference fasta
    @param kmer_size:           NOTE unpredictable behavior may occur if even numbers are used here.
    @param nprocs:              number of processors to use
    @param header:              False (default) if the bed file does not have a header
    @return:
    """
    start_time = time.time()

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
    dc = manager.Value(ModelContainer, ModelContainer())
    pool = mp.Pool(nprocs)

    def sigint_handler(signal_received, frame):
        model = dc.get()
        model.writetofile()
        pool.terminate()
        pool.join()
        exit(0)

    signal(SIGINT, sigint_handler)

    regions = get_bed_regions(bed_path, invert_selection=False, header=header, clean_bed=True,
                              strand_col=strand_col, bed_names_col=bed_names_col)
    # Bundle arguments to pass to 'model_region' function
    arguments = [(dc, vcf_path, fasta_path, kmer_size, region, AC_cutoff) for region in regions]
    # Distribute workload
    pool.starmap(model_region, arguments)
    pool.close()
    pool.join()

    print("Done in {}".format(time.time() - start_time))

    dc.get().writetofile()

    return  # master_ref_counts, transitions_list


if __name__ == "__main__":
    bedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/QUERY_TEST_100lines.bed'
    vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
    fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
    train_kmer_model(bedpath, vcfpath, fastapath, 3, nprocs=12)
