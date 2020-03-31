import eskedit as ek
import sys
import pandas as pd
import argparse
import time


# from signal import signal, SIGINT
#
#
# def sigint_handler(signal_received, frame):
#     # print('SIGINT or CTRL-C detected. Exiting gracefully')
#     exit(0)


def test_train_kmer_model(arguments, test=None):
    if test is None:
        result = ek.train_kmer_model(arguments.bed_path, arguments.vcf_path, arguments.fasta_path, arguments.kmer_size,
                                     nprocs=arguments.nprocs,
                                     clean_bed=True,
                                     invert_selection=arguments.invert, strand_col=arguments.strand_col,
                                     bed_names_col=arguments.bed_name_col, singletons=arguments.check_singletons,
                                     nonsingletons=arguments.check_AF)
        for k, v, in result.items():
            outfile1 = 'regional_' + str(k) + '_' + str(arguments.kmer_size) + 'mer_count.csv'
            outfile2 = 'regional_transitions_' + str(k) + '_' + str(arguments.kmer_size) + 'mer.csv'
            if len(v) == 3:
                outfile2 = 'regional_ACtransitions_' + str(k) + '_' + str(arguments.kmer_size) + 'mer.csv'
                outfile3 = 'regional_ANtransitions_' + str(k) + '_' + str(arguments.kmer_size) + 'mer.csv'
                pd.DataFrame.from_dict(v[2], orient='index', columns=list('ACGT')).to_csv(outfile3)
            pd.DataFrame.from_dict(v[0], orient='index').to_csv(outfile1)
            pd.DataFrame.from_dict(v[1], orient='index', columns=list('ACGT')).to_csv(outfile2)
    else:  # Test mode
        kmer_size = 3
        # big bed
        # bedpath = '/Users/simonelongo/too_big_for_icloud/merged_exons_grch38.bed'
        # small bed
        bedpath = '/Users/simonelongo/too_big_for_icloud/small_test.bed'
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        numprocs = 8
        result = ek.train_kmer_model(bedpath, vcfpath, fastapath, kmer_size, nprocs=numprocs, clean_bed=True,
                                     invert_selection=False, strand_col=None,
                                     bed_names_col=None, singletons=False, nonsingletons=True)
        for k, v in result.items():
            outfile1 = 'TEST_regional_' + str(k) + '_' + str(kmer_size) + 'mer_count.csv'
            outfile2 = 'TEST_regional_transitions_' + str(k) + '_' + str(kmer_size) + 'mer.csv'
            if len(v) == 3:
                outfile2 = 'TEST_regional_ACtransitions_' + str(k) + '_' + str(kmer_size) + 'mer.csv'
                outfile3 = 'TEST_regional_ANtransitions_' + str(k) + '_' + str(kmer_size) + 'mer.csv'
                pd.DataFrame.from_dict(v[2], orient='index', columns=list('ACGT')).to_csv(outfile3)
            pd.DataFrame.from_dict(v[0], orient='index').to_csv(outfile1)
            pd.DataFrame.from_dict(v[1], orient='index', columns=list('ACGT')).to_csv(outfile2)
    return


def test_check_bed_regions_for_expected_mutations(arguments, test=None):
    if test is None:
        ek.check_bed_regions(arguments.bed_path, arguments.vcf_path, arguments.fasta_path, arguments.kmer_size,
                             arguments.nprocs,
                             counts_path=arguments.countspath,
                             strand_col=arguments.strand_col, bed_names_col=arguments.bed_name_col,
                             singletons=arguments.check_singletons)
    else:
        kmer_size = 7
        # bedpath = '/Users/simonelongo/too_big_for_icloud/merged_exons_grch38.bed'
        # bedpath = '/Users/simonelongo/Downloads/3primeUTR.refseq.grch38.hg38.bed'
        bedpath = '/Users/simonelongo/Downloads/hg38-tRNAs/hg38-tRNAs.bed'
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        numprocs = 6
        countspath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/input_data/counts_data/7mer_relative_freq_noncoding.csv'
        strand = 5
        bednames = 3
        ek.check_bed_regions(bedpath, vcfpath, fastapath, kmer_size, numprocs, counts_path=countspath,
                             strand_col=strand, bed_names_col=bednames, singletons=False)

    return


def test_chrom_bin_mutability(arguments, test=None):  # vcfpath, fastapath, kmer_size, nbins, chroms=None, numprocs=1):
    if test is None and arguments.bed_path is None:
        mut_table = ek.chrom_bin_mutability(arguments.vcf_path, arguments.fasta_path, arguments.kmer_size,
                                            arguments.nbins,
                                            chroms=arguments.chrom_list, nprocs=int(arguments.nprocs),
                                            af=arguments.check_AF)
        mut_table.to_csv('chrom_%sbins_%smers.csv' % (arguments.nbins, arguments.kmer_size))
    elif arguments.bed_path is not None and test is None:
        mut_table = ek.region_mutability_from_bed(arguments.vcf_path, arguments.fasta_path, arguments.bed_path,
                                                  arguments.kmer_size,
                                                  nprocs=int(arguments.nprocs),
                                                  af=arguments.check_AF)
        # mut_table.to_csv('%s_%smers.csv' % ("".join(arguments.bed_path.split('.')[:-1]), arguments.kmer_size))
    # TEST Below
    else:
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        kmer_size = 3
        nbins = 24
        chroms = ['chr22']
        numprocs = 6
        AF = True
        bedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/BINCHROM_TEST.bed'
        # mut_table = ek.chrom_bin_mutability(vcfpath, fastapath, kmer_size, nbins, chroms=chroms, nprocs=numprocs, af=AF)
        outname = bedpath.split('/')[-1]
        mut_table = ek.region_mutability_from_bed(vcfpath, fastapath, bedpath,
                                                  kmer_size,
                                                  nprocs=numprocs,
                                                  af=AF)
        # mut_table.to_csv('%s_%smers.csv' % ("".join(outname.split('.')[:-1]), kmer_size))

        # mut_table.to_csv('TEST_chrom_%dbins_%dmers.csv' % (nbins, kmer_size))
    return


def test_ktrain(arguments, test=None):
    if test is None:
        # run regular
        ek.train_kmer_model(arguments.bed_path, arguments.vcf_path, arguments.fasta_path,
                            int(arguments.kmer_size), nprocs=arguments.nprocs,
                            strand_col=arguments.strand_col, bed_names_col=arguments.bed_name_col)
    else:
        # runtest
        bedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/QUERY_TEST_100lines.bed'
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        ksize = 7
        numprocs = 6
        ek.train_kmer_model(bedpath, vcfpath, fastapath, ksize, nprocs=numprocs, strand_col=5, bed_names_col=3)
    return


if __name__ == "__main__":
    # register SIGINT handler
    # signal(SIGINT, sigint_handler)
    start = time.time()

    # keywords for script functionality
    FUNCTION_MAP = {
        'query': test_check_bed_regions_for_expected_mutations,
        'train': test_ktrain,
        'binchrom': test_chrom_bin_mutability,
    }

    parser = argparse.ArgumentParser(prog='ESKeDiT',
                                     usage="Execute functions with subcommands:\n$ python3 test.py <COMMAND> -[ARGUMENTS]\n\n$ python3 test.py query -k 3 -b bedfile.bed -v variants.vcf -f hg38.fa -N 32\n\n")
    parser.add_argument('command', choices=FUNCTION_MAP.keys())
    parser.add_argument('--version', action='version', version='%s %s' % ('ESKeDiT', ek.VERSION))
    parser.add_argument('--test', action='store_true', dest='loctest', help=argparse.SUPPRESS)
    parser.add_argument('--kmer_size', '-k', action='store', dest='kmer_size', help='Length of k-mer motif')
    parser.add_argument('--bedpath', '-b', action='store', dest='bed_path',
                        help='Path to bed file containing genomic regions')
    parser.add_argument('--vcfpath', '-v', action='store', dest='vcf_path', help='Path to vcf file containing variants')
    parser.add_argument('--fastapath', '-f', action='store', dest='fasta_path', help='Path to reference genome fasta')
    parser.add_argument('--countspath', '-c', action='store', dest='countspath', default=None,
                        help='Path to counts table')
    parser.add_argument('--nbins', action='store', dest='nbins', help='Number of bins to split each chromosome into')
    parser.add_argument('--chrom_list', action='store', dest='chrom_list',
                        help='Comma-separated list of chromosomes to evaluate (eg. \'--chrom_list chr1,chr2,chr3\'). Default is all autosomes.')
    parser.add_argument('--nprocs', '-N', action='store', dest='nprocs', default=1,
                        help='Number of processes to use (default=1)')
    parser.add_argument('--invert', action='store_true',
                        help='If flag is present, will invert regions given in bed file.')
    parser.add_argument('--strand', '-S', nargs='?', const=5, type=int, dest='strand_col',
                        help='Enter (zero-based) integer value of column in bed file with strand information if it\'s not 5')
    parser.add_argument('--bed_names', nargs='?', const=3, type=int, dest='bed_name_col',
                        help='Enter (zero-based) integer value of column in bed file with region/gene name information')
    parser.add_argument('--singletons', action='store_true', dest='check_singletons')
    parser.add_argument('--AF', action='store_true', dest='check_AF')
    # add_argument('-t', '--test', nargs='?', const=5, type=int)
    args = parser.parse_args()

    if args.loctest:
        func = FUNCTION_MAP[args.command]
        func(args, test='test')
    else:
        try:
            args.kmer_size = int(args.kmer_size)
            args.nprocs = int(args.nprocs)
            if args.nprocs < 1 or args.kmer_size < 1:
                raise ValueError
        except ValueError:
            print('nprocs and kmer_size must be positive integers')
            exit(1)
        func = FUNCTION_MAP[args.command]
        func(args)
    print('\nFinished in %s.\n' % str(time.time() - start))
    exit(0)

"""
various filepaths that may be useful:
-----------------------------------------------------
vcf_path = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
fasta_path = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
vcf_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/data/samp_build38.vcf.bgz'

vcf_path = '/scratch/general/lustre/u0319040/gnomadv3/gnomad.genomes.r3.0.sites.vcf.bgz'
fasta_path = '/scratch/general/lustre/u0319040/ref_genome/hg38/hg38.fa'
"""
