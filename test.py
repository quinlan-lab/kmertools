import eskedit as ek
import sys
import pandas as pd
import argparse


def test_train_kmer_model(arguments):
    if len(arguments) == 5:
        kmer_size = int(arguments[0])
        bedpath = arguments[1]
        vcfpath = arguments[2]
        fastapath = arguments[3]
        numprocs = int(arguments[4])
        invert_bed_selection = False
    elif len(arguments) == 6:
        kmer_size = int(arguments[0])
        bedpath = arguments[1]
        vcfpath = arguments[2]
        fastapath = arguments[3]
        numprocs = int(arguments[4])
        invert_bed_selection = bool(arguments[5])
    elif len(arguments) == 1:
        kmer_size = 3
        # big bed
        #bedpath = '/Users/simonelongo/too_big_for_icloud/merged_exons_grch38.bed'
        # small bed
        bedpath = '/Users/simonelongo/too_big_for_icloud/small_test.bed'
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        numprocs = 1
        invert_bed_selection = False
    else:
        return
    #     print("""Please enter arguments in this order:
    #                 1. kmer size
    #                 2. path to bed file
    #                 3. path to vcf file
    #                 4. path to reference fasta file
    #                 5. number of processors to use
    #                 6. invert bed selection (True/False)
    #             """)
    #     exit(1)
    outfile1 = 'TEST2regional_' + str(kmer_size) + 'mer_count.csv'
    outfile2 = 'TEST2regional_transitions_' + str(kmer_size) + 'mer.csv'
    result = ek.train_kmer_model(bedpath, vcfpath, fastapath, kmer_size, nprocs=numprocs, clean_bed=True,
                                 invert_selection=invert_bed_selection)
    print(result)
    pd.DataFrame.from_dict(result[0], orient='index').to_csv(outfile1)
    pd.DataFrame.from_dict(result[1], orient='index', columns=list('ACGT')).to_csv(outfile2)


def test_check_bed_regions_for_expected_mutations(arguments):
    if len(arguments) == 1:
        kmer_size = 7
        bedpath = '/Users/simonelongo/too_big_for_icloud/merged_exons_grch38.bed'
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        numprocs = 1
        countspath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/input_data/counts_data/7mer_relative_freq_noncoding.csv'
    elif len(arguments) == 6:
        kmer_size = int(arguments[0])
        bedpath = arguments[1]
        vcfpath = arguments[2]
        fastapath = arguments[3]
        countspath = arguments[4]
        numprocs = int(arguments[5])
    else:
        print("""Please enter arguments in this order:
            1. kmer size
            2. path to bed file
            3. path to vcf file
            4. path to reference fasta file
            5. path to counts file
            6. number of processors to use        
        """)
        exit(1)
    print('Querying regions for %d-mers using %d processsors\n' % (kmer_size, numprocs))
    ek.check_bed_regions(bedpath, vcfpath, fastapath, kmer_size, numprocs, counts_path=countspath)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', action='store_true', help='Query regions specified in bedfile for an expected number of mutations based on provided counts data.')
    parser.add_argument('--train', action='store_true', help='Build a counts table based on a k-mer model')
    parser.add_argument('--loctest', action='store_true')
    parser.add_argument('--kmer_size', '-k', action='store', dest='kmer_size', help='Length of k-mer motif')
    parser.add_argument('--bedpath', '-b', action='store', dest='bed_path',
                        help='Path to bed file containing genomic regions')
    parser.add_argument('--vcfpath', '-v', action='store', dest='vcf_path', help='Path to vcf file containing variants')
    parser.add_argument('--fastapath', '-f', action='store', dest='fasta_path', help='Path to reference genome fasta')
    parser.add_argument('--countspath', '-c', action='store', dest='countspath', default=None,
                        help='Path to counts table')
    parser.add_argument('--nprocs', '-N', action='store', dest='nprocs', default=1,
                        help='Number of processes to use (default=1)')
    parser.add_argument('--invert', action='store_true', help='If flag is present, will invert regions given in bed file.')

    args = parser.parse_args()
    if not args.query and not args.train:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if args.loctest:
        ksize = 7
        nprocs = 12
    else:
        ksize = int(args.kmer_size)
        nprocs = int(args.nprocs)
    bpath = args.bed_path
    vpath = args.vcf_path
    fpath = args.fasta_path
    cpath = args.countspath
    if args.query:
        if args.loctest:
            test_check_bed_regions_for_expected_mutations(['loctest'])
        else:
            test_check_bed_regions_for_expected_mutations([ksize, bpath, vpath, fpath, cpath, nprocs])
    if args.train:
        if args.loctest:
            test_train_kmer_model(['loctest'])
        else:
            test_train_kmer_model([ksize, bpath, vpath, fpath, nprocs, args.invert])
    # try:
    #     if args.loctest:
    #         ksize = 7
    #         nprocs = 12
    #     else:
    #         ksize = int(args.kmer_size)
    #         nprocs = int(args.nprocs)
    #     bpath = args.bed_path
    #     vpath = args.vcf_path
    #     fpath = args.fasta_path
    #     cpath = args.countspath
    #     if args.query:
    #         if args.loctest:
    #             test_check_bed_regions_for_expected_mutations(['loctest'])
    #         else:
    #             test_check_bed_regions_for_expected_mutations([ksize, bpath, vpath, fpath, cpath, nprocs])
    #     if args.train:
    #         if args.loctest:
    #             test_train_kmer_model(['loctest'])
    #         else:
    #             test_train_kmer_model([ksize, bpath, vpath, fpath, nprocs, args.invert])
    # except IndexError:
    #     print("Invalid parameters. Exiting...")
    #     parser.print_help()
    #     exit(1)


"""
various filepaths that may be useful:
-----------------------------------------------------
vcf_path = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
fasta_path = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
vcf_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/data/samp_build38.vcf.bgz'

vcf_path = '/scratch/general/lustre/u0319040/gnomadv3/gnomad.genomes.r3.0.sites.vcf.bgz'
fasta_path = '/scratch/general/lustre/u0319040/ref_genome/hg38/hg38.fa'
"""
