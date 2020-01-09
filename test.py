import eskedit as ek
import sys
import pandas as pd
import argparse


def test_train_kmer_model(args):
    kmer_size = int(args[1])
    bedpath = args[2]
    vcfpath = args[3]
    fastapath = args[4]
    nprocs = int(args[5])
    outfile1 = 'regional_' + str(kmer_size) + 'mer_count.csv'
    outfile2 = 'regional_transitions_' + str(kmer_size) + 'mer.csv'
    result = ek.train_kmer_model(bedpath, vcfpath, fastapath, kmer_size, nprocs=nprocs, clean_bed=True,
                                 invert_selection=False)
    pd.DataFrame.from_dict(result[0], orient='index').to_csv(outfile1)
    pd.DataFrame.from_dict(result[1], orient='index').to_csv(outfile2)


def test_check_bed_regions_for_expected_mutations(args):
    if len(args) == 1:
        kmer_size = 7
        bedpath = '/Users/simonelongo/too_big_for_icloud/merged_exons_grch38.bed'
        vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
        fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
        numprocs = 12
        countspath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/input_data/counts_data/7mer_relative_freq_noncoding.csv'
    elif len(args) == 6:
        kmer_size = int(args[0])
        bedpath = args[1]
        vcfpath = args[2]
        fastapath = args[3]
        countspath = args[4]
        numprocs = int(args[5])
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
    parser.add_argument('--query', action='store_true')
    parser.add_argument('train', action='store_true')
    parser.add_argument('--kmer_size', '-k', action='store', dest='kmer_size', help='Length of k-mer motif',
                        required=True)
    parser.add_argument('--bedpath', '-b', action='store', dest='bed_path',
                        help='Path to bed file containing genomic regions', required=True)
    parser.add_argument('--vcfpath', '-v', action='store', dest='vcf_path', help='Path to vcf file containing variants',
                        required=True)
    parser.add_argument('--fastapath', '-f', action='store', dest='fasta_path', help='Path to reference genome fasta',
                        required=True)
    parser.add_argument('--countspath', '-c', action='store', dest='countspath', default=None,
                        help='Path to counts table')
    parser.add_argument('--nprocs', '-N', action='store', dest='nprocs', default=1,
                        help='Number of processes to use (default=1)')
    if len(sys.argv) < 6:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    try:
        ksize = int(args.kmer_size)
        bpath = args.bed_path
        vpath = args.vcf_path
        fpath = args.fasta_path
        cpath = args.countspath
        nprocs = int(args.nprocs)
        if args.query:
            test_check_bed_regions_for_expected_mutations([ksize, bpath, vpath, fpath, cpath, nprocs])
    except ValueError:
        print("Invalid parameters. Exiting...")
        parser.print_help()
        exit(1)
    # test_train_kmer_model(sys.argv)
    test_check_bed_regions_for_expected_mutations(sys.argv)

    # results = ek.train_model_old(bedpath, vcfpath, fastapath, kmer_size, invert_selection=True, clean_bed=True)

    # pd.DataFrame.from_dict(results['transitions'], orient='index').to_csv(outfile1)
    # pd.DataFrame.from_dict(results['ref_count'], orient='index').to_csv(outfile2)

# # vcf_path = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
# fasta_path = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
# vcf_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/data/samp_build38.vcf.bgz'
#
# # vcf_path = '/scratch/general/lustre/u0319040/gnomadv3/gnomad.genomes.r3.0.sites.vcf.bgz'
# # fasta_path = '/scratch/general/lustre/u0319040/ref_genome/hg38/hg38.fa'
#
# try:
#     ksize = int(sys.argv[1])
# except (ValueError, IndexError):
#     ksize = 3
#
# try:
#     n_threads = int(sys.argv[2])
# except (ValueError, IndexError):
#     n_threads = 6
#
# if '-save-vars' in sys.argv:
#     var_out = 'singletons_test.bed'
# else:
#     var_out = None
#
# t_fpath = 'counts_' + str(ksize) + 'mer.bed'
# results = ek.get_kmer_context(vcf_path, fasta_path, ksize, nprocs=n_threads)
# # results should be  list of dictionaries
# transitions = []
# positions = []
# for coll in results:
#     for k, v in coll:
#         if 'position' in k:
#             positions.append(v)
#         else:
#             transitions.append(v)
#
# # s_transitions = defaultdict(Counter)
# # s_positions = defaultdict(Variant)
# ek.merge_positions_dd(positions, outfile=var_out)
# ek.merge_transitions_ddc(transitions, outfile=t_fpath)

# /scratch/general/lustre/u0319040/gnomadv3/gnomad.genomes.r3.0.sites.vcf.bgz
# /scratch/general/lustre/u0319040/ref_genome/hg38/hg38.fa

# /uufs/chpc.utah.edu/common/home/u0319040/kmertools/test.py
