import argparse
import time

if __name__ == "__main__":
    print('Running')
    starttime = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf-path', '-v', action='store', dest='vcf_path', help='Enter path to VCF file',
                        default=None)
    parser.add_argument('--fasta-ref', '-fr', action='store', dest='fasta_ref', default=None,
                        help='Enter path to reference fasta')

    # start = time.time()
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--destination', '-df', action='store', dest='destination_dir', default=os.getcwd())
    # parser.add_argument('--fullscan', '-fs', action='store', dest='full_scan', default=None)
    # parser.add_argument('--numthreads', '-N', action='store', dest='N_THREADS', help='Number of threads available',
    #                     default=4)
    # parser.add_argument('--ref-directory', '-rd', action='store', dest='ref_dir_path', default=None,
    #                     help='Enter path to directory containing indexed VCF and Fasta files.')
    # parser.add_argument('--vcf-path', '-v', action='store', dest='vcf_path',
    #                     help='Enter path to VCF file or CSV file containing variants', default=None)
    # parser.add_argument('--variant-csv', '-t', action='store', dest='csv_var_fpath',
    #                     help='Enter filepath for CSV file containing variants if available', default=None)
    # parser.add_argument('--reference-genome', '-rg', action='store', dest='ref_fasta',
    #                     help='Enter path to reference genome fasta file',
    #                     default='/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    # parser.add_argument('--kmer-len', '-k', action='store', dest='kmer_length', help='Enter length of Kmer', default=3)
    # parser.add_argument('--output', '-o', action='store', dest='output_file', help='Enter desired output file name',
    #                     default='variants_samp.csv')
    # if len(sys.argv) == 1:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)
    # args = parser.parse_args()
    # try:
    #     destination_directory = args.destination_dir
    #     is_wkflow = args.full_scan is not None
    #     nthreads = int(args.N_THREADS)
    #     vcf_path = args.vcf_path
    #     kmer_len = int(args.kmer_length)
    #     o_file = args.output_file
    #     csv_fpath = args.csv_var_fpath
    #     fasta_path = args.ref_fasta
    #     ref_dir = args.ref_dir_path
    # except ValueError:
    #     print("Invalid parameters. Exiting...")
    #     exit(1)
    # if is_wkflow:
    #     gv = GenVCF(fasta_path, vcf_path, kmer_len, nthreads)
    #     gv.set_destination(destination_directory)
    #     gv.vcf_scan()
    #     # vcf_process.vcf_singleton_analysis(vcf_path, fasta_path, kmer_len, nprocs=nthreads, var_out_file=o_file)
    #     print("Done in " + str(time.time() - start))
    #     exit(0)
