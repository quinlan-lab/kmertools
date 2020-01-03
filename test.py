import eskedit as ek
import sys
import pandas as pd

if __name__ == "__main__":
    kmer_size = int(sys.argv[1])
    bedpath = sys.argv[2]
    vcfpath = sys.argv[3]
    fastapath = sys.argv[4]
    results = ek.train_model(bedpath, vcfpath, fastapath, invert_selection=True, clean_bed=True)
    outfile1 = 'regional_transitions_' + str(kmer_size) + 'mer.csv'
    outfile2 = 'regional_' + str(kmer_size) + 'mer_count.csv'
    pd.DataFrame.from_dict(results['transitions'], orient='index').to_csv(outfile1)
    pd.DataFrame.from_dict(results['ref_count'], orient='index').to_csv(outfile2)

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
