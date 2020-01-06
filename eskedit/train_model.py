"""
The goal of this module is to take a BED file containing regions of the genome we wish to exclude.
This module will then examine those regions and build a kmer mutability model based on those regions.
"""
from collections import defaultdict, Counter
import time
from cyvcf2 import VCF
from pyfaidx import Fasta
import eskedit as ek
import multiprocessing as mp
from eskedit import VCFRegion, get_autosome_lengths_grch38, get_grch38_chroms, RegionContainer, Variant


def get_bed_regions(bed_path, invert_selection=False, header=False, clean_bed=False):
    regions = RegionContainer()
    with open(bed_path, 'r') as bedfile:
        if header:
            bedfile.readline()
        for line in bedfile.readlines():
            fields = line.split('\t')
            if clean_bed:
                regions.add_distinct_region(VCFRegion(fields[0], fields[1], fields[2]))
            else:
                regions.add_region(VCFRegion(fields[0], fields[1], fields[2]))
    if invert_selection:
        return regions.get_inverse()
    else:
        return regions.get_regions()


def model_region(vcf_path, fasta_path, kmer_size, region):
    start = time.time()
    fasta = Fasta(fasta_path)
    vcf = VCF(vcf_path)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    sequence = fasta.get_seq(region[0], region[1], region[2]).seq.upper()
    region_ref_counts = ek.kmer_search(sequence, kmer_size)  # nprocs=1 due to short region
    r_string = str(region[0]) + ':' + str(region[1]) + '-' + str(region[2])
    transitions = defaultdict(Counter)
    for variant in vcf(r_string):
        if ek.is_quality_singleton(variant):
            new_var = Variant(variant=variant, fields=['vep'])
            # TODO: with less memory overhead variant_positions[new_var.INDEX] = new_var
            # take 7mer around variant. pyfaidx excludes start index and includes end index
            adj_seq = fasta[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
            if ek.complete_sequence(adj_seq):
                transitions[adj_seq.upper()][new_var.ALT[0]] += 1
    print('Finished region %s in %s' % (str(region), str(time.time() - start)))
    return region_ref_counts, transitions


def train_model(bed_path, vcf_path, fasta_path, kmer_size, nprocs=None, invert_selection=False, clean_bed=False,
                header=False):
    regions = get_bed_regions(bed_path, invert_selection=invert_selection, header=header, clean_bed=clean_bed)
    arguments = [(vcf_path, fasta_path, kmer_size, region.flist) for region in regions]
    pool = mp.Pool(nprocs)
    results = pool.starmap_async(model_region, arguments)
    pool.close()
    pool.join()
    master_ref_counts = Counter()
    transitions_list = []
    for result in results.get():
        for k, v in result[0].items():
            master_ref_counts[k] += v
        transitions_list.append(result[1])
    merged_transitions = ek.merge_transitions_ddc(transitions_list, outfile=None)
    return master_ref_counts, merged_transitions


def train_model_old(bed_path, vcf_path, fasta_path, kmer_size, nprocs=None, invert_selection=False, clean_bed=False,
                    header=False):
    regions = get_bed_regions(bed_path, invert_selection=invert_selection, header=header, clean_bed=clean_bed)
    fasta = Fasta(fasta_path)
    vcf = VCF(vcf_path)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)
    total_counts = []
    variant_positions = defaultdict(Variant)
    variant_transitions = defaultdict(Counter)
    for region in regions:
        start = time.time()
        seq = str(fasta.get_seq(region.chrom, region.start, region.stop))
        total_counts.append(ek.get_kmer_count(seq, kmer_size, nprocs=1))
        for variant in vcf(str(region)):
            if ek.is_quality_singleton(variant):
                new_var = Variant(variant=variant, fields=['vep'])
                # TODO: with less memory overhead variant_positions[new_var.INDEX] = new_var
                # take 7mer around variant. pyfaidx excludes start index and includes end index
                adj_seq = fasta[str(new_var.CHROM)][(new_var.POS - start_idx_offset):(new_var.POS + kmer_mid_idx)].seq
                if ek.complete_sequence(adj_seq):
                    variant_transitions[adj_seq.upper()][new_var.ALT[0]] += 1
        print('Finished region %s in %s' % (str(region), str(time.time() - start)))
    master_count = Counter()
    for count in total_counts:
        for k, v in count.items():
            master_count[k] += v
    return {'variants': variant_positions, 'transitions': variant_transitions, 'ref_count': master_count}

# if __name__ == "__main__":
#     # test_path = '/Users/simonelongo/too_big_for_icloud/exons_grch38.bed'
#     # test_path = '/Users/simonelongo/too_big_for_icloud/exons_chr1_grch38.bed'
#     test_path = '/Users/simonelongo/too_big_for_icloud/chr1_test.bed'
#     # test_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/input_data/test_data/test.bed'
#     test_reg = get_bed_regions(test_path, invert_selection=True)
#     for tregion in test_reg:
#         print(tregion)
