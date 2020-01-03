"""
The goal of this module is to take a BED file containing regions of the genome we wish to exclude.
This module will then examine those regions and build a kmer mutability model based on those regions.
"""
from collections import defaultdict, Counter
import time
from cyvcf2 import VCF
from pyfaidx import Fasta
import eskedit as ek
# import multiprocessing as mp
from eskedit import VCFRegion, get_autosome_lengths_grch38, get_grch38_chroms, RegionContainer, Variant


# def gis(bpath, header):
#     regions = []
#     with open(bpath, 'r') as bedfile:
#         if header:
#             bedfile.readline()
#         lines = [line.split('\t') for line in bedfile.readlines()]
#         prev_chrom = ''
#         prev_start = 1
#         prev_stop = -1
#         chrom_names = get_grch38_chroms()
#         new_region = VCFRegion(None, None, None)
#         for i, line in enumerate(lines):
#             current_chrom = line[0]
#             current_start = int(line[1])
#             current_end = int(line[2])
#             if current_chrom != prev_chrom:
#                 if i > 0 and prev_stop < chrom_names[prev_chrom]:
#                     # add remaining old chromosome and start new
#                     regions.append(VCFRegion(prev_chrom, prev_stop, chrom_names[prev_chrom]))
#                     new_region.chrom = current_chrom
#                 if current_start > 1:
#                     regions.append(VCFRegion(current_chrom, 1, current_start))
#                     new_region.start = current_end
#                     new_region.chrom = current_chrom
#                     prev_stop = current_end
#                     prev_start = current_start
#                     prev_chrom = current_chrom
#                 else:
#                     new_region.chrom = current_chrom
#                     new_region.start = current_end
#                     prev_stop = current_end
#                     prev_start = current_start
#                     prev_chrom = current_chrom
#             else:  # if we are on same chromosome
#                 if new_region.start is None:
#                     new_region.start = current_end
#                 else:  # if new_region.stop is None
#                     new_region.stop = current_start
#                 prev_stop = current_end
#                 prev_start = current_start
#                 prev_chrom = current_chrom
#                 # if current_start <= prev_stop:
#                 #     if current_end > prev_stop:
#                 #         new_region.stop = current_end
#             # if new_region.start is None:
#             #     new_region.start = current_end
#             # else:
#             #     if current_start > new_region.start:
#             #         new_region.stop = current_start
#             #     prev_start = current_start
#             #     prev_stop = current_end
#             #     prev_chrom = current_chrom
#             if new_region.is_complete() and regions[-1] != new_region:
#                 regions.append(new_region)
#                 new_region = VCFRegion(current_chrom, None, None)
#     return regions


def get_inverse_selection(bed_path, header=False):
    regions = []
    chrom_info = get_autosome_lengths_grch38()
    prev_chrom = ''
    with open(bed_path, 'r') as bedfile:
        if header:
            bedfile.readline()  # skip header
        lines = [line.split('\t') for line in bedfile.readlines()]
        for i, line in enumerate(lines):
            if i == 0:
                prev_chrom = lines[i][0]
                if int(lines[i][1]) > 1:
                    regions.append(VCFRegion(prev_chrom, 1, lines[i][1]))
            else:
                if lines[i][0] != prev_chrom:
                    if int(lines[i - 1][2]) < chrom_info[prev_chrom]:
                        regions.append(VCFRegion(prev_chrom, lines[i - 1][2], chrom_info[prev_chrom]))
                    if int(lines[i][1]) > 1:
                        regions.append(VCFRegion(lines[i][0], 1, lines[i][1]))
                    prev_chrom = lines[i][0]
                else:
                    regions.append(VCFRegion(prev_chrom, lines[i - 1][2], lines[i][1]))
            if i == len(lines) - 1:
                if int(lines[i][2]) < chrom_info[prev_chrom]:
                    regions.append(VCFRegion(prev_chrom, lines[i][2], chrom_info[prev_chrom]))

    return regions


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


def train_model(bed_path, vcf_path, fasta_path, kmer_size, nprocs=None, invert_selection=False, clean_bed=False,
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
