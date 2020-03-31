import datetime
import re
import sys
import time
from collections import defaultdict
import multiprocessing as mp
from cyvcf2 import VCF
from pyfaidx import FetchError, Fasta

from eskedit import GRegion, KmerWindow, is_singleton_snv, is_quality_snv, QueryWindow


# def count_regional_variants(vcf_region):
#     singletons = 0
#     count = 0
#     for v in vcf_region:
#         if is_singleton_snv(v):
#             singletons += 1
#         if is_quality_snv(v):
#             count += 1
#     return count, singletons


def count_regional_alleles(vcf_region):
    AF = 0.0
    AN = 0
    AC = 0
    singletons = 0
    count = 0
    for v in vcf_region:
        if is_quality_snv(v):
            try:
                newAC = int(v.INFO.get('AC'))
                AF += float(v.INFO.get('AF'))
                AC += newAC
                AN += int(v.INFO.get('AN'))
                count += 1
                if newAC == 1:
                    singletons += 1
            except ValueError:
                AF = 0
                AC = 0
                AN = 0
                print("ERROR: casting AC, AF, or AN in %s" % (str(vcf_region)), file=sys.stderr, flush=True)
    return AF, AC, AN, singletons, count


def query_bed_region(region, vcf_path, fasta, kmer_size, singleton_path, af_path, an_path, ac_path):
    """
    @param ac_path:
    @param an_path:
    @param af_path:
    @param singleton_path:
    @param region:
    @param vcf_path:
    @param fasta:
    @param kmer_size:
    @return:
    """
    # TODO: Add binning somehow (either keep equal size or equal number of bins
    start = time.time()
    vcf = VCF(vcf_path)
    fasta = Fasta(fasta)
    window = QueryWindow(kmer_size, singleton_path=singleton_path, af_path=af_path, an_path=an_path, ac_path=ac_path)
    # The first kmer actually begins centered around first nucleotide in sequence so
    # start position is shifted upstream by half the kmer length
    # end position is shifted downstream by the same
    shift = kmer_size // 2
    try:
        if region.strand is not None:
            if is_dash(region.strand):
                sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).complement.seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq.upper()
        else:
            sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq.upper()
        exp = window.calculate_expected(sequence)  # this does account for strandedness
        AF, AC, AN, singletons, count = count_regional_alleles(vcf(str(region)))
        field1 = count  # 'NumSNVs'
        field2 = singletons  # 'Singletons'
        field3 = AC  # 'AC'
        field4 = AN  # 'AN'
        field5 = AF  # 'AF'
        field6 = exp.get('singleton')  # 'ExpectedSingletons'
        field7 = exp.get('AC')  # 'ExpectedAC'
        field8 = exp.get('AN')  # 'ExpectedAN'
        field9 = exp.get('AF')  # 'ExpectedAF'

    except (KeyError, FetchError):
        field1 = 0  # 'NumSNVs'
        field2 = 0  # 'Singletons'
        field3 = 0  # 'AC'
        field4 = 0  # 'AN'
        field5 = 0  # 'AF'
        field6 = 0  # 'ExpectedSingletons'
        field7 = 0  # 'ExpectedAC'
        field8 = 0  # 'ExpectedAN'
        field9 = 0  # 'ExpectedAF'

    # print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
    #     region.printstr(), str(field1), str(field2), str(field3), str(field4), str(field5), str(field6), str(field7),
    #     str(field8), str(field9)), flush=True)
    regname = region.str_name().split('\t')
    print(
        '{: <8} {: <12} {: <12} {: <20} {: <8} {: <10} {: <12} {: <10} {: <10} {: <24} {: <22} {: <20} {: <20} {: <20}'.format(
            str(regname[0]), str(regname[1]), str(regname[2]), str(regname[3]), str(regname[4]), str(field1), str(field2), str(field3), str(field4), str(field5), str(field6), str(field7), str(field8), str(field9)),
        flush=True)
    return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        region.str_name(), str(field1), str(field2), str(field3), str(field4), str(field5), str(field6), str(field7),
        str(field8), str(field9))


def check_bed_regions(bed_path, vcf_path, fasta_path, kmer_size, nprocs=4, singleton_path=None, af_path=None,
                      an_path=None, ac_path=None, outfile=None,
                      strand_col=None, bed_names_col=None):
    field1 = 'NumSNVs'
    field2 = 'Singletons'
    field3 = 'AC'
    field4 = 'AN'
    field5 = 'AF'
    field6 = 'ExpectedSingletons'
    field7 = 'ExpectedAC'
    field8 = 'ExpectedAN'
    field9 = 'ExpectedAF'

    try:
        kmer_size = int(kmer_size)
        nprocs = int(nprocs)
    except ValueError:
        print('ERROR: kmer_size and nprocs must be integers', file=sys.stderr, flush=True)
        exit(1)
    additional_fields = defaultdict(int)
    if strand_col is not None:
        try:
            strand_col = int(strand_col)
            additional_fields['strand'] = strand_col
        except ValueError:
            strand_col = None
    if bed_names_col is not None:
        try:
            bed_names_col = int(bed_names_col)
            additional_fields['name'] = bed_names_col
        except ValueError:
            bed_names_col = None
    regions = []
    if outfile is None:
        outfile = bed_path.split('/')[-1]
        outfile = "".join(outfile.split('.')[:-1])
        outfile = '%s_%s.%dmer.dat' % (outfile, str(datetime.date.today()), kmer_size)

    with open(bed_path, 'r') as bedfile:
        kwargs = defaultdict(str)
        for line in bedfile.readlines():
            fields = line.split('\t')
            for k, v in additional_fields.items():
                kwargs[k] = fields[v]
            regions.append(GRegion(*[fields[0], fields[1], fields[2]], **kwargs))

    if (bed_names_col is not None) and (strand_col is not None):
        name_header = regions[0].name_header().split('\t')
    else:
        name_header = 'Region\t \t \t \t '.split('\t')

    with open(outfile, 'w') as output:
        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            "".join(name_header), field1, field2, field3, field4, field5, field6, field7, field8))

    # print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
    #     name_header, field1, field2, field3, field4, field5, field6, field7, field8, field9), flush=True)
    print(
        '{: <8} {: <12} {: <12} {: <20} {: <8} {: <10} {: <12} {: <10} {: <10} {: <24} {: <22} {: <20} {: <20} {: <20}'.format(
            *name_header, field1, field2, field3, field4, field5, field6, field7, field8, field9), flush=True)
    arguments = [(region, vcf_path, fasta_path, kmer_size, singleton_path, af_path, an_path, ac_path) for region in
                 regions]
    pool = mp.Pool(nprocs)
    results = pool.starmap_async(query_bed_region, arguments)
    pool.close()
    pool.join()
    with open(outfile, 'w') as output:
        for result in results.get():
            output.write(result)
    return


def is_dash(pdash):
    regex = '[\u002D\u058A\u05BE\u1400\u1806\u2010-\u2015\u2E17\u2E1A\u2E3A\u2E3B\u2E40\u301C\u3030\u30A0\uFE31\uFE32\uFE58\uFE63\uFF0D]'
    if re.match(regex, pdash) is not None:
        return True
    return False


if __name__ == "__main__":
    # def check_bed_regions(bed_path, vcf_path, fasta_path, kmer_size, nprocs=4, singleton_path=None, af_path=None,
    #                       an_path=None, ac_path=None, outfile=None,
    #                       strand_col=None, bed_names_col=None):
    # small bed
    bedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/kmertools/QUERY_TEST_1klines.bed'
    vcfpath = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
    fastapath = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'
    check_bed_regions(bedpath, vcfpath, fastapath, 7, nprocs=12, strand_col=5, bed_names_col=3,
                      outfile='QUERY_TEST_RUN1.dat')
