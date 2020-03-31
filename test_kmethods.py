from unittest import TestCase
from eskedit.kmethods import *


class Test(TestCase):

    def test_generate_kmers(self):
        print('Testing %s' % 'test_generate_kmers')
        for i in range(1, 8):
            self.assertEqual(len(generate_kmers(i)), 4 ** i)

    def test_gen_random_sequence(self):
        print('Testing %s' % 'test_gen_random_sequence')
        self.assertTrue(len(gen_random_sequence(7)) == 7)
        self.assertTrue(True)

    def test_ref_genome_as_string(self):
        print('Testing %s' % 'test_ref_genome_as_string')
        self.assertTrue(True)

    def test_complement(self):
        print("Testing %s" % "test_complement")
        self.assertTrue(True)

    def test_get_complementary_sequence(self):
        print("Testing %s" % "test_get_complementary_sequence")
        self.assertTrue(True)

    def test_is_quality_snv(self):
        print("Testing %s" % "test_is_quality_snv")
        self.assertTrue(True)

    def test_is_quality_nonsingleton(self):
        print("Testing %s" % "test_is_quality_nonsingleton")
        self.assertTrue(True)

    def test_is_quality_singleton(self):
        print("Testing %s" % "test_is_quality_singleton")
        self.assertTrue(True)

    def test_is_singleton_snv(self):
        print("Testing %s" % "test_is_singleton_snv")
        self.assertTrue(True)

    def test_complete_sequence(self):
        print("Testing %s" % "test_complete_sequence")
        self.assertTrue(True)

    def test_kmer_search(self):
        print("Testing %s" % "test_kmer_search")
        self.assertTrue(True)

    def test_get_vcf_info_fields(self):
        print("Testing %s" % "test_get_vcf_info_fields")
        self.assertTrue(True)

    def test_get_kmer_count(self):
        print("Testing %s" % "test_get_kmer_count")
        self.assertTrue(True)

    def test_merge_transitions_ddc(self):
        print("Testing %s" % "test_merge_transitions_ddc")
        self.assertTrue(True)

    def test_merge_positions_dd(self):
        print("Testing %s" % "test_merge_positions_dd")
        self.assertTrue(True)

    def test_combine_df_kmer_indices(self):
        print("Testing %s" % "test_combine_df_kmer_indices")
        self.assertTrue(True)

    def test_clean_counts_df(self):
        print("Testing %s" % "test_clean_counts_df")
        self.assertTrue(True)

    def test_file_len(self):
        print("Testing %s" % "test_file_len")
        self.assertTrue(True)

    def test_get_counts_from_file(self):
        print("Testing %s" % "test_get_counts_from_file")
        self.assertTrue(True)

    def test_get_counts_dict(self):
        print("Testing %s" % "test_get_counts_dict")
        self.assertTrue(True)

    def test_count_regional_variants(self):
        print("Testing %s" % "test_count_regional_variants")
        self.assertTrue(True)

    def test_count_regional_af(self):
        print("Testing %s" % "test_count_regional_af")
        self.assertTrue(True)

    def test_query_bed_region(self):
        print("Testing %s" % "test_query_bed_region")
        self.assertTrue(True)

    def test_check_bed_regions(self):
        print("Testing %s" % "test_check_bed_regions")
        self.assertTrue(True)

    def test_is_dash(self):
        print("Testing %s" % "test_is_dash")
        self.assertTrue(True)

    def test_check_clinvar(self):
        print("Testing %s" % "test_check_clinvar")
        self.assertTrue(True)
