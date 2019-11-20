import eskedit as ek
vcf_path = '/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz'
fasta_path = '/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa'

results = ek.get_kmer_context(vcf_path, fasta_path, 3, nprocs=6)
# results should be  list of dictionaries
transitions = []
positions = []
for coll in results:
    for k, v in coll:
        if 'position' in k:
            positions.append(v)
        else:
            transitions.append(v)

# s_transitions = defaultdict(Counter)
# s_positions = defaultdict(Variant)
ek.merge_positions_dd(positions, outfile='singletons_test.bed')
ek.merge_transitions_ddc(transitions)
