from pkg_resources import resource_filename


def get_autosome_names_grch38():
    """
    :return: a list of the primary autosome names for build 38 (eg. chr1)
    """
    return ['chr1',
            'chr2',
            'chr3',
            'chr4',
            'chr5',
            'chr6',
            'chr7',
            'chr8',
            'chr9',
            'chr10',
            'chr11',
            'chr12',
            'chr13',
            'chr14',
            'chr15',
            'chr16',
            'chr17',
            'chr18',
            'chr19',
            'chr20',
            'chr21',
            'chr22']


def get_grch38_chroms():
    names = []
    lengths = []
    filepath = resource_filename('eskedit', 'chrom_names_lengths.csv')
    with open(filepath, 'r') as cnames:
        for line in cnames.readlines():
            split = line.split(',')
            names.append(split[0])
            lengths.append(int(split[1]))
    return dict(zip(names, lengths))


def get_autosome_lengths_grch38():
    """
    :return: a dictionary mapping chromosome names to number of nucleotides
    """
    return {'chr1': 248956422,
            'chr2': 242193529,
            'chr3': 198295559,
            'chr4': 190214555,
            'chr5': 181538259,
            'chr6': 170805979,
            'chr7': 159345973,
            'chr8': 145138636,
            'chr9': 138394717,
            'chr10': 133797422,
            'chr11': 135086622,
            'chr12': 133275309,
            'chr13': 114364328,
            'chr14': 107043718,
            'chr15': 101991189,
            'chr16': 90338345,
            'chr17': 83257441,
            'chr18': 80373285,
            'chr19': 58617616,
            'chr20': 64444167,
            'chr21': 46709983,
            'chr22': 50818468}


def get_test_path():
    """
    :return: a tuple with vcf path and reference fasta path for testing
    """
    vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz"
    fasta_path = "/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa"
    return vcf_path, fasta_path
