# kmertools
This is currently being developed and has some instability.

Requires the installation of `cyvcf2`, `numpy`, `pandas`, and `pyfaidx`.
```
pip install cyvcf2 numpy pandas pyfaidx 
```

See usage below. Note that the command requires either `--query` or `--train` as a flag.

```
usage: test.py [-h] [--query] [--train] [--kmer_size KMER_SIZE]
               [--bedpath BED_PATH] [--vcfpath VCF_PATH]
               [--fastapath FASTA_PATH] [--countspath COUNTSPATH]
               [--nprocs NPROCS] [--invert]

optional arguments:
  -h, --help            show this help message and exit
  
  --query               Query regions specified in bedfile for an expected
                        number of mutations based on provided counts data.
                        
  --train               Build a counts table based on a k-mer model

  --kmer_size KMER_SIZE, -k KMER_SIZE      Length of k-mer motif
                        
  --bedpath BED_PATH, -b BED_PATH          Path to bed file containing genomic regions
                        
  --vcfpath VCF_PATH, -v VCF_PATH          Path to vcf file containing variants
                        
  --fastapath FASTA_PATH, -f FASTA_PATH    Path to reference genome fasta
                        
  --countspath COUNTSPATH, -c COUNTSPATH   Path to counts table
                        
  --nprocs NPROCS, -N NPROCS               Number of processes to use (default=1)
  
  --invert                                 If flag is present, will invert regions given in bed
                                           file.
                        
```


## Examples

Train a k-mer model:
```
bedfile="merged_exons_grch38.bed"
vcf="gnomad.genomes.r3.0.sites.vcf.bgz"
fasta="hg38.fa"
NTASKS=20

python test.py --train -k 7 -b $bedfile -v $vcf -f $fasta -N $NTASKS --invert
```
This will return a counts file which will be used for the expectation value

Query a set of regions using a trained model:
```
kmer=7
bedfile="bed_regions_you_want_to_calculate_expectation_for.bed"
vcf="gnomad.genomes.r3.0.sites.vcf.bgz"
fasta="hg38.fa"
counts="7mer_relative_freq_noncoding.csv"
NTASKS=20

python test.py --query -k $kmer -b $bedfile -v $vcf -f $fasta -c $counts -N $NTASKS
```
