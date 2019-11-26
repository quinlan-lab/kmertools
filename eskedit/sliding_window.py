import pandas as pd
import multiprocessing as mp
import eskedit as ek

"""
This package will search a genome reference sliding window with variants aligned.

This currently works with 3, 5, and 7mer contexts. The expected input is a tab separated
file with the following fields in this order:
CHROM   POS REF ALT`
    1. CHROM - chromosome should be as reported by grch38
    2. POS - position on chromosome aligned to hg38
    3. REF - reference allele
    4. ALT - alternate allele from VCF
    * Additional fields after this will be ignored
"""
WINDOW_SIZE = 0


def process_chunk(df):
    ref = df.iloc[:, 0].to_numpy()
    alt = df.iloc[:, 1].to_numpy()


def sliding_window_from_csv(var_csv, window_size, nprocs=None):
    global WINDOW_SIZE
    WINDOW_SIZE = window_size
    if nprocs is None:
        nprocs = mp.cpu_count()
    chunk_size = ek.file_len(var_csv) // nprocs
    df_chunks = pd.read_csv(var_csv, sep='\t', chunksize=chunk_size)
    pool = mp.pool(nprocs)
    results = [result.get() for result in [pool.starmap_async(process_chunk, df_chunks)]]
    pool.close()
    pass
