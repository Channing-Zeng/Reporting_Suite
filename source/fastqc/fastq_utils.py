from itertools import izip, product
import os
from os.path import splitext, dirname, join
import random
import gzip
from source.file_utils import open_gzipsafe, file_transaction, file_exists, intermediate_fname


def downsample(cnf, fastq_L_fpath, fastq_R_fpath, N, quick=False):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    quick=True will just grab the first N reads rather than do a true
    downsampling
    """
    N = int(N)
    if quick:
        rand_records = range(N)
    else:
        records = sum(1 for _ in open(fastq_L_fpath)) / 4
        N = records if N > records else N
        rand_records = random.sample(xrange(records), N)

    fh1 = open_gzipsafe(fastq_L_fpath)
    fh2 = open_gzipsafe(fastq_R_fpath) if fastq_R_fpath else None
    outf1 = intermediate_fname(cnf, fastq_L_fpath, 'subset')
    outf2 = intermediate_fname(cnf, fastq_R_fpath, 'subset')

    if file_exists(outf1):
        if not outf2:
            return outf1, outf2
        elif file_exists(outf2):
            return outf1, outf2

    out_files = (outf1, outf2) if outf2 else (outf1)

    with file_transaction(cnf.work_dir, out_files) as tx_out_files:
        if isinstance(tx_out_files, basestring):
            tx_out_f1 = tx_out_files
        else:
            tx_out_f1, tx_out_f2 = tx_out_files
        sub1 = open_gzipsafe(tx_out_f1, "w")
        sub2 = open_gzipsafe(tx_out_f2, "w") if outf2 else None
        rec_no = - 1
        for rr in rand_records:
            while rec_no < rr:
                rec_no += 1
                for i in range(4): fh1.readline()
                if fh2:
                    for i in range(4): fh2.readline()
            for i in range(4):
                sub1.write(fh1.readline())
                if sub2:
                    sub2.write(fh2.readline())
            rec_no += 1
        fh1.close()
        sub1.close()
        if fastq_R_fpath:
            fh2.close()
            sub2.close()

    return outf1, outf2