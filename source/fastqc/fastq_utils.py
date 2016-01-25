from itertools import izip, product
import os
from os.path import splitext, dirname, join
import random
import gzip
from source import info
from source.file_utils import open_gzipsafe, file_transaction, file_exists, intermediate_fname, verify_file


def downsample(cnf, fastq_L_fpath, fastq_R_fpath, N, quick=False):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    quick=True will just grab the first N reads rather than do a true
    downsampling
    """
    outf1 = intermediate_fname(cnf, fastq_L_fpath, 'subset')
    outf2 = intermediate_fname(cnf, fastq_R_fpath, 'subset')
    if cnf.reuse_intermediate and verify_file(outf1, silent=True) and verify_file(outf2, silent=True):
        info(outf1 + ' and ' + outf2 + ' exist, reusing.')
        return outf1, outf2

    N = int(N)
    if quick:
        rand_records = range(N)
    else:
        records_num = sum(1 for _ in open(fastq_L_fpath)) / 4
        if N > records_num:
            N = records_num
            info('Downsampling to ' + str(N))
            rand_records = range(N)
        else:
            rand_records = random.sample(xrange(records_num), N)

    fh1 = open_gzipsafe(fastq_L_fpath)
    fh2 = open_gzipsafe(fastq_R_fpath) if fastq_R_fpath else None

    out_files = (outf1, outf2) if outf2 else (outf1)

    with file_transaction(cnf.work_dir, out_files) as tx_out_files:
        if isinstance(tx_out_files, basestring):
            tx_out_f1 = tx_out_files
        else:
            tx_out_f1, tx_out_f2 = tx_out_files
        sub1 = open_gzipsafe(tx_out_f1, "w")
        sub2 = open_gzipsafe(tx_out_f2, "w") if outf2 else None
        rec_no = -1
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