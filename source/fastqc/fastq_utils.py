from itertools import izip, product
import os
from os.path import splitext, dirname, join, basename
import random
import gzip
from source import info
from source.calling_process import call_check_output
from source.file_utils import open_gzipsafe, file_transaction, file_exists, intermediate_fname, verify_file, add_suffix
from source.logger import critical, debug


LIMIT = 5*1000*1000


def downsample(cnf, fastq_L_fpath, fastq_R_fpath, N, output_dir, suffix=None, quick=False):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    quick=True will just grab the first N reads rather than do a true
    downsampling
    """
    l_out_fpath = join(output_dir, add_suffix(basename(fastq_L_fpath), suffix or 'subset'))
    r_out_fpath = join(output_dir, add_suffix(basename(fastq_R_fpath), suffix or 'subset'))
    if cnf.reuse_intermediate and verify_file(l_out_fpath, silent=True) and verify_file(r_out_fpath, silent=True):
        info(l_out_fpath + ' and ' + r_out_fpath + ' exist, reusing.')
        return l_out_fpath, r_out_fpath

    N = int(N)
    records_num = N
    if quick:
        rand_records = range(N)
    else:
        records_num = sum(1 for _ in open(fastq_L_fpath)) / 4
        if records_num > LIMIT:
            info('The number of reads is higher than ' + str(LIMIT) +
                 ', sampling from only first ' + str(LIMIT))
            records_num = LIMIT
        info(str(records_num) + ' reads in ' + fastq_L_fpath)
        if records_num < N:
            info('...it is less than ' + str(N) + ', so no downsampling.')
            return fastq_L_fpath, fastq_R_fpath
        else:
            info('Downsampling to ' + str(N))
            rand_records = sorted(random.sample(xrange(records_num), N))

    info('Opening ' + fastq_L_fpath)
    fh1 = open_gzipsafe(fastq_L_fpath)
    info('Opening ' + fastq_R_fpath)
    fh2 = open_gzipsafe(fastq_R_fpath) if fastq_R_fpath else None

    out_files = (l_out_fpath, r_out_fpath) if r_out_fpath else (l_out_fpath)

    written_records = 0
    with file_transaction(cnf.work_dir, out_files) as tx_out_files:
        if isinstance(tx_out_files, basestring):
            tx_out_f1 = tx_out_files
        else:
            tx_out_f1, tx_out_f2 = tx_out_files
        info('Opening ' + str(tx_out_f1) + ' to write')
        sub1 = open_gzipsafe(tx_out_f1, "w")
        info('Opening ' + str(tx_out_f2) + ' to write')
        sub2 = open_gzipsafe(tx_out_f2, "w") if r_out_fpath else None
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
            written_records += 1
            rec_no += 1
            if written_records % 100000 == 0:
                info('written ' + str(written_records) + ', rec_no ' + str(rec_no))
            if rec_no > records_num:
                info('Reached the limit of ' + str(records_num), ' read lines, stopping.')
                break
        info('Done, written ' + str(written_records) + ', rec_no ' + str(rec_no))
        fh1.close()
        sub1.close()
        if fastq_R_fpath:
            fh2.close()
            sub2.close()

    info('Done downsampling, saved to ' + l_out_fpath + ' and ' + r_out_fpath + ', total ' + str(written_records) + ' paired reads written')
    return l_out_fpath, r_out_fpath