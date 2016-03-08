#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import sys
from optparse import SUPPRESS_HELP, OptionParser

import source
from source.file_utils import iterate_file, open_gzipsafe, intermediate_fname, verify_file, adjust_path
from source.logger import info, err, warn, step_greetings, critical


def get_args():
    info(' '.join(sys.argv))
    info()
    parser = OptionParser()
    parser.add_option('-o', dest='output_file')

    parser.set_usage('Usage: ' + __file__ + ' dbsnp.vcf.gz -o output_fpath')

    (opts, args) = parser.parse_args()
    if len(args) < 1:
        critical("Provide the first argument - path to dbsnp VCF")

    vcf2txt_res_fpath = verify_file(args[0])

    if not opts.output_file:
        critical('Please, specify the output fpath with -o')

    info()

    return vcf2txt_res_fpath, adjust_path(opts.output_file)


def main():
    dbsnp_fpath, out_fpath = get_args()

    info('-' * 70)
    info('Reading ' + dbsnp_fpath + ', writing to ' + out_fpath)

    with open_gzipsafe(dbsnp_fpath) as dbsnp, open(out_fpath, 'w') as out:
        for l in dbsnp:
            if l.startswith('#'):
                continue
            fs = l.replace('\n', '').split('\t')
            assert len(fs) == 8, str(fs)

            chrom, pos, rsid, ref, alt, _, _, inf = l.replace('\n', '').split('\t')
            alts = alt.split(',')
            if len(alts) > 1:
                caf = next((kv.split('=')[1] for kv in inf.split(';') if kv.split('=')[0] == 'CAF'), None)
                if caf:
                    cafs = caf.replace('[', '').replace(']', '').split(',')[1:]
                    assert len(cafs) == len(alts), l
                    for alt, caf in zip(alts, cafs):
                        if caf != '.':
                            l = '\t'.join([chrom, pos, rsid, ref, alt, caf]) + '\n'
                            out.write(l)

    info()
    info('Saved to ' + out_fpath)


if __name__ == '__main__':
    main()
