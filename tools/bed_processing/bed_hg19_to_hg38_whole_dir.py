#!/usr/bin/env python

import os
from os.path import abspath, dirname, realpath, join, relpath, splitext, isfile, getsize
from site import addsitedir
from source.targetcov.bam_and_bed_utils import count_bed_cols
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys
from source.logger import critical, info, is_local, err
from source.file_utils import adjust_path, safe_mkdir, verify_file

liftover_fpath = '/group/ngs/src/liftOver/liftOver'
hg38_chain_fpath = '/group/ngs/src/liftOver/hg19ToHg38.over.chain.gz'
grch37_chain_fpath = '/group/ngs/src/liftOver/hg19ToGRCh37.over.chain.gz'
if is_local():
    liftover_fpath = '/Users/vladsaveliev/az/liftOver/liftOver'
    hg38_chain_fpath = '/Users/vladsaveliev/az/liftOver/hg19ToHg38.over.chain.gz'

liftover_cmdline = liftover_fpath + ' "{inp_fpath} {chain_fpath} {out_fpath}" "{unlifted_fpath}"'


def main():
    if len(sys.argv) <= 2:
        critical('Usage: ' + __file__ + ' input_root_directory output_root_directory')
        sys.exit(1)

    inp_root = adjust_path(sys.argv[1])
    out_root = adjust_path(sys.argv[2])
    safe_mkdir(out_root)

    chain_fpath = hg38_chain_fpath

    for inp_dirpath, subdirs, files in os.walk(inp_root):
        for fname in files:
            if fname.endswith('.bed'):
                inp_fpath = join(inp_dirpath, fname)
                print inp_fpath, count_bed_cols(fname) + ' columns'
                out_dirpath = join(out_root, relpath(inp_dirpath, inp_root))
                safe_mkdir(out_dirpath)
                out_fpath = join(out_dirpath, fname)
                unlifted_fpath = join(out_dirpath, fname + '.unlifted')

                cmdline = liftover_cmdline.format(**locals())
                info(cmdline)
                os.system(cmdline)
                verify_file(out_fpath)
                if isfile(unlifted_fpath):
                    if getsize(unlifted_fpath) <= 0:
                        os.remove(unlifted_fpath)
                    else:
                        err('Some records were unlifted and saved to ' + unlifted_fpath)


if __name__ == '__main__':
    main()