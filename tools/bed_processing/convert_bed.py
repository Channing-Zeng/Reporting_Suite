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
chains_dirpath = '/group/ngs/src/liftOver'
chains = dict(
    hg38=join(chains_dirpath, 'hg19ToHg38.over.chain.gz'),
    grch37=join(chains_dirpath, 'hg19ToGRCh37.over.chain.gz')
)
if is_local():
    liftover_fpath = '/Users/vladsaveliev/az/liftOver/liftOver'
    chains_dirpath = '/Users/vladsaveliev/az/liftOver'
    chains = dict(
        hg38=join(chains_dirpath, 'hg19ToHg38.over.chain.gz'),
    )


liftover_cmdline = liftover_fpath + ' "{inp_fpath}" {chain_fpath} "{out_fpath}" "{unlifted_fpath}"'


def main(args):
    if len(args) < 2:
        critical('Usage: ' + __file__ + ' InputRootDirectory OutputRootDirectory [Build=hg38]')
        sys.exit(1)

    inp_root = adjust_path(args[0])
    out_root = adjust_path(args[1])

    build = 'hg38'
    if len(args) >= 3:
        build = args[2]

    chain_fpath = chains[build]

    for inp_dirpath, subdirs, files in os.walk(inp_root):
        for fname in files:
            if fname == 'sample1-cn_mops.bed':
                pass
            if fname.endswith('.bed'):
                inp_fpath = adjust_path(join(inp_dirpath, fname))
                print inp_fpath + ': ' + str(count_bed_cols(inp_fpath)) + ' columns'
                out_dirpath = adjust_path(join(out_root, relpath(inp_dirpath, inp_root)))
                safe_mkdir(out_dirpath)
                out_fpath = adjust_path(join(out_dirpath, fname))
                unlifted_fpath = adjust_path(join(out_dirpath, fname + '.unlifted'))

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
    main(sys.argv[1:])