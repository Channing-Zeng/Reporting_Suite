import sys
import os
import shutil

from source.playground.no_conf_annotate import annotate_GRCh37


snp_eff = '/group/ngs/src/snpEff/snpEff3.5/'
gatk_dir = '/opt/az/broadinstitute/gatk/1.6'

if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) < 1:
        print >> sys.stderr, 'Usage: python ' + __file__ + ' input_directory_tree'
        exit(1)

    input_dir = args[0]
    dir_name = os.path.basename(os.path.dirname(input_dir))
    if os.path.exists(dir_name):
        print >> sys.stderr, 'Directory already exists'
        exit(1)

    shutil.copytree(input_dir, dir_name)

    for dirpath, dnames, fnames in os.walk(dir_name):
        for fname in fnames:
            if fname.endswith('.vcf'):
                base_name, ext = os.path.splitext(fname)
                annotated_fpath = os.path.join(dirpath, base_name + 'ANNOTATED' + ext)
                input_fpath = os.path.join(dirpath, fname)

                annotate_GRCh37(input_fpath, snp_eff, gatk_dir)