#!/usr/bin/env python
import os
import traceback
from collections import defaultdict
from genericpath import isfile
from os.path import join, basename, splitext, isdir

from source.calling_process import call
from source.file_utils import verify_dir, safe_mkdir, verify_file, splitext_plus, adjust_path
from source.logger import warn, info
from source.targetcov.bam_and_bed_utils import index_bam
from source.utils import is_uk

zhongwu_replacement_tx_for_genes = {
    'FANCL':   'NM_018062.3',
    'MET':     'NM_000245.2',
    'CDKN2A':  'NM_000077.4',
    'BRCA1':   'NM_007294.3',
    'MYD88':   'NM_002468.4',
    'PPP2R2A': 'NM_002717.3',
    'RAD51D':  'NM_002878.3',
    'RAD54L':  'NM_003579.3',
    'ESR1':    'NM_000125.3',
    'AKT1':    'NM_005163.2',
    'FGFR3':   'NM_000142.4',
    'CD79B':   'NM_000626.2',
    'CHEK2':   'NM_007194.3',
    'CHEK1':   'NM_001274.5'
}
"""
java -Xms750m -Xmx3g -jar /group/ngs/src/snpEff4.2/snpEff.jar eff -d -canon \
-c /group/ngs/src/snpEff4.2/snpEff.config hg19 2> snpeff_verbose_output_hg19.txt
"""

"""
java -Xms750m -Xmx3g -jar /group/ngs/src/snpEff4.2/snpEff.jar eff -d -canon \
-c /group/ngs/src/snpEff4.2/snpEff.config hg38 2> snpeff_verbose_output_hg38.txt
"""

"""
java -Xms750m -Xmx3g -jar /ngs/RDI/PROGRAMS/snpEff4.2/snpEff.jar eff -d -canon \
-c /ngs/RDI/PROGRAMS/snpEff4.2/snpEff.config hg38 2> snpeff_verbose_output_hg38.txt
"""

def later_tx_version(t1, t2):
    return t1.split('.') > t2.split('.')

def eq_tx(t1, t2):
    return t1.split('.')[0] == t2.split('.')[0]

base_dir = adjust_path('~/Dropbox/az/reference_data/canonical_transcripts')

def get_transcripts_from_snpeff_output(snpeff_output_fname):
    tx_per_gene = defaultdict(list)

    snpeff_output_fpath = join(base_dir, snpeff_output_fname)
    with open(snpeff_output_fpath) as f:
        transcripts_started = False
        for l in f:
            if transcripts_started:
                fs = l.strip().split()
                if 'done.' in fs:
                    break
                if len(fs) != 4 or 'geneName' in fs:
                    continue
                else:
                    g = fs[0]
                    t = fs[2]
                    is_met = False
                    for prev_t in tx_per_gene[g]:
                        if eq_tx(t, prev_t):
                            is_met = True
                            # info(t + ' has the same base tx id as ' + prev_t)
                            if later_tx_version(t, prev_t):
                                # info(t + ' is later than ' + prev_t + ', replacing')
                                tx_per_gene[g][tx_per_gene[g].index(prev_t)] = t
                                break
                            else:
                                # info(t + ' is earlier than ' + prev_t + ', skipping')
                                break
                    if not is_met:
                        tx_per_gene[g].append(t)

            elif 'Canonical transcripts:' in l:
                transcripts_started = True
    return tx_per_gene

info('Getting transcripts for hg19')
hg19_trs_by_gene = get_transcripts_from_snpeff_output('snpeff_verbose_output_hg19.txt')
hg19_genes = set(hg19_trs_by_gene.keys())

info('Getting transcripts for hg38')
hg38_trs_by_gene = get_transcripts_from_snpeff_output('snpeff_verbose_output_hg38.txt')
hg38_genes = set(hg38_trs_by_gene.keys())

info('Genes present only in hg19: ' + str(len(hg19_genes - hg38_genes)))
info('Genes present only in hg38: ' + str(len(hg38_genes - hg19_genes)))

canon_tx_hg19_fpath = join(base_dir, 'canonical_transcripts_hg19.txt')
canon_tx_hg38_fpath = join(base_dir, 'canonical_transcripts_hg38.txt')

not_matching_tr_count = 0
mult_hg19_tx_count = 0
mult_hg38_tx_count = 0
with open(canon_tx_hg19_fpath, 'w') as hg19, open(canon_tx_hg38_fpath, 'w') as hg38:
    for g in hg38_genes:
        if g in zhongwu_replacement_tx_for_genes:
            repl_t = zhongwu_replacement_tx_for_genes[g]

            print g + ' in Zhongwu\'s list as ' + repl_t
            print '  hg19: replacing ' + ', '.join(hg19_trs_by_gene[g]) + ' -> ' + repl_t
            hg19_trs_by_gene[g] = [repl_t]
            print '  hg38:'
            is_met = False
            for t in hg38_trs_by_gene[g]:
                if eq_tx(t, repl_t):
                    print '    ' + t + ' eq to ' + repl_t
                    is_met = True
            if not is_met:
                print '  not met in ' + ', '.join(hg38_trs_by_gene[g])
            print '  replacing ' + ', '.join(hg38_trs_by_gene[g]) + ' -> ' + repl_t
            hg38_trs_by_gene[g] = [repl_t]

        if not set(hg19_trs_by_gene[g]) == set(hg38_trs_by_gene[g]):
            # print 'Transcripts do not match:\n  hg19: ' + ', '.join(hg19_trs_by_gene[g]) + '\n  hg38: ' + ', '.join(hg38_trs_by_gene[g])
            not_matching_tr_count += 1
        if len(hg38_trs_by_gene[g]) > 1:
            # print 'Multiple hg38 tx for ' + g + ': ' + ', '.join(hg38_trs_by_gene[g])
            mult_hg38_tx_count += 1
        if len(hg19_trs_by_gene[g]) > 1:
            # print 'Multiple hg19 tx for ' + g + ': ' + ', '.join(hg38_trs_by_gene[g])
            mult_hg19_tx_count += 1



        if hg19_trs_by_gene[g]:
            hg19.write('\n'.join(hg19_trs_by_gene[g]) + '\n')
        hg38.write('\n'.join(hg38_trs_by_gene[g]) + '\n')

print '---'
print 'total genes:', len(hg38_genes)
print 'not_matching_tr_count:', not_matching_tr_count
print 'mult_hg38_tx_count:', mult_hg38_tx_count
print 'mult_hg19_tx_count:', mult_hg19_tx_count







