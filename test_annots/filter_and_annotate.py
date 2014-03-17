# Annotation script that takes 1-3 inputs, first being the vcf file name,
# second being an indicator if the vcf is from bcbio's ensemble pipeline ('true' if true) and
# third being 'RNA' if the vcf is from the rna-seq mutect pipeline
import os
import subprocess
import sys
import shutil


def log_print(msg='', fpath=None):
    print msg
    if fpath:
        open(fpath, 'a').write(msg + '\n')


def _call_and_rename(cmdline, input_fpath, suffix, log_fpath=None, save_prev=False, stdout=True):
    basepath, ext = os.path.splitext(input_fpath)
    output_fpath = basepath + '.' + suffix + ext

    log_print('', log_fpath)
    log_print('*' * 70, log_fpath)
    log_print(cmdline, log_fpath)
    res = subprocess.call(cmdline.split(),
                          stdout=open(output_fpath, 'w') if stdout else open(log_fpath, 'a') if log_fpath else None,
                          stderr=open(log_fpath, 'a') if log_fpath else None)
    log_print('', log_fpath)
    if res != 0:
        log_print('Command returned status ' + str(res) + ('. Log in ' + log_fpath if log_fpath else ''),
                  log_fpath)
        return input_fpath
    else:
        log_print('Saved to ' + output_fpath, log_fpath)
        if log_fpath:
            print 'Log in ' + log_fpath

    if not save_prev:
        os.remove(input_fpath)
    log_print('Now processing ' + output_fpath, log_fpath)
    return output_fpath


def snpsift_annotate(snpsift_jar, db, suffix, vcf_fpath, save_prev):
    cmdline = 'java -jar %s annotate -v %s %s' % (snpsift_jar, db, vcf_fpath)
    return _call_and_rename(cmdline, vcf_fpath, suffix, log_fpath, save_prev, stdout=True)


def snpsift_dbnsfp(snpsift_jar, db, vcf_fpath, save_prev):
    annots = 'SIFT_score,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,' \
             'MutationTaster_score,MutationTaster_pred,MutationAssessor_score,' \
             'MutationAssessor_pred,FATHMM_score,ESP6500_AA_AF,ESP6500_EA_AF,' \
             'Ensembl_geneid,Ensembl_transcriptid'

    cmdline = 'java -jar %s dbnsfp -f %s -v %s %s' % (snpsift_jar, annots, db, vcf_fpath)
    return _call_and_rename(cmdline, vcf_fpath, 'db_nsfp', log_fpath, save_prev, stdout=True)


def snpeff(snpeff_jar, datadir, ref, vcf_fpath, save_prev):
    cmdline = 'java -Xmx4g -jar %s eff -dataDir %s -cancer ' \
              '-noLog -1 -i vcf -o vcf %s %s' % \
              (snpeff_jar, datadir, ref, vcf_fpath)
    return _call_and_rename(cmdline, vcf_fpath, 'snpEff', log_fpath, save_prev, stdout=True)


def rna_editing_sites(db, vcf_fpath, save_prev):
    cmdline = 'vcfannotate -b %s -k RNA_editing_site %s' % (db, vcf_fpath)
    return _call_and_rename(cmdline, vcf_fpath, 'edit', log_fpath, save_prev, stdout=True)


def gatk(gatk_jar, ref_path, vcf_fpath, save_prev):
    base_name, ext = os.path.splitext(vcf_fpath)
    output_fpath = base_name + '.gatk' + ext

    cmdline = 'java -Xmx2g -jar %s -R %s -T VariantAnnotator ' \
              '-o %s --variant %s' % \
              (gatk_jar, ref_path, output_fpath, vcf_fpath)

    annotations = [
        "DepthOfCoverage", "BaseQualityRankSumTest", "FisherStrand",
        "GCContent", "HaplotypeScore", "HomopolymerRun",
        "MappingQualityRankSumTest", "MappingQualityZero",
        "QualByDepth", "ReadPosRankSumTest", "RMSMappingQuality",
        "DepthPerAlleleBySample"]
    for ann in annotations:
        cmdline += " -A " + ann

    return _call_and_rename(cmdline, vcf_fpath, 'gatk', log_fpath, save_prev, stdout=False)


def annotate_hg19(sample_fpath, snp_eff_dir, snp_eff_scritps, gatk_dir, save_intermediate=False,
                  log_fpath=None, is_rna=False, is_ensemble=False):
    ref_name = 'hg19'
    ref_path = '/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa'
    dbsnp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/dbsnp_137.vcf'
    cosmic_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cosmic-v67_20131024-hg19.vcf'
    db_nsfp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/dbNSF/dbNSFP2.3/dbNSFP2.3.txt.gz'
    snpeff_datadir = '/ngs/reference_data/genomes/Hsapiens/hg19/snpeff'
    annot_track = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/Human_AG_all_hg19_INFO.bed'

    annotate(sample_fpath,
             snp_eff_dir, snp_eff_scritps, gatk_dir,
             ref_name, ref_path,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track,
             log_fpath, save_intermediate, is_rna, is_ensemble)


def annotate_GRCh37(sample_fpath, snp_eff_dir, snp_eff_scripts, gatk_dir, save_intermediate=False,
                    log_fpath=None, is_rna=False, is_ensemble=False):
    ref_name = 'GRCh37'
    ref_path = '/ngs/reference_data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
    dbsnp_db = '/ngs/reference_data/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf'
    cosmic_db = '/ngs/reference_data/genomes/Hsapiens/GRCh37/variation/cosmic-v67_20131024-GRCh37.vcf'
    db_nsfp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/dbNSF/dbNSFP2.3/dbNSFP2.3.txt.gz'
    snpeff_datadir = '/ngs/reference_data/genomes/Hsapiens/GRCh37/snpeff'
    annot_track = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/Human_AG_all_hg19_INFO.bed'

    annotate(sample_fpath,
             snp_eff_dir, snp_eff_scripts, gatk_dir,
             ref_name, ref_path,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track,
             log_fpath, save_intermediate, is_rna, is_ensemble)


def annotate(sample_fpath,
             snp_eff_dirpath, snp_eff_scripts, gatk_dirpath,
             ref_name, ref_path,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track,
             log_fpath, save_intermediate, is_rna, is_ensemble):
    # sample_dbsnp_fpath = sample_basepath + '.dbsnp' + ext
    # sample_cosmic_fpath = sample_basepath + '.cosmic' + ext
    # sample_snpeff_fpath = sample_basepath + '.snpeff' + ext
    # sample_dbnsfp_fpath = sample_basepath + '.dbnsf' + ext
    # sample_rna_edit_sites_fpath = sample_basepath + '.RNAeditSites' + ext

    if is_ensemble:
        sample_fname = os.path.basename(sample_fpath)
        sample_basename, ext = os.path.splitext(sample_fname)
        cmdline = 'vcf-subset -c %s -e %s' % (sample_basename.replace('-ensemble', ''), sample_fpath)

        sample_fpath = _call_and_rename(cmdline, sample_fpath, '.ensm', log_fpath, save_intermediate, True)

    if is_rna:
        sample_basepath, ext = os.path.splitext(sample_fpath)
        pass_sample_fpath = sample_basepath + '.pass' + ext
        with open(sample_fpath) as sample, open(pass_sample_fpath, 'w') as pass_sample:
            for line in sample.readlines():
                if 'REJECT' not in line:
                    pass_sample.write(line)
        if save_intermediate:
            sample_fpath = pass_sample_fpath
        else:
            os.remove(sample_fpath)
            os.rename(pass_sample_fpath, sample_fpath)

    snpsift_jar = os.path.join(snp_eff_dirpath, 'SnpSift.jar')
    snpeff_jar = os.path.join(snp_eff_dirpath, 'snpEff.jar')
    vcfoneperline = os.path.join(snp_eff_scripts, 'vcfEffOnePerLine.pl')
    gatk_jar = os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

    sample_fpath = snpsift_annotate(snpsift_jar, dbsnp_db, 'dbsnp', sample_fpath, save_intermediate)
    sample_fpath = snpsift_annotate(snpsift_jar, cosmic_db, 'cosmic', sample_fpath, save_intermediate)
    sample_fpath = snpsift_dbnsfp(snpsift_jar, db_nsfp_db, sample_fpath, save_intermediate)
    if is_rna:
        sample_fpath = rna_editing_sites(annot_track, sample_fpath, save_intermediate)
    sample_fpath = gatk(gatk_jar, ref_path, sample_fpath, save_intermediate)
    sample_fpath = snpeff(snpeff_jar, snpeff_datadir, ref_name, sample_fpath, save_intermediate)


    cmdline = 'cat ' + sample_fpath + ' | ' \
              'perl ' + vcfoneperline + ' | ' \
              'java -jar ' + snpsift_jar + ' extractFields - ' \
              'CHROM POS ID CNT GMAF REF ALT QUAL FILTER TYPE ' \
              '"EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].CODON" ' \
              '"EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" ' \
              '"EFF[*].FUNCLASS" "EFF[*].BIOTYPE" "EFF[*].CODING" ' \
              '"EFF[*].TRID" "EFF[*].RANK" ' \
              'dbNSFP_SIFT_score dbNSFP_Polyphen2_HVAR_score ' \
              'dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred ' \
              'dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred ' \
              'dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred ' \
              'dbNSFP_FATHMM_score dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid ' \
              'dbNSFP_Uniprot_acc dbNSFP_1000Gp1_AC dbNSFP_1000Gp1_AF ' \
              'dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF KGPROD PM PH3 ' \
              'AB AC AF DP FS GC HRun HaplotypeScore ' \
              'G5 CDA GMAF GENEINFO OM DB GENE AA CDS ' \
              'MQ0 QA QD ReadPosRankSum '

    sample_fpath = _call_and_rename(cmdline, sample_fpath, 'extract', log_fpath, save_intermediate, stdout=True)
    os.rename(sample_fpath, os.path.splitext(sample_fpath)[0] + '.tsv')


def remove_quotes(str):
    if str and str[0] == '"':
        str = str[1:]
    if str and str[-1] == '"':
        str = str[:-1]
    return str


def split_genotypes(sample_fpath, result_fpath, save_intermediate):
    with open(sample_fpath) as vcf, open(result_fpath, 'w') as out:
        for i, line in enumerate(vcf):
            clean_line = line.strip()
            if not clean_line or clean_line[0] == '#':
                out.write(line)
            else:
                tokens = line.split()
                alt_field = remove_quotes(tokens[4])
                alts = alt_field.split(',')
                if len(alts) > 1:
                    for alt in set(alts):
                        line = '\t'.join(tokens[:2] + ['.'] + [tokens[3]] + [alt] + tokens[5:]) + '\n'
                        out.write(line)
                else:
                    line = '\t'.join(tokens[:2] + ['.'] + tokens[3:]) + '\n'
                    out.write(line)

    if save_intermediate:
        return result_fpath
    else:
        os.remove(sample_fpath)
        os.rename(result_fpath, sample_fpath)
        return sample_fpath


snpeff_dir = '/group/ngs/src/snpEff/snpEff3.5/'
snpeff_scripts = '/group/ngs/src/snpEff/snpEff3.5/scripts'
gatk_dir = '/opt/az/broadinstitute/gatk/1.6'

if __name__ == '__main__':
    args = sys.argv[1:]

    flags = ['-rna', '-ensemble', '-split', '-to-valid']
    rna = '-rna' in args
    ensemble = '-ensemble' in args
    do_split_genotypes = '-split' in args
    # save_intermediate = '-intermediate' in args
    to_valid = '-to-valid' in args

    if len(args) < 1:
        print >> sys.stderr, \
            'Usage: python ' + __file__ + ' sample.vcf [result_dir] ' \
            '[-split] [-ensemble] [-rna] [-to_valid]'
        exit(1)

    sample_fpath = os.path.realpath(args[0])
    assert os.path.isfile(sample_fpath), sample_fpath + ' does not exists or is not a file.'

    if len(args) > 1 and args[1] not in flags:
        result_dir = os.path.realpath(args[1])
    else:
        result_dir = os.getcwd()

    sample_fname = os.path.basename(sample_fpath)
    sample_basename, ext = os.path.splitext(sample_fname)

    if result_dir != os.path.realpath(os.path.dirname(sample_fpath)):
        new_sample_fpath = os.path.join(result_dir, sample_fname)
        if os.path.exists(new_sample_fpath):
            os.remove(new_sample_fpath)
        shutil.copyfile(sample_fpath, new_sample_fpath)
        sample_fpath = new_sample_fpath

    log_fpath = os.path.join(os.path.dirname(sample_fpath), sample_basename + '.log')
    if os.path.isfile(log_fpath):
        os.remove(log_fpath)

    log_print('Writing into ' + result_dir, log_fpath)

    print 'Note: please, load modules before start:'
    print '   source /etc/profile.d/modules.sh'
    print '   module load java'
    print '   module load perl'
    print ''
    print 'In Waltham, run this as well:'
    print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl'

    if do_split_genotypes:
        sample_basepath, ext = os.path.splitext(sample_fpath)
        result_fpath = sample_basepath + '.split' + ext
        log_print('', log_fpath)
        log_print('*' * 70, log_fpath)
        log_print('Splitting genotypes.', log_fpath)
        sample_fpath = split_genotypes(sample_fpath, result_fpath, save_intermediate=True)
        log_print('Saved to ' + result_fpath, log_fpath)
        log_print('', log_fpath)

    annotate_GRCh37(sample_fpath, snpeff_dir, snpeff_scripts, gatk_dir, save_intermediate=True,
                  log_fpath=log_fpath, is_rna=rna, is_ensemble=ensemble)