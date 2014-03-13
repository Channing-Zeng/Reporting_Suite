# Annotation script that takes 1-3 inputs, first being the vcf file name,
# second being an indicator if the vcf is from bcbio's ensemble pipeline ('true' if true) and
# third being 'RNA' if the vcf is from the rna-seq mutect pipeline
import os
import subprocess
import sys
import shutil


def _call_and_rename(cmdline, save_prev, input_fpath, suffix, stdout=True):
    basepath, ext = os.path.splitext(input_fpath)
    output_fpath = basepath + '.' + suffix + ext

    print ''
    print '*' * 70
    print cmdline
    res = subprocess.call(cmdline.split(), open(output_fpath, 'w') if stdout else None)
    if res != 0:
        print ''
        print '*' * 70
        print 'Command returned status ' + str(res)
        exit(1)

    if save_prev:
        return output_fpath
    else:
        os.remove(input_fpath)
        os.rename(output_fpath, input_fpath)
        return input_fpath


# def _call_and_rename(cmdline, fpath):
#     base_name, ext = os.path.splitext(fpath)
#     output_fpath = base_name + '.tmp' + ext
#     _call(cmdline, open(output_fpath, 'w'))
#     if os.path.isfile(fpath):
#         os.remove(fpath)
#     os.rename(output_fpath, fpath)


def snpsift_annotate(snpsift_jar, db, suffix, vcf_fpath, save_prev):
    cmdline = 'java -jar %s annotate -v %s %s' % (snpsift_jar, db, vcf_fpath)
    return _call_and_rename(cmdline, save_prev, vcf_fpath, suffix, stdout=True)


def snpsift_dbnsfp(snpsift_jar, db, vcf_fpath, save_prev):
    annots = 'SIFT_score,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,' \
             'MutationTaster_score,MutationTaster_pred,MutationAssessor_score,' \
             'MutationAssessor_pred,FATHMM_score,ESP6500_AA_AF,ESP6500_EA_AF,' \
             'Ensembl_geneid,Ensembl_transcriptid'

    cmdline = 'java -jar %s dbnsfp -f %s -v %s %s' % (snpsift_jar, annots, db, vcf_fpath)
    return _call_and_rename(cmdline, save_prev, vcf_fpath, 'db_nsfp', stdout=True)


def snpeff(snpeff_jar, datadir, ref, vcf_fpath, save_prev):
    cmdline = 'java -Xmx4g -jar %s eff -dataDir %s -cancer ' \
              '-noLog -1 -i vcf -o vcf %s %s' % \
              (snpeff_jar, datadir, ref, vcf_fpath)
    return _call_and_rename(cmdline, save_prev, vcf_fpath, 'snpEff', stdout=True)


def rna_editing_sites(db, vcf_fpath, save_prev):
    cmdline = 'vcfannotate -b %s -k RNA_editing_site %s' % (db, vcf_fpath)
    return _call_and_rename(cmdline, save_prev, vcf_fpath, 'edit', stdout=True)


def gatk(gatk_jar, ref_path, vcf_fpath, save_prev):
    base_name, ext = os.path.splitext(vcf_fpath)
    output_fpath = base_name + '.gatk' + ext

    cmdline = 'java -Xmx2g -jar %s -R %s -T VariantAnnotator ' \
              '-o %s --useAllAnnotations --variant %s' % \
              (gatk_jar, ref_path, output_fpath, vcf_fpath)

    annotations = [
        "Coverage", "BaseQualityRankSumTest", "FisherStrand",
        "GCContent", "HaplotypeScore", "HomopolymerRun",
        "MappingQualityRankSumTest", "MappingQualityZero",
        "QualByDepth", "ReadPosRankSumTest", "RMSMappingQuality",
        "DepthPerAlleleBySample"]
    for ann in annotations:
        cmdline += " -A " + ann

    return _call_and_rename(cmdline, save_prev, vcf_fpath, 'gatk', stdout=False)


def annotate_hg19(sample_fpath, save_intermediate, snp_eff, gatk_dir, is_rna=False, is_ensemble=False):
    ref_name = 'hg19'
    ref_path = '/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa'
    dbsnp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/dbsnp_137.vcf'
    cosmic_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cosmic-v67_20131024-hg19.vcf'
    db_nsfp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/dbNSF/dbNSFP2.3/dbNSFP2.3.txt.gz'
    snpeff_datadir = '/ngs/reference_data/genomes/Hsapiens/hg19/snpeff'
    annot_track = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/Human_AG_all_hg19_INFO.bed'

    annotate(sample_fpath, save_intermediate, is_rna, is_ensemble,
             snp_eff, gatk_dir,
             ref_name, ref_path,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track)


def annotate_GRCh37(sample_fpath, save_intermediate, snp_eff, gatk_dir, is_rna=False, is_ensemble=False):
    ref_name = 'GRCh37'
    ref_path = '/ngs/reference_data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
    dbsnp_db = '/ngs/reference_data/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf'
    cosmic_db = '/ngs/reference_data/genomes/Hsapiens/GRCh37/variation/cosmic-v67_20131024-GRCh37.vcf'
    db_nsfp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/dbNSF/dbNSFP2.3/dbNSFP2.3.txt.gz'
    snpeff_datadir = '/ngs/reference_data/genomes/Hsapiens/GRCh37/snpeff'
    annot_track = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/Human_AG_all_hg19_INFO.bed'

    annotate(sample_fpath, save_intermediate, is_rna, is_ensemble,
             snp_eff, gatk_dir,
             ref_name, ref_path,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track)


def annotate(sample_fpath, save_intermediate, is_rna, is_ensemble,
             snpeff_dirpath, gatk_dirpath,
             ref_name, ref_path,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track):
    # sample_dbsnp_fpath = sample_basepath + '.dbsnp' + ext
    # sample_cosmic_fpath = sample_basepath + '.cosmic' + ext
    # sample_snpeff_fpath = sample_basepath + '.snpeff' + ext
    # sample_dbnsfp_fpath = sample_basepath + '.dbnsf' + ext
    # sample_rna_edit_sites_fpath = sample_basepath + '.RNAeditSites' + ext

    if is_ensemble:
        sample_fname = os.path.basename(sample_fpath)
        sample_basename, ext = os.path.splitext(sample_fname)
        cmdline = 'vcf-subset -c %s -e %s' % (sample_basename.replace('-ensemble', ''), sample_fpath)

        sample_fpath = _call_and_rename(cmdline, save_intermediate, sample_fpath, '.ensm', True)

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

    snpsift_jar = os.path.join(snpeff_dirpath, 'SnpSift.jar')
    snpeff_jar = os.path.join(snpeff_dirpath, 'snpEff.jar')
    gatk_jar = os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

    sample_fpath = snpsift_annotate(snpsift_jar, dbsnp_db, 'dbsnp', sample_fpath, save_intermediate)
    sample_fpath = snpsift_annotate(snpsift_jar, cosmic_db, 'cosmic', sample_fpath, save_intermediate)
    sample_fpath = snpsift_dbnsfp(snpsift_jar, db_nsfp_db, sample_fpath, save_intermediate)
    sample_fpath = snpeff(snpeff_jar, snpeff_datadir, ref_name, sample_fpath, save_intermediate)
    if is_rna:
        sample_fpath = rna_editing_sites(annot_track, sample_fpath, save_intermediate)
    sample_fpath = gatk(gatk_jar, ref_path, sample_fpath, save_intermediate)


    cmdline = 'cat %s | ' \
              'vcfEffOnePerLine.pl | ' \
              'java -jar %s extractFields - ' \
              'CHROM POS ID CNT GMAF REF ALT QUAL FILTER TYPE ' \
              '"EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].CODON" ' \
              '"EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" ' \
              'dbNSFP_SIFT_score dbNSFP_Polyphen2_HVAR_score ' \
              'dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred ' \
              'dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred ' \
              'dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred ' \
              'dbNSFP_FATHMM_score dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid ' \
              'dbNSFP_Uniprot_acc dbNSFP_1000Gp1_AC dbNSFP_1000Gp1_AF ' \
              'dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF KGPROD PM PH3 ' \
              'AB AC AF DP FS GC HRun HaplotypeScore MQ0 QA QD ReadPosRankSum ' \
              'set' % (sample_fpath, snpsift_jar)

    sample_fpath = _call_and_rename(cmdline, save_intermediate, sample_fpath, 'extract', stdout=True)
    os.rename(sample_fpath, os.path.splitext(sample_fpath)[0] + '.tsv')


def remove_quotes(str):
    if str and str[0] == '"':
        str = str[1:]
    if str and str[-1] == '"':
        str = str[:-1]
    return str


def split_genotypes(sample_fpath, save_intermediate):
    print 'Splitting genotypes'

    sample_basepath, ext = os.path.splitext(sample_fpath)
    result_fpath = sample_basepath + '.split' + ext

    with open(sample_fpath) as vcf, open(result_fpath, 'w') as out:
        for i, line in enumerate(vcf):
            clean_line = line.strip()
            if not clean_line or clean_line[0] == '#':
                out.write(line)
            else:
                tokens = line.split()
                id_field = remove_quotes(tokens[2])
                alt_field = remove_quotes(tokens[4])

                ids = id_field.split(',')
                alts = alt_field.split(',')
                if len(ids) != len(alts):
                    print 'Number of IDs is not equal to the number of ALTs: ' + str(i) + '. ' + line
                    continue
                if len(ids) > 1:
                    print 'Splitting ' + str(i) + '. ' + line
                    for id, alt in zip(ids, alts):
                        line = '\t'.join(tokens[:2] + [id] + [tokens[3]] + [alt] + tokens[5:]) + '\n'
                        out.write(line)
                else:
                    out.write(line)

    if save_intermediate:
        return result_fpath
    else:
        os.remove(sample_fpath)
        os.rename(result_fpath, sample_fpath)
        return sample_fpath


snp_eff = '/group/ngs/src/snpEff/snpEff3.5/'
gatk_dir = '/opt/az/broadinstitute/gatk/1.6'

if __name__ == '__main__':
    args = sys.argv[1:]

    flags = ['-rna', '-ensemble', '-split', '-intermediate', '-to-valid']
    rna = '-rna' in args
    ensemble = '-ensemble' in args
    do_split_genotypes = '-split' in args
    save_intermediate = '-intermediate' in args
    to_valid = '-to-valid' in args

    if len(args) < 1:
        print >> sys.stderr, \
            'Usage: python ' + __file__ + ' sample.vcf [result.vcf] ' \
            '[-split] [-ensemble] [-rna] [-intermediate] [-to_valid]'
        exit(1)

    sample_fpath = args[0]
    assert os.path.isfile(sample_fpath), \
        os.path.realpath(sample_fpath) + ' does not exists or is not a file'

    if len(args) > 1 and args[1] not in flags:
        result_fpath = args[1]
    else:
        result_fpath = os.path.join(os.getcwd(), os.path.basename(sample_fpath))

    if result_fpath != sample_fpath:
        if os.path.exists(result_fpath):
            os.remove(result_fpath)
        shutil.copyfile(sample_fpath, result_fpath)
        sample_fpath = result_fpath

    sample_basedir = os.path.dirname(sample_fpath)
    sample_basepath, ext = os.path.splitext(sample_fpath)

    if do_split_genotypes:
        sample_fpath = split_genotypes(sample_fpath, save_intermediate)

    print 'Please, run this before start:'
    print '   source /etc/profile.d/modules.sh'
    print '   module load bcbio-nextgen/0.7.6'
    print ''
    print 'In Waltham, run this as well:'
    print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl'

    annotate_hg19(sample_fpath, save_intermediate, snp_eff, gatk_dir, rna, ensemble)