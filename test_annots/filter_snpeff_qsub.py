# Annotation script that takes 1-3 inputs, first being the vcf file name,
# second being an indicator if the vcf is from bcbio's ensemble pipeline ('true' if true) and
# third being 'RNA' if the vcf is from the rna-seq mutect pipeline
import os
import subprocess
import sys
import shutil


def _call(cmdline, stdout):
    print ''
    print '*' * 70
    print cmdline
    subprocess.call(cmdline.split(), stdout=stdout)


def _call_and_rename(cmdline, fpath):
    output_fpath = fpath + '_output'
    _call(cmdline, open(output_fpath, 'w'))
    if os.path.isfile(fpath):
        os.remove(fpath)
    os.rename(output_fpath, fpath)


def snpsift_annotate(snpsift_jar, db, vcf_fpath):
    cmdline = 'java -jar %s annotate -v %s %s' % (snpsift_jar, db, vcf_fpath)
    _call_and_rename(cmdline, vcf_fpath)


def snpsift_dbnsfp(snpsift_jar, db, vcf_fpath):
    annots = 'SIFT_score,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,' \
             'MutationTaster_score,MutationTaster_pred,MutationAssessor_score,' \
             'MutationAssessor_pred,FATHMM_score,ESP6500_AA_AF,ESP6500_EA_AF,' \
             'Ensembl_geneid,Ensembl_transcriptid'

    cmdline = 'java -jar %s dbnsfp -f %s -v %s %s' % \
              (snpsift_jar, annots, db, vcf_fpath)
    _call_and_rename(cmdline, vcf_fpath)


def snpeff(snpeff_jar, datadir, ref, vcf_fpath):
    cmdline = 'java -Xmx4g -jar %s eff -dataDir %s -cancer ' \
              '-noLog -1 -i vcf -o vcf %s %s' % \
              (snpeff_jar, datadir, ref, vcf_fpath)
    _call_and_rename(cmdline, vcf_fpath)


def rna_editing_sites(db, vcf_fpath):
    cmdline = 'vcfannotate -b %s -k RNA_editing_site %s' % (db, vcf_fpath)
    _call_and_rename(cmdline, vcf_fpath)


def gatk(gatk_jar, ref_path, sample_fpath):
    cmdline = 'java -Xmx2g -jar %s -R %s -T VariantAnnotator ' \
              '-o %s --useAllAnnotations --variant %s' % \
              (gatk_jar, ref_path, sample_fpath, sample_fpath)
    _call_and_rename(cmdline, sample_fpath)


def annotate(sample_fpath, is_rna, is_ensemble,
             ref_name, ref_path, snpeff_dirpath, gatk_dirpath,
             dbsnp_db, cosmic_db, db_nsfp_db, snpeff_datadir,
             annot_track):
    sample_fname = os.path.basename(sample_fpath)
    sample_basepath, ext = os.path.splitext(sample_fpath)
    sample_basename, ext = os.path.splitext(sample_fname)

    # sample_dbsnp_fpath = sample_basepath + '.dbsnp' + ext
    # sample_cosmic_fpath = sample_basepath + '.cosmic' + ext
    # sample_snpeff_fpath = sample_basepath + '.snpeff' + ext
    # sample_dbnsfp_fpath = sample_basepath + '.dbnsf' + ext
    # sample_rna_edit_sites_fpath = sample_basepath + '.RNAeditSites' + ext
    # -------

    if is_ensemble:
        tmp_fpath = sample_basepath + '.tmp' + ext, 'w'
        with open(tmp_fpath) as tmp_f:
            cmdline = 'vcf-subset -c %s -e %s' % \
                      (sample_basename.replace('-ensemble', ''), sample_fpath)
            print cmdline
            subprocess.call(cmdline.split(), stdout=tmp_f)

        os.rename(sample_fpath, sample_basepath + '.combined' + ext)
        os.rename(tmp_fpath, sample_fpath)

    if is_rna:
        pass_sample_fpath = sample_basepath + '.PASS' + ext
        with open(sample_fpath) as sample, \
             open(pass_sample_fpath, 'w') as pass_sample:
            for line in sample.readlines():
                if 'REJECT' not in line:
                    pass_sample.write(line)

        sample_fpath = pass_sample_fpath
        sample_fname = os.path.basename(sample_fpath)
        sample_basename, ext = os.path.split(sample_fname)
        sample_basepath, ext = os.path.split(sample_fpath)

    snpsift_jar = os.path.join(snpeff_dirpath, 'SnpSift.jar')
    snpeff_jar = os.path.join(snpeff_dirpath, 'snpEff.jar')
    gatk_jar = os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

    snpsift_annotate(snpsift_jar, dbsnp_db, sample_fpath)
    snpsift_annotate(snpsift_jar, cosmic_db, sample_fpath)
    snpsift_dbnsfp(snpsift_jar, db_nsfp_db, sample_fpath)
    snpeff(snpeff_jar, snpeff_datadir, ref_name, sample_fpath)
    if is_rna:
        rna_editing_sites(annot_track, sample_fpath)
    gatk(gatk_jar, ref_path, sample_fpath)

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

    print '*' * 70
    print cmdline
    tsv_fpath = sample_basepath + '.tsv'
    subprocess.call(
        cmdline,
        shell=True,
        stdout=open(tsv_fpath, 'w'))


def remove_quotes(str):
    if str and str[0] == '"':
        str = str[1:]
    if str and str[-1] == '"':
        str = str[:-1]
    return str


def split_genotypes(vcf_fpath):
    output = vcf_fpath + '_output'
    with open(vcf_fpath) as vcf, open(output, 'w') as out:
        for line in vcf:
            clean_line = line.strip()
            if not clean_line or clean_line[0] == '#':
                out.write(line)
            else:
                tokens = line.split()
                id_field = remove_quotes(tokens[2])
                alt_field = remove_quotes(tokens[4])

                ids = id_field.split(',')
                alts = alt_field.split(',')
                assert len(ids) == len(alts), 'Number of IDs is not equal to the number of ALTs'
                if len(ids) > 1:
                    for id, alt in zip(ids, alts):
                        line = '\t'.join(tokens[:2] + [id] + [tokens[3]] + [alt] + tokens[5:]) + '\n'
                        out.write(line)
                else:
                    out.write(line)

    os.remove(vcf_fpath)
    os.rename(output, vcf_fpath)
    return sample


if __name__ == '__main__':
    args = sys.argv[1:]

    rna = len(args) > 4 and args[4].lower() == 'rna'
    ensemble = len(args) > 3 and args[3].lower() == 'true'
    do_split_genotypes = len(args) > 2 and args[2].lower() == 'true'
    if len(args) < 2:
        print >> sys.stderr, 'Usage: python filter_snpeff_qsub.py sample.vcf result.vcf [true] [true] [RNA]'
        exit(1)
    result = args[1]
    sample = args[0]

    ref_name = 'hg19'
    ref_path = '/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa'
    snp_eff = '/group/ngs/src/snpEff/snpEff3.5/'
    gatk_dir = '/opt/az/broadinstitute/gatk/1.6'
    dbsnp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/dbsnp_137.vcf'
    cosmic_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cosmic-v67_20131024-hg19.vcf'
    db_nsfp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/dbNSF/dbNSFP2.3/dbNSFP2.3.txt.gz'
    snpeff_datadir = '/ngs/reference_data/genomes/Hsapiens/hg19/snpeff'
    annot_track = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/Human_AG_all_hg19_INFO.bed'

    if result != sample:
        if os.path.exists(result):
            os.remove(result)
        assert os.path.isfile(os.path.realpath(sample)), \
            os.path.realpath(sample) + ' does not exists or is not a file'
        shutil.copyfile(sample, result)
        sample = result

    if do_split_genotypes:
        sample = split_genotypes(sample)
    exit(1)

    print 'Please, run this before start:'
    print '   source /etc/profile.d/modules.sh'
    print '   module load bcbio-nextgen/0.7.6'
    print ''
    print 'In Waltham, run this as well:'
    print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl'

    annotate(sample, rna, ensemble,
             ref_name, ref_path, snp_eff, gatk_dir,
             dbsnp_db, cosmic_db, db_nsfp_db,
             snpeff_datadir, annot_track)