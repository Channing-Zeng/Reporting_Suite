#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from os.path import abspath, dirname, realpath, join, exists, splitext, relpath, basename
from collections import defaultdict
import os
import sys
import traceback
import pysam
from source.targetcov.bam_and_bed_utils import verify_bam
from source import verify_file
from source.logger import err
from source.file_utils import verify_dir, safe_mkdir, adjust_path


class Sample:
    def __init__(self):
        self.bam_base_name = None
        self.bam = None
        self.description = None
        self.patient = None
        self.is_normal = False
        self.is_blood = False
        self.batches = set()

    def __repr__(self):
        return self.description


class Batch:
    def __init__(self, name=None):
        self.name = name
        self.normal = None
        self.tumour = None


def get_bam_version(sample_name):
    parts = sample_name.split('.')
    if len(parts) > 1:
        try:
            version = int(parts[-1])
        except ValueError:
            return None
        else:
            return version


def main(args):
    if len(args) < 2:
        sys.exit('Usage ' + __file__ + ' input.tsv bcbio.csv [dir_with_bams] [bina_dir]')

    inp_fpath = args[0]
    verify_file(args[0], is_critical=True)

    out_fpath = args[1]
    verify_dir(dirname(adjust_path(out_fpath)), is_critical=True)

    bam_dirpath = None
    if len(args) > 2:
        bam_dirpath = args[2]
        verify_dir(adjust_path(bam_dirpath), is_critical=True)

    # bam_opt = args[2]
    # try:
    #     bam_col = int(bam_opt)
    #     bam_dirpath = None
    # except ValueError:
    #     bam_col = None
    #     verify_dir(bam_opt, is_critical=True)
    #     bam_dirpath = args[2]

    bina_dirpath = None
    if len(args) > 3:
        bina_dirpath = args[3]
        verify_dir(dirname(adjust_path(bina_dirpath)), is_critical=True)

    # filtered_bams_dirpath = adjust_path(sys.argv[3])
    # verify_dir(join(filtered_bams_dirpath, os.pardir), is_critical=True)

    columns_names = 'study	barcode	disease	disease_name	sample_type	sample_type_name	analyte_type	library_type	center	center_name	platform	platform_name	assembly	filename	 files_size 	checksum	analysis_id	aliquot_id	participant_id	sample_id	tss_id	sample_accession	published	uploaded	modified	state	reason'

    samples_by_patient = defaultdict(list)

    delim = '\t'
    barcode_col = 1
    bam_col = 13
    is_tcga_tsv = True

    with open(inp_fpath) as fh:
        for i, l in enumerate(fh):
            if not l.strip():
                continue

            if i == 0:
                if len(l.split('\t')) == 27:
                    err('Interpreting as TCGA tsv')
                    if l.split('\t')[0] != 'TCGA': continue  # skipping header
                else:
                    delim = None
                    for j, f in enumerate(l.split()):
                        if f.startswith('TCGA'):
                            barcode_col = j
                            err('barcode col is ' + str(j))
                        if f.endswith('bam'):
                            bam_col = j
                            err('bam col is ' + str(j))
                    is_tcga_tsv = False

            fs = l.split(delim)

            barcode = fs[barcode_col].split('-')  # TCGA-05-4244-01A-01D-1105-08

            sample = Sample()
            sample.bam = fs[bam_col]
            sample.bam_base_name = basename(os.path.splitext(fs[bam_col])[0])
            sample.description = fs[barcode_col]
            sample.patient = '-'.join(barcode[:3])
            if is_tcga_tsv:
                sample.reason = fs[26]

            sample_type = int(barcode[3][:2])
            if sample_type >= 20 or sample_type <= 0:
                continue
            sample.is_normal = 10 <= sample_type < 20
            sample.is_blood = sample_type in [3, 4, 9, 10]  # https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm

            if any(s.description == sample.description for s in samples_by_patient[sample.patient]):
                prev_sample = next(s for s in samples_by_patient[sample.patient] if s.description == sample.description)

                # comp reason
                # if 'Fileset modified' not in prev_sample.reason and 'Fileset modified' in sample.reason:
                #     err('Duplicated sample: ' + sample.description + '  Fileset modified not in old ' + prev_sample.name + ' over ' + sample.name)
                #     pass
                # elif 'Fileset modified' in prev_sample.reason and 'Fileset modified' not in sample.reason:
                #     samples_by_patient[sample.patient].remove(prev_sample)
                #     samples_by_patient[sample.patient].append(sample)
                #     err('Duplicated sample: ' + sample.description + '  Fileset modified not in new ' + sample.name + ' over ' + prev_sample.name)
                # else:
                # comp version
                prev_version = get_bam_version(prev_sample.bam_base_name)
                version = get_bam_version(sample.bam_base_name)
                err('Duplicated sample: ' + sample.description + '  Resolving by version (' + ' over '.join(map(str, sorted([prev_version, version])[::-1])) + ')')
                if version > prev_version:
                    samples_by_patient[sample.patient].remove(prev_sample)
                    samples_by_patient[sample.patient].append(sample)
            else:
                samples_by_patient[sample.patient].append(sample)

    batches = []
    final_samples = set()

    if bina_dirpath:
        safe_mkdir(bina_dirpath)

    for patient, patient_samples in samples_by_patient.iteritems():
        tumours = [s for s in patient_samples if not s.is_normal]
        normals = [s for s in patient_samples if s.is_normal]

        main_normal = None
        if len(normals) >= 1:
            if any(n.is_blood for n in normals):
                main_normal = next(n for n in normals if n.is_blood)
            else:
                main_normal = normals[0]
                if tumours:
                    for n in normals[1:]:
                        b = Batch(n.description + '-batch')
                        b.tumour = n
                        batches.append(b)

        for t in tumours:
            b = Batch(t.description + '-batch')
            b.tumour = t
            t.batches.add(b)
            final_samples.add(t)
            if main_normal:
                b.normal = main_normal
                main_normal.batches.add(b)
                final_samples.add(main_normal)
            batches.append(b)

        ##################
        ###### Bina ######
        if bina_dirpath:
            bina_patient_dirpath = join(bina_dirpath, patient)
            safe_mkdir(bina_patient_dirpath)
            normals_csv_fpath = join(bina_patient_dirpath, 'normals.csv')
            tumours_csv_fpath = join(bina_patient_dirpath, 'tumors.csv')

            if main_normal:
                with open(normals_csv_fpath, 'w') as f:
                    f.write('name,bam\n')
                    bam_fpath = join(bam_dirpath, main_normal.bam) if bam_dirpath else main_normal.bam
                    f.write(main_normal.description + ',' + bam_fpath + '\n')

            with open(tumours_csv_fpath, 'w') as f:
                f.write('name,bam\n')
                for t in tumours:
                    bam_fpath = join(bam_dirpath, t.bam) if bam_dirpath else t.bam
                    f.write(t.description + ',' + bam_fpath + '\n')

    if bina_dirpath:
        err('Saved bina CSVs to ' + bina_dirpath)

    ###########################
    ######## Bcbio CSV ########
    print 'bcbio_nextgen.py -w template bcbio.yaml', out_fpath,
    with open(out_fpath, 'w') as out:
        out.write('sample,description,batch,phenotype\n')
        for s in sorted(final_samples, key=lambda s: s.bam_base_name):
            out.write(','.join([s.bam_base_name, s.description, ';'.join(sorted(b.name for b in s.batches)),
                ('normal' if s.is_normal else 'tumor')]) + '\n')
            bam_fpath = join(bam_dirpath, s.bam) if bam_dirpath else s.bam

            if verify_bam(bam_fpath, is_critical=False):
                try:
                    bam = pysam.Samfile(bam_fpath, "rb")
                except ValueError:
                    err(traceback.format_exc())
                    err('Cannot read ' + bam_fpath)
                    err()
                    # n_rgs = max(1, len(bam.header.get("RG", [])))
                else:
                    print bam_fpath,

    ############################
    ###### Fixed BAM list ######
    # safe_mkdir(filtered_bams_dirpath)
    # for fname in os.listdir(bams_dirpath):
    #     if fname.endswith('.bam'):
    #         basefname = fname.split('.bam')[0]
    #
    #         if basefname in [s.name for s in final_samples]:
    #             src_fpath = join(bams_dirpath, fname)
    #             dst_fpath = join(filtered_bams_dirpath, fname)
    #             if not exists(dst_fpath):
    #                 err('Symlinking ' + src_fpath + ' -> ' + dst_fpath)
    #                 os.symlink(src_fpath, dst_fpath)
    #             if not exists(dst_fpath + '.bai'):
    #                 os.symlink(src_fpath + '.bai', dst_fpath + '.bai')
    #
    # err('Saved symlinks to BAMs to ' + filtered_bams_dirpath)


if __name__ == '__main__':
    main(sys.argv[1:])