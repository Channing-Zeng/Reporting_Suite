#!/usr/bin/env python
from os.path import abspath, dirname, realpath, join, exists, splitext, relpath
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict
import os
import sys
from source import verify_file
from source.logger import err
from source.file_utils import verify_dir, safe_mkdir, adjust_path


class Sample:
    def __init__(self):
        self.name = None
        # self.bam = None
        self.patient = None
        self.is_normal = False
        self.is_blood = False
        self.batches = set()

    def __repr__(self):
        return self.name


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


def main():
    if len(sys.argv) <= 3:
        sys.exit('Usage ' + __file__ + ' input.tsv dir_with_bams dir_with_filtered_bam_symlinks > output.tsv')

    inp_fpath = verify_file(sys.argv[1], is_critical=True)
    bams_dirpath = verify_dir(sys.argv[2], is_critical=True)
    filtered_bams_dirpath = adjust_path(sys.argv[3])
    verify_dir(join(filtered_bams_dirpath, os.pardir), is_critical=True)

    columns_names = 'study	barcode	disease	disease_name	sample_type	sample_type_name	analyte_type	library_type	center	center_name	platform	platform_name	assembly	filename	 files_size 	checksum	analysis_id	aliquot_id	participant_id	sample_id	tss_id	sample_accession	published	uploaded	modified	state	reason'

    samples_by_patient = defaultdict(list)

    with open(inp_fpath) as f:
        for i, l in enumerate(f):
            if not l.strip():
                continue

            fs = l.split('\t')
            if fs[0] == 'study':
                pass
            else:
                barcode = fs[1].split('-')  # TCGA-05-4244-01A-01D-1105-08

                sample = Sample()
                sample.name = os.path.splitext(fs[13])[0]
                sample.description = fs[1]
                sample.patient = '-'.join(barcode[:3])
                sample.reason = fs[26]

                sample_type = int(barcode[3][:2])
                if sample_type >= 20 or sample_type <= 0:
                    continue
                sample.is_normal = 10 <= sample_type < 20
                sample.is_blood = sample_type in [3, 4, 9, 10]  # https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm

                if sample.description == 'TCGA-64-1676-10A-01D-0969-08':
                    pass
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
                    prev_v = get_bam_version(prev_sample.name)
                    v = get_bam_version(sample.name)
                    err('Duplicated sample: ' + sample.description + '  Resolving by version (' + ' over '.join(map(str, sorted([prev_v, v])[::-1])) + ')')
                    if v > prev_v:
                        samples_by_patient[sample.patient].remove(prev_sample)
                        samples_by_patient[sample.patient].append(sample)
                else:
                    samples_by_patient[sample.patient].append(sample)

    batches = []
    final_samples = set()

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


    safe_mkdir(filtered_bams_dirpath)

    for fname in os.listdir(bams_dirpath):
        if fname.endswith('.bam'):
            basefname = fname.split('.bam')[0]

            if basefname in [s.name for s in final_samples]:
                src_fpath = join(bams_dirpath, fname)
                dst_fpath = join(filtered_bams_dirpath, fname)
                if not exists(dst_fpath):
                    err('Symlinking ' + src_fpath + ' -> ' + dst_fpath)
                    os.symlink(src_fpath, dst_fpath)
                if not exists(dst_fpath + '.bai'):
                    os.symlink(src_fpath + '.bai', dst_fpath + '.bai')

    print 'sample,description,batch,phenotype'
    for s in sorted(final_samples, key=lambda s: s.name):
        print ','.join([s.name, s.description, ';'.join(sorted(b.name for b in s.batches)), ('normal' if s.is_normal else 'tumor')])


if __name__ == '__main__':
    main()