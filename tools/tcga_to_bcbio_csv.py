#!/usr/bin/env python
from os.path import abspath, dirname, realpath, join, exists
from site import addsitedir
project_dir = abspath(dirname(dirname(dirname(realpath(__file__)))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict
import os
import sys
from source import verify_file


class Sample:
    def __init__(self):
        self.name = None
        self.bam = None
        self.patient = None
        self.is_normal = False
        self.is_blood = False
        self.batch = None

    def __repr__(self):
        return self.name


class Batch:
    def __init__(self, name=None):
        self.name = name
        self.normal = None
        self.tumour = None


def main():
    inp = sys.stdin
    if len(sys.argv) > 1:
        inp_fpath = sys.argv[1]
        verify_file(inp_fpath, is_critical=True)
        inp = open(inp_fpath)

    columns_names = 'study	barcode	disease	disease_name	sample_type	sample_type_name	analyte_type	library_type	center	center_name	platform	platform_name	assembly	filename	 files_size 	checksum	analysis_id	aliquot_id	participant_id	sample_id	tss_id	sample_accession	published	uploaded	modified	state	reason'

    samples_by_patient = defaultdict(list)

    for i, l in enumerate(inp):
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

            sample_type = int(barcode[3][:2])
            if sample_type >= 20 or sample_type <= 0:
                continue
            sample.is_normal = 10 <= sample_type < 20
            sample.is_blood = sample_type in [3, 4, 9, 10]

            samples_by_patient[sample.patient].append(sample)

    batches = []

    for patient, samples in samples_by_patient.iteritems():
        tumours = [s for s in samples if not s.is_normal]
        normals = [s for s in samples if s.is_normal]

        main_normal = None
        if len(normals) >= 1:
            if any(n.is_blood for n in normals):
                main_normal = next(n for n in normals if n.is_blood)
            else:
                main_normal = normals[0]

        for t in tumours:
            b = Batch(t.description + '-batch')
            b.tumour = t
            if main_normal:
                b.normal = main_normal
            batches.append(b)

    print 'sample\tdescription\tbatch\tphenotype'
    for b in batches:
        print b.tumour.name + '\t' + b.tumour.description + '\t' + b.name + '\ttumor'
        if b.normal:
            print b.normal.name + '\t' + b.normal.description + '\t' + b.name + '\tnormal'


if __name__ == '__main__':
    main()