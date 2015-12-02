from collections import OrderedDict, defaultdict
from itertools import dropwhile
import re
import os
from os.path import join, isfile, isdir, basename
import shutil
from source import TargQC_Sample
from source.logger import critical, err, info, warn
from source.file_utils import verify_dir, verify_file, splitext_plus, safe_mkdir, file_transaction


class DatasetStructure:
    pre_fastqc_repr =        'Preproc FastQC'
    downsample_targqc_repr = 'TargQC downsampled'

    @staticmethod
    def create(dir_path, project_name):
        if 'datasets/miseq/' in dir_path.lower():
            return MiSeqStructure(dir_path, project_name)

        elif 'datasets/hiseq/' in dir_path.lower():
            return HiSeqStructure(dir_path, project_name)

        elif 'datasets/hiseq4000/' in dir_path.lower():
            return HiSeq4000Structure(dir_path, project_name)

        else:
            critical('Directory must be datasets/miseq/, datasets/hiseq/, or datasets/hiseq4000/. Found ' + dir_path)

    def __init__(self, dirpath, az_project_name):
        self.az_project_name = az_project_name

        illumina_project_name = None
        if '/Unalign/' in dirpath:
            self.dirpath = dirpath.split('/Unalign/')[0]
            self.unaligned_dirpath = self.__find_unaligned_dir()
            verify_dir(self.unaligned_dirpath, description='Unalign dir', is_critical=True)
            illumina_project_name = dirpath.split('/Unalign/')[1]  # something like AURA.FFPE.AZ300, in contast with project_name which is something like Bio_123_AURA_FFPE_AZ300
        else:
            self.dirpath = dirpath
            self.unaligned_dirpath = self.__find_unaligned_dir()

        self.basecalls_dirpath = join(self.dirpath, 'Data/Intensities/BaseCalls')
        verify_dir(self.basecalls_dirpath, is_critical=True)

        self.bcl2fastq_dirpath = None
        self.source_fastq_dirpath = None

        self.samplesheet_fpath = self.__find_sample_sheet()
        self.project_by_name = self._parse_sample_sheet(self.samplesheet_fpath)

        if illumina_project_name:  # we want only a specific project
            if illumina_project_name not in self.project_by_name:
                info()
                critical('ERROR: Project ' + illumina_project_name + ' not in the SampleSheet ' + self.samplesheet_fpath)
            self.project_by_name = {illumina_project_name: self.project_by_name[illumina_project_name]}

    def __find_unaligned_dir(self):
        unaligned_dirpath = join(self.dirpath, 'Unalign')
        if verify_dir(unaligned_dirpath, description='"Unalign" directory', silent=True):
            unaligned_dirpath = unaligned_dirpath
        else:
            unaligned_dirpath = None
            warn('No unalign directory')
        return unaligned_dirpath

    def __find_sample_sheet(self):
        ss_fpath = join(self.dirpath, 'SampleSheet.csv')
        if not isfile(ss_fpath):
            ss_fpath = join(self.basecalls_dirpath, 'SampleSheet.csv')
        verify_file(ss_fpath, description='Sample sheet', is_critical=True)
        return ss_fpath

    def _parse_sample_sheet(self, sample_sheet_fpath):
        info('Parsing sample sheet ' + sample_sheet_fpath)
        with open(sample_sheet_fpath) as f:
            def check_if_header(l):
                return any(l.startswith(w) for w in [
                    'Sample_ID,',  # MiSeq
                    'FCID,',       # HiSeq
                    'Lane,',       # HiSeq4000
                ])

            sample_lines = dropwhile(lambda l: not check_if_header(l), f)
            sample_infos = []
            keys = []
            for l in sample_lines:
                if check_if_header(l):
                    keys = l.strip().split(',')
                else:
                    fs = l.strip().split(',')
                    sample_infos.append(dict(zip(keys, fs)))

        project_by_name = OrderedDict()

        for i, info_d in enumerate(sample_infos):
            proj_name = info_d.get('Sample_Project', info_d.get('SampleProject'))
            if not proj_name:
                critical('No SampleProject or Sample_Project field in the SampleSheet ' + sample_sheet_fpath)
            if proj_name not in project_by_name:
                project_by_name[proj_name] = DatasetProject(proj_name)
            project = project_by_name[proj_name]

            sname = info_d.get('Sample_ID', info_d.get('SampleID'))
            if sname in project.sample_by_name:
                s = project.sample_by_name[sname]
                s.lane_numbers.add(info_d.get('Lane', 1))  # lanes are in HiSeq and HiSeq4000 (not in MiSeq!)
            else:
                s = DatasetSample(sname, index=info_d.get('index', info_d.get('Index')))
                s.lane_numbers.add(info_d.get('Lane', 1))  # lanes are in HiSeq and HiSeq4000 (not in MiSeq!)
                if 'FCID' in info_d:
                    s.fcid = info_d['FCID']  # HiSeq only
                info('Sample ' + sname)
                project.sample_by_name[sname] = s
                # sample_names.append(info_d[key].replace(' ', '-') + '_' + info_d['Index'] + '_L%03d' % lane)
                # sample_names.append(info_d[key].replace(' ', '-').replace('_', '-') + '_S' + str(i + 1) + '_L001')

        return project_by_name


class HiSeqStructure(DatasetStructure):
    def __init__(self, dirpath, az_project_name=None):
        info('Parsing the HiSeq project structure')
        self.kind = 'hiseq'
        DatasetStructure.__init__(self, dirpath, az_project_name)

        verify_dir(self.unaligned_dirpath, is_critical=True)

        for pname, project in self.project_by_name.items():
            proj_dirpath = join(self.unaligned_dirpath, 'Project_' + pname.replace(' ', '-'))  #.replace('-', '_').replace('.', '_'))
            project.set_dirpath(proj_dirpath, self.az_project_name)
            for sname, sample in project.sample_by_name.items():
                sample.source_fastq_dirpath = join(project.dirpath, 'Sample_' + sname.replace(' ', '-'))  #.replace('-', '_').replace('.', '_'))
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

        self.basecall_stat_html_reports = self.__get_basecall_stats_reports()

    # def __get_bcl2fastq_dirpath(self):
    #     # Reading project name
    #     bcl2fastq_dirpath = None
    #     try:
    #         bcl2fastq_dirpath = join(self.unaligned_dirpath, next(fn for fn in os.listdir(self.unaligned_dirpath) if fn.startswith('Project_')))
    #     except StopIteration:
    #         critical('Could not find directory starting with Project_ in ' + self.unaligned_dirpath)
    #     return bcl2fastq_dirpath

    def __get_basecall_stats_reports(self):
        basecall_stats_dirnames = [fname for fname in os.listdir(self.unaligned_dirpath) if fname.startswith('Basecall_Stats_')]
        if len(basecall_stats_dirnames) > 1:
            err('More than 1 Basecall_Stats_* dirs found in unalign_dirpath')
        if len(basecall_stats_dirnames) == 0:
            err('No Basecall_Stats_* dirs found in unalign_dirpath')
        if len(basecall_stats_dirnames) == 1:
            basecall_stats_dirpath = join(self.unaligned_dirpath, basecall_stats_dirnames[0])
            basecall_reports = [verify_file(join(basecall_stats_dirpath, html_fname)) for html_fname in ['Demultiplex_Stats.htm', 'All.htm', 'IVC.htm']]
            return filter(None, basecall_reports)


class MiSeqStructure(DatasetStructure):
    def __init__(self, dirpath, az_project_name):
        info('Parsing the MiSeq project structure')
        self.kind = 'miseq'
        DatasetStructure.__init__(self, dirpath, az_project_name)

        base_dirpath = self.unaligned_dirpath
        if not verify_dir(base_dirpath, silent=True):
            base_dirpath = self.basecalls_dirpath
        verify_dir(base_dirpath, description='Source fastq dir')

        for proj_name, project in self.project_by_name.items():
            project.set_dirpath(base_dirpath, self.az_project_name)
            for sample in project.sample_by_name.values():
                sample.source_fastq_dirpath = project.dirpath
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

        self.basecall_stat_html_reports = []

    def __find_fastq_dir(self):
        for dname in os.listdir(self.unaligned_dirpath):
            dpath = join(self.unaligned_dirpath, dname)
            if isdir(dpath) and any(f.endswith('.fastq.gz') for f in os.listdir(dpath)):
                return dpath

class HiSeq4000Structure(DatasetStructure):
    def __init__(self, dirpath, az_project_name):
        info('Parsing the HiSeq4000 project structure - same as MiSeq')
        self.kind = 'hiseq4000'
        DatasetStructure.__init__(self, dirpath, az_project_name)

        verify_dir(self.unaligned_dirpath, is_critical=True)

        for proj_name, project in self.project_by_name.items():
            proj_dirpath = join(self.unaligned_dirpath, proj_name)
            az_project_name = self.az_project_name
            if len(self.project_by_name) > 1:
                az_project_name += '_' + proj_name.replace(' ', '_').replace('-', '_').replace('.', '_')
            project.set_dirpath(proj_dirpath, az_project_name)
            for sample in project.sample_by_name.values():
                sample.source_fastq_dirpath = project.dirpath
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

        self.basecall_stat_html_reports = self.__get_basecall_stats_reports()

    def __find_fastq_dir(self):
        for dname in os.listdir(self.unaligned_dirpath):
            dpath = join(self.unaligned_dirpath, dname)
            if isdir(dpath) and any(f.endswith('.fastq.gz') for f in os.listdir(dpath)):
                return dpath

    def __get_basecall_stats_reports(self):
        dirpath = join(self.unaligned_dirpath, 'Reports', 'html')
        index_html_fpath = join(dirpath, 'index.html')
        if verify_dir(dirpath) and verify_file(index_html_fpath):
            return [index_html_fpath]

class DatasetProject:
    def __init__(self, name):
        self.name = name
        self.sample_by_name = OrderedDict()
        self.dirpath = None

    def set_dirpath(self, dirpath, az_project_name):
        self.dirpath = dirpath
        self.az_project_name = az_project_name
        verify_dir(self.dirpath, is_critical=True)

        merged_dirpath = join(self.dirpath, 'merged')
        if verify_dir(merged_dirpath, silent=True):
            self.mergred_dir_found = True
            self.fastq_dirpath = self.fastqc_dirpath = merged_dirpath
        else:
            self.mergred_dir_found = False
            self.fastq_dirpath = join(self.dirpath, 'fastq')
            self.fastqc_dirpath = join(self.fastq_dirpath, 'FastQC')
        info()

        self.comb_fastqc_fpath = join(self.fastqc_dirpath, 'FastQC.html')
        self.downsample_targqc_report_fpath = None
        self.project_report_html_fpath = None

        self.downsample_metamapping_dirpath = join(self.dirpath, 'Downsample_MetaMapping')
        self.downsample_targqc_dirpath = join(self.dirpath, 'Downsample_TargQC')
        self.downsample_targqc_report_fpath = join(self.downsample_targqc_dirpath, 'targQC.html')
        self.project_report_html_fpath = join(self.dirpath, az_project_name + '.html')

    def concat_fastqs(self, cnf):
        info('Concatenating fastqc files for ' + self.name)
        if self.mergred_dir_found:
            info('  found already merged fastq dir, skipping.')
            return
        if not self.sample_by_name:
            err('  no samples found.')
            return
        safe_mkdir(self.fastq_dirpath)
        for s in self.sample_by_name.values():
            _concat_fastq(cnf, s.find_raw_fastq('R1'), s.l_fpath)
            _concat_fastq(cnf, s.find_raw_fastq('R2'), s.r_fpath)
        info()

class DatasetSample:
    def __init__(self, name, index=None, source_fastq_dirpath=None):
        self.name = name
        self.index = None
        self.lane_numbers = set()
        self.fcid = None  # for HiSeq

        self.source_fastq_dirpath = source_fastq_dirpath

    def set_up_out_dirs(self, fastq_dirpath, fastqc_dirpath, downsample_targqc_dirpath):
        self.fastq_dirpath = fastq_dirpath
        self.fastqc_dirpath = fastqc_dirpath
        self.downsample_targqc_dirpath = downsample_targqc_dirpath

        self.l_fpath = join(fastq_dirpath, self.name + '_R1.fastq.gz')
        self.r_fpath = join(fastq_dirpath, self.name + '_R2.fastq.gz')

        self.sample_fastqc_dirpath = join(fastqc_dirpath, self.name + '.fq_fastqc')
        self.fastqc_html_fpath = join(fastqc_dirpath, self.name + '.fq_fastqc.html')
        self.l_fastqc_base_name = splitext_plus(basename(self.l_fpath))[0]
        self.r_fastqc_base_name = splitext_plus(basename(self.r_fpath))[0]
        # self.l_fastqc_html_fpath = None  # join(ds.fastqc_dirpath,  + '_fastqc.html')
        # self.r_fastqc_html_fpath = None  # join(ds.fastqc_dirpath, splitext_plus(self.r_fpath)[0] + '_fastqc.html')

        if not isfile(self.fastqc_html_fpath):
            self.fastqc_html_fpath = join(self.sample_fastqc_dirpath, 'fastqc_report.html')

        self.targqc_sample = TargQC_Sample(self.name, join(downsample_targqc_dirpath, self.name))
        self.targetcov_html_fpath = self.targqc_sample.targetcov_html_fpath
        self.ngscat_html_fpath    = self.targqc_sample.ngscat_html_fpath
        self.qualimap_html_fpath  = self.targqc_sample.qualimap_html_fpath

    def find_raw_fastq(self, suf='R1'):
        fastq_fpaths = [
            join(self.source_fastq_dirpath, fname)
                for fname in os.listdir(self.source_fastq_dirpath)
                if re.match(self.name.replace('-', '.').replace('_', '.').replace(' ', '.') +
                            '.*_' + suf + '.*\.fastq\.gz', fname)]
        fastq_fpaths = sorted(fastq_fpaths)
        if not fastq_fpaths:
            critical('Error: no fastq files for the sample ' + self.name +
                     ' were found inside ' + self.source_fastq_dirpath)
        info(self.name + ': found raw fastq files ' + ', '.join(fastq_fpaths))
        return fastq_fpaths

    def find_fastqc_html(self, end_name):
        sample_fastqc_dirpath = join(self.fastqc_dirpath, end_name + '_fastqc')
        fastqc_html_fpath = join(self.fastqc_dirpath, end_name + '_fastqc.html')
        if isfile(fastqc_html_fpath):
            return fastqc_html_fpath
        else:
            fastqc_html_fpath = join(sample_fastqc_dirpath, 'fastqc_report.html')
            if isfile(fastqc_html_fpath):
                return fastqc_html_fpath
            else:
                return None


def _concat_fastq(cnf, fastq_fpaths, output_fpath):
    info('  merging ' + ', '.join(fastq_fpaths))
    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info(output_fpath + ' exists, reusing')
    else:
        with file_transaction(cnf.work_dir, output_fpath) as tx:
            with open(tx, 'w') as out:
                for fq_fpath in fastq_fpaths:
                    with open(fq_fpath, 'r') as inp:
                        shutil.copyfileobj(inp, out)
    return output_fpath
