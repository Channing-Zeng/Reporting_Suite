from os.path import join
from source.reporting import summarize, write_summary_reports, Metric, Record
from source.logger import step_greetings, info
from source.bcbio_structure import BCBioStructure


def summary_reports(cnf, bcbio_structure):
    step_greetings('ngsCAT statistics for all samples')

    html_by_sample = bcbio_structure.get_ngscat_html_by_sample()
    sum_report = summarize(html_by_sample, _parse_ngscat_sample_report)

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        sum_report,
        join(cnf.output_dir, BCBioStructure.ngscat_name),
        'ngsCAT statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)

    return final_summary_report_fpaths


METRICS = Metric.to_dict([
    Metric('Number reads',                       'Reads',              'Number of mapped reads'),
    Metric('% target bases with coverage >= 1x', 'Target covered',     '% target bases with coverage >= 1x'),
    Metric('Coverage saturation',                'Saturation',         'Coverage saturation (slope at the end of the curve)',           quality='Less is better'),
    Metric('% reads on target',                  '% reads on target',  '% reads on target'),
    Metric('Duplicated reads on/off target',     'Duplicated reads',   '% duplicated reads on/off target. '
                                                                       'Percentage of duplicated on-target reads normally should be greater '
                                                                       'than the percentage of duplicated off-target reads',           quality='Equal'),
    Metric('mean coverage',                      'Mean cov.',          'Coverage distribution (mean target coverage)'),
    Metric('Coverage per position',              'Cov. per position',  'Coverage per position (consecutive bases with coverage <= 6x)', quality='Less is better'),
    Metric('Standard deviation of coverage',     'Cov. std. dev.',     'Standard deviation of coverage within regions',                 quality='Less is better')
])

ALLOWED_UNITS = ['%']


def _parse_ngscat_sample_report(report_fpath):
    records = []

    def __parse_cell(metric_name, line):
        record = Record(METRICS[metric_name])

        crop_left = line.split('>')
        if len(crop_left) < 2:
            record.value = None
            return
        crop_right = crop_left[1].split('<')
        val = crop_right[0].strip()
        val = val.replace(' ', '').replace(';', '/')

        num_chars = []
        unit_chars = []
        i = 0
        while i < len(val) and (val[i].isdigit() or val[i] in ['.', 'e', '+', '-']):
            num_chars += val[i]
            i += 1
        while i < len(val):
            unit_chars += val[i]
            i += 1
        val_num = ''.join(num_chars)
        val_unit = ''.join(unit_chars)

        if val_unit and val_unit in ALLOWED_UNITS:
            record.metric.unit = val_unit
        try:
            val = int(val_num)
        except ValueError:
            try:
                val = float(val_num)
            except ValueError:  # it is a string
                val = val_num + val_unit
        record.value = val
        return record

    with open(report_fpath) as f:
        parsing_summary_table = False
        header_id = -1
        cell_id = -1
        column_id_to_metric_name = dict()
        for line in f:
            # looking for table's beginning
            if not parsing_summary_table:
                if line.find("table id='summary_table'") != -1:
                    parsing_summary_table = True
                continue
            # checking table's end
            if line.find('Overall status') != -1:
                break

            # parsing header
            if line.find('class="table-header"') != -1:
                header_id += 1
                for metric_name in METRICS.keys():
                    if all(word in line for word in metric_name.split()):
                        column_id_to_metric_name[header_id] = metric_name
                        break
            # parsing content
            if line.find('class="table-cell"') != -1:
                for subline in line.split('class="table-cell"')[1:]:  # for old fashion ngsCAT reports
                    cell_id += 1
                    if cell_id in column_id_to_metric_name.keys():
                        metric_name = column_id_to_metric_name[cell_id]
                        records.append(__parse_cell(metric_name, subline))
    return records