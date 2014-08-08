from source.reporting import get_sample_report_fpaths_for_bcbio_final_dir, \
    summarize, write_summary_reports, Metric, Record
from source.logger import step_greetings, info
from source.utils import OrderedDefaultDict

def summary_reports(cnf, sample_names):
    step_greetings('ngsCAT statistics for all samples')

    sample_sum_reports, sample_names = get_sample_report_fpaths_for_bcbio_final_dir(
        cnf['bcbio_final_dir'], sample_names, cnf['base_name'], 'captureQC.html', raw_ending=True)

    sum_report = summarize(sample_names, sample_sum_reports, _parse_ngscat_sample_report)

    sum_report_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], sum_report, 'ngscat', 'ngsCAT statistics')

    return sum_report_fpaths


def __to_dict(metrics):
    return {m.name: m for m in metrics}

METRICS = __to_dict([
    Metric('Number reads',                       'reads',              'Number of mapped reads'),
    Metric('% target bases with coverage >= 1x', 'target covered',     '% target bases with coverage >= 1x'),
    Metric('Coverage saturation',                'saturation',         'Coverage saturation (slope at the end of the curve)',           quality='Less is better'),
    Metric('% reads on target',                  '% reads on target',  '% reads on target'),
    Metric('Duplicated reads on/off target',     'duplicated reads',   '% duplicated reads on/off target'),
    Metric('mean coverage',                      'mean cov.',          'Coverage distribution (mean target coverage)'),
    Metric('Coverage per position',              'cov. per position',  'Coverage per position (consecutive bases with coverage <= 6x)', quality='Less is better'),
    Metric('Standard deviation of coverage',     'cov. std. dev.',     'Standard deviation of coverage within regions',                 quality='Less is better')
])

ALLOWED_UNITS = ['%']


def _parse_ngscat_sample_report(report_fpath):
    records = OrderedDefaultDict(Record)

    def __parse_cell(metric_name, line):
        records[metric_name].metric = METRICS[metric_name]

        crop_left = line.split('>')
        if len(crop_left) < 2:
            records[metric_name].value = None
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
            records[metric_name].metric.unit = val_unit
        try:
            val = int(val_num)
        except ValueError:
            try:
                val = float(val_num)
            except ValueError: # it is a string
                val = val_num + val_unit
        records[metric_name].value = val

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
                test = METRICS.keys()
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
                        __parse_cell(metric_name, subline)

    info("report_fpath is " + report_fpath)
    return records