database = 'cosmic'
main_novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
sample_header = 'Sample name:'


def parse_qc_sample_report(report_fpath):
    """ returns row_per_sample =
            dict(metricName=None, value=None,
            isMain=True, quality='More is better')
    """
    row_per_sample = []

    with open(report_fpath) as f:
        # parsing Sample name and Database columns
        database_col_id = None
        novelty_col_id = None
        for line in f:
            if line.startswith(sample_header):
                sample_name = line[len(sample_header):].strip()
            elif line.startswith(metrics_header):
                if database in line:
                    database_col_id = line.split().index(database)
                if novelty_header in line:
                    novelty_col_id = line.split().index(novelty_header)
                break

        if database_col_id:
            # parsing rest of the report
            for line in f:
                is_main = True
                novelty = line.split()[novelty_col_id]
                if novelty_col_id and novelty != main_novelty:
                    is_main = False

                cur_metric_name = line.split()[0] + ' ' + novelty
                cur_value = line.split()[database_col_id]

                row_per_sample.append(dict(
                    metricName=cur_metric_name, value=cur_value,
                    isMain=is_main, quality='More is better'))

    return row_per_sample