from itertools import repeat, izip
from os.path import join

from source.targetcov.copy_number import run_copy_number


_BCBIO_DIR = "test/bcbio_test_dir/"
_FILE_MAPPED_READS_PATH = _BCBIO_DIR + "reads_mapped.txt"
_FILE_COV = _BCBIO_DIR + "cov.txt"


def summarize_copy_number(report_summary_fpath, fpath_cov):

    gene_summary_lines = _get_lines_by_region_type(report_summary_fpath, 'Whole-Gene')

    cov_by_sample = _get_lines_from_mapped_reads(fpath_cov)

    report_data = run_copy_number(cov_by_sample, gene_summary_lines)

    return report_data


def _get_lines_by_region_type(report_fpath, region_type):
    gene_summary_lines = []

    with open(report_fpath, 'r') as f:
        for line in f:
            if region_type in line:
                #hack
                new_array = []
                tokens = line.split()

                new_array.append(tokens[0])
                new_array.extend(tokens[2:5])
                new_array.append(tokens[1])
                new_array.extend(tokens[5:8])
                gene_summary_lines.append(new_array)

    return gene_summary_lines


def _get_lines_from_mapped_reads(fpath_cov):
    cov_by_sample = dict()
    with open(fpath_cov, 'r') as f:
        it = iter(f)
        next(it)

        for line in it:
            tokens = line.split()
            cov_by_sample[tokens[0]] = float(tokens[1])
    return cov_by_sample


def write_txt(rows, output_dirpath, base_fname):
    output_fpath = join(output_dirpath, base_fname + '.txt')

    col_widths = repeat(0)

    for row in rows:
        col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

    with open(output_fpath, 'w') as out:
        for row in rows:
            for val, w in izip(row, col_widths):
                out.write(val + (' ' * (w - len(val) + 2)))
            out.write('\n')

    return output_fpath


if __name__ == '__main__':
    data = summarize_copy_number(_FILE_COV, _FILE_MAPPED_READS_PATH)
    write_txt(data, _BCBIO_DIR , "z_test")
