import matplotlib
import matplotlib.ticker as ticker
import matplotlib.pyplot
from os.path import join
from source import verify_file, info
from source.utils import get_chr_lengths

cnv_plot_ending = '.cnv.png'


def draw_seq2c_plot(cnf, seq2c_tsv_fpath, sample_name, output_dir, key_gene_names):
    info('Seq2C plot builder')
    plot_fpath = join(output_dir, sample_name + cnv_plot_ending)
    if cnf.reuse_intermediate and verify_file(plot_fpath, silent=True):
        info('Seq2C plot ' + plot_fpath + ' exists, reusing...')
        return plot_fpath

    if not verify_file(seq2c_tsv_fpath, 'Seq2C.tsv'):
        return None

    chr_names_lengths = get_chr_lengths(cnf)
    chr_names = [chr for chr in chr_names_lengths.keys()]
    chr_short_names = [chr[3:] for chr in chr_names_lengths.keys()]
    chr_lengths = [chr for chr in chr_names_lengths.values()]

    fig = matplotlib.pyplot.figure(figsize=(25, 5))
    matplotlib.pyplot.xlim([0, len(chr_lengths)+2])
    chr_cum_lengths = [sum(chr_lengths[:i]) for i in range(len(chr_lengths)+1)]
    matplotlib.pyplot.xticks(chr_cum_lengths, [])

    ax = matplotlib.pyplot.gca()
    chr_names_coords = [chr_cum_lengths[i+1] - chr_lengths[i]/2 for i in range(len(chr_lengths))]
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_names_coords))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(chr_short_names))

    max_y = 0
    min_y = 0

    with open(seq2c_tsv_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            fs = l.strip().split('\t')
            if fs[0] == sample_name and len(fs) > 12 and fs[1] in key_gene_names:
                sample, gene, chr, start, end, length, log2r, sig, type, amp_del, ab_seg, total_seg, ab_log2r = fs [:13]
                x_vals = [chr_cum_lengths[chr_names.index(chr)] + (int(start) + int(end))/2]
                point_y = float(ab_log2r)
                y_vals = [point_y]
                max_y = max(max_y, point_y)
                min_y = min(min_y, point_y)
                matplotlib.pyplot.plot(x_vals, y_vals, 'ro', markersize=2)
            elif fs[1] in key_gene_names:
                sample, gene, chr, start, end, length, log2r = fs [:7]
                x_vals = [chr_cum_lengths[chr_names.index(chr)] + (int(start) + int(end))/2]
                point_y = float(log2r)
                y_vals = [point_y]
                max_y = max(max_y, point_y)
                min_y = min(min_y, point_y)
                matplotlib.pyplot.plot(x_vals, y_vals, 'bo', markersize=2)

    matplotlib.pyplot.ylim([min_y * 1.05, max_y * 1.05])
    matplotlib.pyplot.tick_params(
        axis='x',
        which='minor',
        bottom='off',
        top='off',
        labelbottom='on')
    info('Saving plot to ' + plot_fpath)
    matplotlib.pyplot.tight_layout()
    fig.savefig(plot_fpath, bbox_inches='tight')
    matplotlib.pyplot.close(fig)

    info('Done')
    info('-' * 70)
    return plot_fpath

