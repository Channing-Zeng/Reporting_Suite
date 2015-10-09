from collections import OrderedDict
import matplotlib
from source.logger import warn
matplotlib.use('Agg')
import matplotlib.ticker as ticker
import matplotlib.pyplot
from os.path import join
from source import verify_file, info
from source.utils import get_chr_lengths

cnv_plot_ending = '.cnv.png'


def draw_seq2c_plot(cnf, seq2c_tsv_fpath, sample_name, output_dir, key_gene_names=None):
    info('Seq2C plot builder')
    plot_fpath = join(output_dir, sample_name + cnv_plot_ending)
    if cnf.reuse_intermediate and verify_file(plot_fpath, silent=True):
        info('Seq2C plot ' + plot_fpath + ' exists, reusing...')
        return plot_fpath

    if not verify_file(seq2c_tsv_fpath, 'Seq2C.tsv'):
        return None

    chr_names_lengths = OrderedDict((chr_, l) for chr_, l in get_chr_lengths(cnf).items() if '_' not in chr_) # not drawing extra chromosomes chr1_blablabla
    chr_names = chr_names_lengths.keys()
    chr_short_names = [chr_[3:] for chr_ in chr_names_lengths.keys()]
    chr_lengths = [chr_ for chr_ in chr_names_lengths.values()]

    fig = matplotlib.pyplot.figure(figsize=(25, 5))
    matplotlib.pyplot.xlim([0, len(chr_lengths) + 2])
    chr_cum_lengths = [sum(chr_lengths[:i]) for i in range(len(chr_lengths)+1)]
    matplotlib.pyplot.xticks(chr_cum_lengths, [])

    ax = matplotlib.pyplot.gca()
    chr_names_coords = [chr_cum_lengths[i+1] - chr_lengths[i]/2 for i in range(len(chr_lengths))]
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_names_coords))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(chr_short_names))

    max_y = 0
    min_y = 0
    def add_rec_to_plot(chr_, start, end, log2r, marker, max_y, min_y):
        x_vals = [chr_cum_lengths[chr_names.index(chr_)] + (int(start) + int(end))/2]
        point_y = float(log2r)
        y_vals = [point_y]
        max_y = max(max_y, point_y)
        min_y = min(min_y, point_y)
        matplotlib.pyplot.plot(x_vals, y_vals, marker, markersize=2)
        return max_y, min_y

    with open(seq2c_tsv_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            fs = l.strip().split('\t')
            sname, gname = fs[0], fs[1]
            if key_gene_names and gname not in key_gene_names:
                continue
            if sname != sample_name:
                continue

            if len(fs) > 12:
                marker = 'o'
                sname, gname, chr_, start, end, length, log2r, sig, type_, amp_del, ab_seg, total_seg, ab_log2r = fs[:13]
                if ab_log2r:
                    log2r = ab_log2r
                    marker = 'x'
            else:
                sname, gname, chr_, start, end, length, log2r = fs[:7]
                marker = '.'

            log2r = float(log2r)
            if -2 <= log2r <= 2:
                color = 'b'
            else:
                color = 'r'
            max_y, min_y = add_rec_to_plot(chr_, start, end, log2r, color + marker, max_y, min_y)

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