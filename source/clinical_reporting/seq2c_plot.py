from collections import OrderedDict
from itertools import chain
from source.logger import warn

import matplotlib
matplotlib.use('Agg')
import matplotlib.ticker as ticker
import matplotlib.pyplot
from os.path import join
from source import verify_file, info
from source.utils import get_chr_lengths

cnv_plot_ending = '.seq2c.png'


def draw_seq2c_plot(cnf, seq2c_tsv_fpath, sample_name, output_dir, key_gene_names=None, chr_lens=None):
    info('Seq2C plot builder')
    plot_fpath = join(output_dir, sample_name + cnv_plot_ending)
    if cnf.reuse_intermediate and verify_file(plot_fpath, silent=True):
        info('Seq2C plot ' + plot_fpath + ' exists, reusing...')
        return plot_fpath

    if not verify_file(seq2c_tsv_fpath, 'Seq2C.tsv'):
        return None

    chr_names_lengths = OrderedDict((chr_, l) for chr_, l in (chr_lens or get_chr_lengths(cnf))
                                    if '_' not in chr_)  # not drawing extra chromosomes chr1_blablabla
    chr_names = chr_names_lengths.keys()
    chr_short_names = [chrom[3:] for chrom in chr_names_lengths.keys()]
    chr_lengths = [chrom for chrom in chr_names_lengths.values()]

    fig = matplotlib.pyplot.figure(figsize=(25, 5))
    matplotlib.pyplot.xlim([0, len(chr_lengths) + 1])
    chr_cum_lens = [sum(chr_lengths[:i]) for i in range(len(chr_lengths)+1)]
    matplotlib.pyplot.xticks(chr_cum_lens, [])

    ax = matplotlib.pyplot.gca()
    chr_names_coords = [chr_cum_lens[i+1] - chr_lengths[i]/2 for i in range(len(chr_lengths))]
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_names_coords))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(chr_short_names))

    # def add_rec_to_plot(chr_, start, end, log2r, max_y, min_y, marker, color, label=None):
    #     x_vals = [chr_cum_lengths[chr_names.index(chr_)] + (int(start) + int(end))/2]
    #     point_y = float(log2r)
    #     y_vals = [point_y]
    #     max_y = max(max_y, point_y)
    #     min_y = min(min_y, point_y)
    #     if label:
    #         matplotlib.pyplot.plot(x_vals, y_vals, marker, markersize=2, label=label)
    #     else:
    #         matplotlib.pyplot.plot(x_vals, y_vals, marker, markersize=2)
    #     return max_y, min_y

    chr_cum_len_by_chrom = dict(zip(chr_names, chr_cum_lens))
    nrm_xs = []
    nrm_ys = []
    amp_xs = []
    amp_ys = []
    amp_gs = []
    del_xs = []
    del_ys = []
    del_gs = []
    with open(seq2c_tsv_fpath) as f:
        for i, l in enumerate(f):
            if i == 0: continue
            fs = l.replace('\n', '').split('\t')
            sname, gname = fs[0], fs[1]
            if key_gene_names and gname not in key_gene_names: continue
            if sname != sample_name: continue

            sname, gname, chrom, start, end, length, log2r, sig, type_, amp_del, ab_seg, total_seg, \
                ab_log2r, log2r_diff, ab_seg_loc, ab_samples, ab_samples_pcnt = fs[:17]
            x = chr_cum_len_by_chrom[chrom] + (int(start)+int(end))/2

            if not ab_log2r or type_ == 'BP':  # breakpoint, meaning part of exon is not amplified
                nrm_xs.append(x)
                nrm_ys.append(float(log2r))
                # add_rec_to_plot(chrom, start, end, log2r, max_y, min_y, marker='b.')

            if ab_log2r:
                y = float(ab_log2r)
                if amp_del == 'Amp':
                    amp_xs.append(x)
                    amp_ys.append(y)
                    amp_gs.append(gname)
                elif amp_del == 'Del':
                    del_xs.append(x)
                    del_ys.append(y)
                    del_gs.append(gname)
                else:
                    warn('Event is not Amp or Del, it\'s ' + amp_del)

                # max_y, min_y = add_rec_to_plot(chrom, start, end, log2r, max_y, min_y, marker=color + 'o', label=gname)

                # log2r = float(log2r)
                # if -0.5 < log2r < 0.5:
                #     color = 'k'
                # elif -1.5 < log2r < 1.5:
                #     color = 'g'
                # else:
                #     color = 'r'

    matplotlib.pyplot.scatter(nrm_xs, nrm_ys, marker='.', color='k', s=1)
    matplotlib.pyplot.scatter(amp_xs, amp_ys, marker='o', color='b', s=2)
    matplotlib.pyplot.scatter(del_xs, del_ys, marker='o', color='r', s=2)
    if len(amp_xs) <= 10 or len(amp_xs) + len(del_xs) < 40:
        for x, y, g in zip(amp_xs, amp_ys, amp_gs):
            ax.text(x, y, g, fontsize=9, color='g',
                    verticalalignment='center',
                    horizontalalignment='center')
    if len(del_xs) <= 10 or len(amp_xs) + len(del_xs) < 40:
        for x, y, g in zip(del_xs, del_ys, del_gs):
            ax.text(x, y, g, fontsize=9, color='r',
                    verticalalignment='center',
                    horizontalalignment='center')

    matplotlib.pyplot.ylim(
        ymax=max(chain(nrm_ys, amp_ys, del_ys, [2])) * 1.05,
        ymin=min(chain(nrm_ys, amp_ys, del_ys, [-2])) * 1.05)
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
