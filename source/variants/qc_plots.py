import textwrap
from os.path import join
from source.file_utils import open_gzipsafe

from source.utils import human_sorted, get_chr_lengths
from source.logger import step_greetings, info
from ext_modules import vcf_parser as vcf

import matplotlib
import matplotlib.pyplot
from collections import OrderedDict
import numpy

distr_plot_ending = '.variant_distribution.png'
indels_plot_ending = '.indels.png'
substs_plot_ending = '.substitutions.png'


def draw_plots(cnf, vcf_fpath):
    step_greetings('Quality control plots')

    chr_lengths = get_chr_lengths(cnf)
    qc_cnf = cnf['quality_control']
    variants_per_kbp = qc_cnf.get('variant_distribution_scale')
    plot_scale = 1000 * variants_per_kbp

    info()
    info('Subsitutions and indel stats...')
    variants_distribution, substituitions, indel_lengths = _get_subs_and_indel_stats(vcf_fpath, chr_lengths, plot_scale)
    substs_plot_fpath = _draw_substitutions(cnf, substituitions)
    indels_plot_fpath = _draw_indel_lengths(cnf, indel_lengths)
    if substs_plot_fpath:
        info('  Substitutions: ' + substs_plot_fpath)
    if indels_plot_fpath:
        info('  Indels:        ' + indels_plot_fpath)
    variants_distribution_plot_fpath = _draw_variants_distribution(cnf, variants_distribution, chr_lengths, variants_per_kbp)
    if variants_distribution_plot_fpath:
        info('  Variant distr: ' + variants_distribution_plot_fpath)
    return [variants_distribution_plot_fpath, substs_plot_fpath, indels_plot_fpath]


def _get_subs_and_indel_stats(vcf_fpath, chr_lengths, plot_scale):
    reader = vcf.Reader(open_gzipsafe(vcf_fpath, 'r'))

    variants_distribution = dict()
    for chr_name, chr_length in chr_lengths.items():
        variants_distribution[chr_name] = [0] * max(1, chr_length / plot_scale)
    variants_distribution['OTHER'] = 0

    substituitions = OrderedDict()
    nucleotides = ['A', 'C', 'G', 'T']
    for nucl1 in nucleotides:
        substituitions[nucl1] = OrderedDict()
        for nucl2 in nucleotides:
            if nucl1 != nucl2:
                substituitions[nucl1][nucl2] = 0
    indel_lengths = []

    for rec in reader:
        # for variants distribution plot
        if rec.CHROM not in variants_distribution:
            variants_distribution['OTHER'] += 1
        else:
            region_id = min((rec.POS - 1) / plot_scale, len(variants_distribution[chr_name]) - 1)
            variants_distribution[rec.CHROM][region_id] += 1
        # for substitution and indel plots
        for alt in rec.ALT:
            if rec.is_snp:
                substituitions[rec.REF][str(alt)] += 1
            elif rec.is_indel:
                if alt is None:
                    indel_lengths.append(-1)
                else:
                    indel_lengths.append(len(alt) - len(rec.REF))

    # the last region in each chromosome is not exactly equal to plot_scale
    for chr_name, chr_length in chr_lengths.items():
        last_region_length = chr_length % plot_scale + (0 if chr_length < plot_scale else plot_scale)
        variants_distribution[chr_name][-1] = int(variants_distribution[chr_name][-1] * plot_scale /
                                                  float(last_region_length))
    return variants_distribution, substituitions, indel_lengths


def _draw_variants_distribution(cnf, variants_distribution, chr_lengths, variants_per_kbp):
    plot_fpath = join(cnf.output_dir, cnf.name + '-' + cnf.caller + distr_plot_ending)

    not_counted = variants_distribution['OTHER']
    del variants_distribution['OTHER']
    if not_counted:
        info('Warning: some variants were not counted (chromosome names not found): ' + str(not_counted))
    empty_chr = []
    for chr_name in chr_lengths.keys():
        if sum(variants_distribution[chr_name]) == 0:
            empty_chr.append(chr_name)
            del variants_distribution[chr_name]
    if empty_chr:
        info('  Chromosomes without variants: ' + ', '.join(human_sorted(empty_chr)))

    nplots = len(variants_distribution.keys())
    ncols = min(4, nplots)
    nrows = 1 + (nplots - 1) / 4
    fontsize = 15
    mbp = 10**6
    fig = matplotlib.pyplot.figure(figsize=(6 * ncols, 3 * nrows))
    for id, chr_name in enumerate(human_sorted(variants_distribution.keys())):
        ax = fig.add_subplot(nrows, ncols, id + 1)
        ax.plot(range(len(variants_distribution[chr_name])), variants_distribution[chr_name], color='#46a246', linewidth=3.0)
        ax.axhline(linewidth=3.0)
        ax.axvline(linewidth=3.0)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xticklabels('')
        ax.set_ylim(bottom=0)
        if chr_lengths[chr_name] < mbp / 10:
            chr_size = '<0.1 Mbp'
        else:
            chr_size = '%.1f Mbp' % (float(chr_lengths[chr_name]) / mbp)
        ax.set_xlabel(chr_name + ', ' + chr_size)

        for item in ([ax.xaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            item.set_weight('bold')

        #TODO: individual plots for each chromosome separately
        #fig.savefig(join(qc_dir, cnf['name'] + '_qc_' + chr_name + '.png'))
        #matplotlib.pyplot.close(fig)

    fig.suptitle('Variants per %d kbp' % variants_per_kbp, y=0.93, fontsize=int(1.5 * fontsize), fontweight='bold')

    if empty_chr:
        if len(empty_chr) > 1:
            chr_without_variants = ', '.join(human_sorted(empty_chr))
            renderer = fig.canvas.get_renderer()
            aspect_ratio = 0.5 # This varies with the font!!
            pixels_per_char = aspect_ratio * renderer.points_to_pixels(fontsize)
            new_width = int(fig.dpi * fig.get_figwidth() * 0.7)
            wrap_width = max(1, new_width // pixels_per_char)
            chr_without_variants_wrapped = textwrap.fill(chr_without_variants, wrap_width)
        else:
            chr_without_variants_wrapped = empty_chr[0]
        fig.suptitle('Variants on other chromosomes (not found in reference): ' + str(not_counted),
                     x=0.12, y=0.07, fontsize=fontsize, fontweight='bold', ha='left')
        fig.suptitle('Chromosomes without variants:', x=0.12, y=0.05, fontsize=fontsize, fontweight='bold', ha='left')
        fig.suptitle(chr_without_variants_wrapped, x=0.12, y=0.03, fontsize=fontsize, ha='left')

    fig.savefig(plot_fpath, bbox_inches='tight')
    matplotlib.pyplot.close(fig)
    return plot_fpath


def _draw_substitutions(cnf, substituitions):
    plot_fpath = join(cnf.output_dir, cnf.name + ('-' + cnf.caller if cnf.caller else '') + substs_plot_ending)

    colors = ['#CC0000', '#CC6600', '#CCCC00', '#66CC00']
    # params of bars
    width = 0.3
    interval = width / 3
    start_pos = interval / 2

    counts = []
    labels = []
    for nucl1 in substituitions.iterkeys():
        for nucl2, count in substituitions[nucl1].iteritems():
            counts.append(count)
            labels.append('%s>%s' % (nucl1, nucl2))
    total = sum(counts)

    def __to_percents(x):
        return float(x) * 100 / total

    positions = []
    ax = matplotlib.pyplot.gca()
    alty = ax.twinx()
    for i, val in enumerate(counts):
        positions.append(start_pos + (width + interval) * i)
        ax.bar(positions[-1], val, width, color=colors[i * len(colors) / len(counts)])
        alty.bar(positions[-1], __to_percents(val), width, color=colors[i * len(colors) / len(counts)])

    matplotlib.pyplot.title('Substitutions')
    ax.set_ylabel('Count')
    alty.set_ylabel('Rate, %')

    ax.yaxis.grid(True)
    ax.set_ylim(0, ax.get_ylim()[1])
    alty.set_ylim(0, __to_percents(ax.get_ylim()[1]))
    ax.set_xticks([position + float(width) / 2 for position in positions])
    ax.set_xticklabels(labels)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    #alty.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f %%'))

    matplotlib.pyplot.savefig(plot_fpath)
    matplotlib.pyplot.close()
    return plot_fpath


def _draw_indel_lengths(cnf, indel_lengths):
    plot_fpath = join(cnf.output_dir, cnf.name + '-' + cnf.caller + indels_plot_ending)

    bins = list(numpy.arange(min(indel_lengths) - 0.5, max(indel_lengths) + 1.0, 1.0))
    matplotlib.pyplot.hist(indel_lengths, bins=bins, color='#CC0000')
    matplotlib.pyplot.title('Indel length distribution')
    matplotlib.pyplot.ylim(0, matplotlib.pyplot.ylim()[1])
    ax = matplotlib.pyplot.gca()
    ax.set_ylabel('Count')
    ax.set_xlabel('Indel length, bp')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    matplotlib.pyplot.savefig(plot_fpath)
    matplotlib.pyplot.close()
    return plot_fpath

