import shutil
import textwrap
from os import mkdir, makedirs
from os.path import basename, join, isdir, dirname, expanduser

from source.utils import info, verify_file, human_sorted
from source.logger import step_greetings

import matplotlib
import matplotlib.pyplot


distr_plot_ending = '.variant_distribution.png'


def variants_distribution_plot(cnf, vcf_fpath):
    step_greetings('Quality control variant distribution plots')

    # step 1: get chr lengths
    chr_len_fpath = cnf.get('chr_lengths')
    if not chr_len_fpath:
        chr_len_fpath = expanduser(cnf['genome'].get('chr_lengths'))
    if not verify_file(chr_len_fpath):
        exit(1)
        # no chromosome lengths file for the genome!
        # TODO: process reference fasta and get lengths from it
    chr_lengths = _get_chr_lengths(chr_len_fpath)

    # step 2: get variants distribution (per chromosome)
    qc_cnf = cnf['quality_control']
    variants_per_kbp = qc_cnf.get('variant_distribution_scale')
    plot_scale = 1000 * variants_per_kbp
    variants_distribution, not_counted = _get_variants_distribution(vcf_fpath, chr_lengths, plot_scale)
    if not_counted:
        info('Warning: some variants were not counted (chromosome names not found): ' + str(not_counted))
    empty_chr = []
    for chr_name in chr_lengths.keys():
        if sum(variants_distribution[chr_name]) == 0:
            empty_chr.append(chr_name)
            del variants_distribution[chr_name]
    if empty_chr:
        info('Chromosomes without variants: ' + ', '.join(human_sorted(empty_chr)))

    # step 3: plotting
    nplots = len(variants_distribution.keys())
    ncols = min(4, nplots)
    nrows = 1 + (nplots - 1) / 4
    fontsize = 15
    mbp = 1000000
    fig = matplotlib.pyplot.figure(figsize=(6 * ncols, 3 * nrows))
    for id, chr_name in enumerate(human_sorted(variants_distribution.keys())):
        #fig = matplotlib.pyplot.figure(figsize=(25,6))

        ax = fig.add_subplot(nrows, ncols, id + 1)
        ax.plot(range(len(variants_distribution[chr_name])), variants_distribution[chr_name], color='#46a246', linewidth=3.0)
        ax.axhline(linewidth=3.0)
        ax.axvline(linewidth=3.0)

        #Axe style###############################
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
        fig.suptitle('Number of variants from chromosomes not found in reference: ' + str(not_counted),
                     x=0.12, y=0.07, fontsize=fontsize, fontweight='bold', ha='left')
        fig.suptitle('Chromosomes without variants:', x=0.12, y=0.05, fontsize=fontsize, fontweight='bold', ha='left')
        fig.suptitle(chr_without_variants_wrapped, x=0.12, y=0.03, fontsize=fontsize, ha='left')

    variants_distribution_plot_fpath = join(cnf.output_dir, cnf.name + distr_plot_ending)
    fig.savefig(variants_distribution_plot_fpath, bbox_inches='tight')
    matplotlib.pyplot.close(fig)
    return variants_distribution_plot_fpath


def _get_chr_lengths(chr_len_fpath):
    chr_lengths = dict()
    with open(chr_len_fpath, 'r') as f:
        for line in f:
            if len(line.split()) == 2:
                chr_name = line.split()[0]
                chr_length = int(line.split()[1])
                chr_lengths[chr_name] = chr_length
    return chr_lengths


def _get_variants_distribution(vcf_fpath, chr_lengths, plot_scale):
    variants_distribution = dict()
    for chr_name, chr_length in chr_lengths.items():
        variants_distribution[chr_name] = [0] * max(1, chr_length / plot_scale)
    not_counted = 0

    with open(vcf_fpath, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.split()) < 8:
                continue
            chr_name = line.split()[0]
            chr_pos = int(line.split()[1])
            if chr_name not in variants_distribution:
                not_counted += 1
            else:
                region_id = min((chr_pos - 1) / plot_scale, len(variants_distribution[chr_name]) - 1)
                variants_distribution[chr_name][region_id] += 1

    # the last region in each chromosome is not exactly equal to plot_scale
    for chr_name, chr_length in chr_lengths.items():
        last_region_length = chr_length % plot_scale + (0 if chr_length < plot_scale else plot_scale)
        variants_distribution[chr_name][-1] = int(variants_distribution[chr_name][-1] * plot_scale /
                                                  float(last_region_length))
    return variants_distribution, not_counted

