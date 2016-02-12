import shutil
import os
from os import listdir
from os.path import relpath, join, exists, dirname, basename, abspath
from collections import OrderedDict, defaultdict

import source
from source.clinical_reporting.clinical_reporting import gray
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import verify_bed
from source.reporting.reporting import FullReport, Metric, MetricStorage, ReportSection, write_tsv_rows, PerRegionSampleReport, BaseReport, \
    write_static_html_report, build_report_html
from source.logger import step_greetings, info, warn, err
from source.targetcov.bam_and_bed_utils import prepare_beds
from source.targetcov.cov import get_detailed_metric_storage, get_header_metric_storage



# from source.targetcov.flag_regions import DepthsMetric
from source.targetcov.flag_regions import DepthsMetric, _intersect_with_tricky_regions
from source.tools_from_cnf import get_system_path, get_qualimap_type
from source.calling_process import call
from source.file_utils import safe_mkdir, verify_file, verify_dir, intermediate_fname, symlink_plus, adjust_path, \
    file_transaction, add_suffix
from source.bcbio.bcbio_structure import BCBioStructure
from source.variants.vcf_processing import bgzip_and_tabix


def _run_multisample_qualimap(cnf, output_dir, samples, targqc_full_report):
    """ 1. Generates Qualimap2 plots and put into plots_dirpath
        2. Adds records to targqc_full_report.plots
    """
    plots_dirpath = join(output_dir, 'plots')
    if cnf.reuse_intermediate and verify_dir(plots_dirpath) and [f for f in listdir(plots_dirpath) if not f.startswith('.')]:
        info('Qualimap miltisample plots exist - ' + plots_dirpath + ', reusing...')
    else:
        # Qualimap2 run for multi-sample plots
        if len([s.qualimap_html_fpath for s in samples if s.qualimap_html_fpath]) > 0:
            qualimap = get_system_path(cnf, interpreter_or_name=None, name='qualimap')

            if qualimap is not None and get_qualimap_type(qualimap) == 'full':
                qualimap_output_dir = join(cnf.work_dir, 'qualimap_multi_bamqc')

                _correct_qualimap_genome_results(cnf, samples)
                _correct_qualimap_insert_size_histogram(cnf, samples)

                safe_mkdir(qualimap_output_dir)
                rows = []
                for sample in samples:
                    if sample.qualimap_html_fpath:
                        rows += [[sample.name, sample.qualimap_html_fpath]]

                data_fpath = write_tsv_rows(rows, join(qualimap_output_dir, 'qualimap_results_by_sample.tsv'))
                qualimap_plots_dirpath = join(qualimap_output_dir, 'images_multisampleBamQcReport')
                cmdline = '{qualimap} multi-bamqc --data {data_fpath} -outdir {qualimap_output_dir}'.format(**locals())
                res = call(cnf, cmdline, exit_on_error=False, return_err_code=True, env_vars=dict(DISPLAY=None),
                           output_fpath=qualimap_plots_dirpath, output_is_dir=True)
                if res is None or not verify_dir(qualimap_plots_dirpath):
                    warn('Warning: Qualimap for multi-sample analysis failed to finish. TargQC will not contain plots.')
                    return None
                else:
                    if exists(plots_dirpath):
                        shutil.rmtree(plots_dirpath)
                    shutil.move(qualimap_plots_dirpath, plots_dirpath)
            else:
                warn('Warning: Qualimap for multi-sample analysis was not found. TargQC will not contain plots.')
                return None

    targqc_full_report.plots = []
    for plot_fpath in listdir(plots_dirpath):
        plot_fpath = join(plots_dirpath, plot_fpath)
        if verify_file(plot_fpath) and plot_fpath.endswith('.png'):
            targqc_full_report.plots.append(relpath(plot_fpath, output_dir))


def _make_targetcov_symlinks(samples):
    for sample in samples:
        new_link = join(
            dirname(dirname(sample.targetcov_detailed_txt)),
            basename(sample.targetcov_detailed_txt))
        if exists(new_link):
            os.unlink(new_link)
        symlink_plus(sample.targetcov_detailed_txt, new_link)
        info('TargetCov TXT symlink saved to ' + new_link)


def _make_tarqc_html_report(cnf, output_dir, samples, bed_fpath=None, tag_by_sample=None):
    header_storage = get_header_metric_storage(cnf.coverage_reports.depth_thresholds,
                                               is_wgs=bed_fpath is not None,
                                               padding=cnf.coverage_reports.padding)

    # targqc_metric_storage = _get_targqc_metric_storage([
    #     ('targetcov', header_storage),
    #     ('ngscat', ngscat_report_parser.metric_storage)])
    #     # ('qualimap', qualimap_report_parser.metric_storage)])

    jsons_by_sample = {s.name: s.targetcov_json_fpath for s in samples if verify_file(s.targetcov_json_fpath)}
    htmls_by_sample = {s.name: s.targetcov_html_fpath for s in samples if verify_file(s.targetcov_html_fpath)}

    targqc_full_report = FullReport.construct_from_sample_report_jsons(samples, output_dir, jsons_by_sample, htmls_by_sample)

    # source.targqc_repr, [], metric_storage=targqc_metric_storage)

    for sample_report in targqc_full_report.sample_reports:
    #     records_by_report_type = []
    #     if (verify_file(sample.targetcov_json_fpath, True) or
    #             verify_file(sample.ngscat_html_fpath, True) or
    #             verify_file(sample.qualimap_html_fpath, True)):
    #         records_by_report_type.append(('targetcov', load_records(sample.targetcov_json_fpath) if verify_file(
    #             sample.targetcov_json_fpath, silent=True) else []))
    #         records_by_report_type.append(('ngscat', ngscat_report_parser.parse_ngscat_sample_report(
    #             sample.ngscat_html_fpath) if verify_file(sample.ngscat_html_fpath, silent=True) else []))
    #         records_by_report_type.append(('qualimap', qualimap_report_parser.parse_qualimap_sample_report(
    #             sample.qualimap_html_fpath) if verify_file(sample.qualimap_html_fpath, silent=True) else []))
    #
    #     sample_report = SampleReport(
    #         sample,
    #         records=_get_targqc_records(records_by_report_type, header_storage),
    #         html_fpath=dict(
    #             targetcov=relpath(sample.targetcov_html_fpath, output_dir) if sample.targetcov_html_fpath else None,
    #             ngscat=relpath(sample.ngscat_html_fpath, output_dir) if sample.ngscat_html_fpath else None,
    #             qualimap=relpath(sample.qualimap_html_fpath, output_dir) if sample.qualimap_html_fpath else None
    #         ),
    #         metric_storage=targqc_metric_storage)
        if tag_by_sample:
            sample_report.set_project_tag(tag_by_sample[sample_report.sample.name])
        if verify_file(sample_report.sample.ngscat_html_fpath):
            sample_report.add_record(metric_name='ngsCAT', value='ngsCAT', html_fpath=sample_report.sample.ngscat_html_fpath, silent=True)
        if verify_file(sample_report.sample.qualimap_html_fpath):
            url = relpath(sample_report.sample.qualimap_html_fpath, output_dir)
            r = sample_report.find_record(sample_report.records, 'Qualimap')
            if r:
                r.url = url
            else:
                sample_report.add_record(metric_name='Qualimap', value='Qualimap', url=url, silent=True)
        # targqc_full_report.sample_reports.append(sample_report)

    _run_multisample_qualimap(cnf, output_dir, samples, targqc_full_report)

    # orig_bed_rec = next((r for r in targqc_full_report.get_common_records() if r.metric.name == 'Target'), None)
    # ready_bed_rec = next((r for r in targqc_full_report.get_common_records() if r.metric.name == 'Target ready'), None)

    # if not ready_bed_rec:
    #     ready_bed_rec = orig_bed_rec

    # if ready_bed_rec:
    #     ready_bed = ready_bed_rec.path
    #     if verify_bed(ready_bed, 'ready_bed_rec.value'):
    #         project_ready_bed = join(output_dir, 'target.bed')
    #         shutil.copy(ready_bed, project_ready_bed)
    #         ready_bed_rec.value = project_ready_bed

    # if orig_bed_rec and ready_bed_rec:
    #     orig_bed_rec.value = bed_fpath

    txt_fpath = targqc_full_report.save_txt(join(output_dir, BCBioStructure.targqc_name + '.txt'))
    tsv_fpath = targqc_full_report.save_tsv(join(output_dir, BCBioStructure.targqc_name + '.tsv'))
    html_fpath = targqc_full_report.save_html(cnf, join(output_dir, BCBioStructure.targqc_name + '.html'), 'TargQC')

    return txt_fpath, tsv_fpath, html_fpath


def summarize_targqc(cnf, summary_threads, output_dir, samples, bed_fpath=None, exons_fpath=None, tag_by_sample=None):
    step_greetings('Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    correct_samples = []

    for sample in samples:
        # if not sample.targetcov_done():
        #     err('Error: target coverage is not done (json, html, or detail tsv are not there)')
        # else:
        correct_samples.append(sample)
        # if not sample.ngscat_done():
        # sample.ngscat_html_fpath = None
        # if not sample.qualimap_done():
        # sample.qualimap_html_fpath = None
    samples = correct_samples

    # _make_targetcov_symlinks(samples)

    txt_fpath, tsv_fpath, html_fpath = _make_tarqc_html_report(
            cnf, output_dir, samples, bed_fpath, tag_by_sample=tag_by_sample)

    best_for_regions_fpath = None
    if any(verify_file(s.targetcov_detailed_tsv, silent=True) for s in samples):
        best_for_regions_fpath = _save_best_details_for_each_gene(cnf.coverage_reports.depth_thresholds, samples, output_dir)
    ''' 1. best_regions = get_best_regions()
        2. best_for_regions_fpath = save_per_region_report()
        3. calc median coverage across best regions
        4. flagged_regions_report_fpath = _generate_flagged_regions_report(
             output_dir, 'Best', average_coverage, genes, depth_threshs)
    '''

    if cnf.extended:
        if not exons_fpath or not bed_fpath:
            err('For the extended analysis, capture and exons beds are required!')
        else:
            exons_bed, exons_no_genes_cut_bed, target_bed, _ = prepare_beds(cnf, exons_fpath, bed_fpath)

            #norm_best_var_fpath, norm_comb_var_fpath = _report_normalize_coverage_for_variant_sites(
            #    cnf, summary_threads, output_dir, samples, 'oncomine', bed_fpath)

    info()
    info('*' * 70)
    if not html_fpath and not txt_fpath:
        info('TargQC summary was not generated, because there were no reports generated for individual samples.')
    else:
        info('TargQC summary saved in: ')
        for fpath in [txt_fpath, html_fpath]:
            if fpath: info('  ' + fpath)

    if best_for_regions_fpath:
        info()
        info('Best stats for regions saved in:')
        info('  ' + best_for_regions_fpath)

    # if cnf.extended:
    #     if norm_best_var_fpath:
    #         info()
    #         info('Normalized depths for oncomine saved in:')
    #         info('        ' + norm_comb_var_fpath)
    #         info('  Best: ' + norm_best_var_fpath)

    return html_fpath


def _generate_summary_flagged_regions_report(cnf, bcbio_structure, samples, mutations, key_or_target_genes):
    region_types = ['exons', 'amplicons']
    coverage_types = ['low', 'high']
    flagged_regions_metrics = [
            Metric('Gene', min_width=50, max_width=70),
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Position', td_class='td_position', min_width=70, max_width=120),
            Metric('Ave depth', td_class='long_expanded_line right_aligned', max_width=100, with_heatmap=False),
            Metric('#HS', quality='Less is better', align='right', max_width=30),
            Metric('Hotspots & Deleterious', td_class='long_expanded_line',  min_width=100, max_width=150),
            Metric('Found mutations', td_class='long_expanded_line', min_width=150, max_width=200),
            Metric('Samples', td_class='long_expanded_line', min_width=100, max_width=120),
            Metric('Possible reasons', td_class='long_expanded_line', max_width=120)]
    flagged_regions_metric_storage = MetricStorage(sections=[ReportSection(metrics=flagged_regions_metrics)])
    flagged_regions_report_dirpath = bcbio_structure.flagged_regions_dirpath
    safe_mkdir(flagged_regions_report_dirpath)
    if key_or_target_genes == 'target':
        genes_description = 'genes'
    else:
        genes_description = 'genes that have been previously implicated in various cancers'
    for region_type in region_types:
        regions_dict = {}
        total_regions = 0
        info()
        info('Preparing report for ' + region_type)
        for coverage_type in coverage_types:
            regions_by_gene = {}
            for sample in samples:
                selected_regions_bed_fpath = join(sample.flagged_regions_dirpath, coverage_type + '_cov_' + region_type + '.bed')
                regions_by_reasons = {}
                if verify_file(selected_regions_bed_fpath, is_critical=False):
                    regions_by_reasons = _intersect_with_tricky_regions(cnf, selected_regions_bed_fpath, sample.name)
                total_report_fpath = add_suffix(add_suffix(sample.flagged_tsv, region_type), coverage_type)
                if verify_file(total_report_fpath, is_critical=False):
                    with open(total_report_fpath) as f:
                        for l in f:
                            l = l.strip()
                            if not l or l.startswith('#'):
                                continue
                            fs = l.split('\t')
                            (chrom, start, end, size, gene, strand, feature, biotype, min_depth, avg_depth) = fs[:10]
                            start, end = int(start), int(end)
                            regions_by_gene.setdefault(gene, [])
                            cur_region = Region(sample_name=[sample.name], avg_depth=[avg_depth], gene_name=gene, strand=strand, feature=feature,
                                                biotype=biotype, chrom=chrom, start=start, end=end)
                            for r in regions_by_reasons:
                                if r[0] <= start and end <= r[1]:
                                    cur_region.extra_fields = regions_by_reasons[r]
                            cur_region.missed_by_db = []
                            was_added = False
                            for r in regions_by_gene[gene]:
                                if r.start <= cur_region.start <= r.end and r.start <= cur_region.end <= r.end:
                                    was_added = True
                                    if sample.name not in r.sample_name:
                                        r.sample_name.append(sample.name)
                                        r.avg_depth.append(avg_depth)
                            if not was_added:
                                regions_by_gene[gene].append(cur_region)
                report_fpath = join(sample.flagged_regions_dirpath, coverage_type + '_cov_' + region_type + '.oncomine.tsv')
                if verify_file(report_fpath, is_critical=False):
                    with open(report_fpath) as f:
                        for l in f:
                            l = l.strip()
                            if not l or l.startswith('#'):
                                continue
                            fs = l.split('\t')
                            hotspots = []
                            (gene, chrom, start, end, strand, feature, biotype, id_, num_hotspots) = fs[:9]
                            start, end = int(start), int(end)
                            if int(num_hotspots) != 0:
                                hotspots = fs[9].split()

                            regions_by_gene.setdefault(gene, [])
                            cur_region = Region(sample_name=[sample.name], gene_name=gene, strand=strand, feature=feature,
                                                biotype=biotype, chrom=chrom, start=start, end=end)
                            for r in regions_by_gene[gene]:
                                if r.start <= cur_region.start <= r.end and r.start <= cur_region.end <= r.end:
                                    if sample.name not in r.sample_name:
                                        r.sample_name.append(sample.name)
                                        r.avg_depth.append('.')
                                    new_hotspots = [hs for hs in hotspots if hs not in r.missed_by_db]
                                    r.missed_by_db.extend(new_hotspots)
            flagged_regions_report = PerRegionSampleReport(name='Flagged regions', metric_storage=flagged_regions_metric_storage)
            num_regions = 0
            non_hs_class = ' no_hotspots'
            slash_with_zero_space = '/&#x200b;'
            for gene in regions_by_gene.keys():
                if regions_by_gene[gene]:
                    num_regions += len(regions_by_gene[gene])
                    row_class = ' expandable_row collapsed'
                    if len(regions_by_gene[gene]) > 1:
                        reg = flagged_regions_report.add_row()
                        reg.class_ = ' expandable_gene_row collapsed'
                        chr = regions_by_gene[gene][0].chrom
                        num_hotspots = [len(r.missed_by_db) for r in regions_by_gene[gene]]
                        all_samples = [sample for r in regions_by_gene[gene] for sample in r.sample_name]
                        all_unique_samples = []
                        all_unique_samples = [sample for sample in all_samples if sample not in all_unique_samples and not all_unique_samples.append(sample)]
                        all_tricky_regions = sorted(set([tricky_region for r in regions_by_gene[gene] for tricky_region in r.extra_fields]))
                        all_depths = [[] for x in range(len(all_unique_samples))]
                        for r in regions_by_gene[gene]:
                            for sample_num, sample in enumerate(all_unique_samples):
                                if sample in r.sample_name:
                                    cur_sample_index = r.sample_name.index(sample)
                                    if r.avg_depth[cur_sample_index] != '.':
                                        all_depths[sample_num].append(float(r.avg_depth[cur_sample_index]))
                        avg_depth_per_samples = [sum(all_depths[i])/len(all_depths[i]) if len(all_depths[i]) > 0 else 0 for i in range(len(all_depths)) ]
                        reg.add_record('Gene', gene)
                        reg.add_record('Chr', chr.replace('chr', ''))
                        reg.add_record('#HS', sum(num_hotspots))
                        reg.add_record('Position', str(len(regions_by_gene[gene])) + ' regions')
                        reg.add_record('Ave depth', slash_with_zero_space.join([format(depth, '.2f') if depth != '.' else '.' for depth in avg_depth_per_samples]),
                                       num=sum(avg_depth_per_samples)/len(avg_depth_per_samples))
                        reg.add_record('Hotspots & Deleterious', '')
                        reg.add_record('Possible reasons', ', '.join(all_tricky_regions))
                        reg.add_record('Samples', ',\n'.join(all_unique_samples))
                        reg.add_record('Found mutations', '')
                        if sum(num_hotspots) == 0:
                            reg.class_ += non_hs_class
                        row_class += ' row_to_hide row_hidden'
                    else:
                        row_class += ' not_to_hide'
                    for r in regions_by_gene[gene]:
                        reg = flagged_regions_report.add_row()
                        reg.class_ = row_class
                        reg.add_record('Gene', r.gene_name)
                        reg.add_record('Chr', r.chrom.replace('chr', ''))
                        avg_depths = [float(depth) for depth in r.avg_depth if depth != '.']
                        reg.add_record('Ave depth', slash_with_zero_space.join([format(depth, '.2f') if depth != '.'
                                                                                else depth for depth in avg_depths]), num=sum(avg_depths)/len(avg_depths))
                        reg.add_record('Position', Metric.format_value(r.start, human_readable=True, is_html=True) + '-' +
                                       Metric.format_value(r.end, human_readable=True, is_html=True))
                        reg.add_record('#HS', len(r.missed_by_db))
                        if len(r.missed_by_db) == 0:
                            reg.class_ += non_hs_class
                        uniq_hs_positions = sorted(set([hotspot.split(':')[0] for hotspot in r.missed_by_db]))
                        hs_by_pos = {pos: [h.split(':')[1] for h in r.missed_by_db if h.split(':')[0] == pos] for pos in uniq_hs_positions}
                        hs_breakable = [gray(Metric.format_value(int(pos.replace(',','')), human_readable=True, is_html=True))
                                        + ': ' + ','.join([h.replace('/', slash_with_zero_space) for h in hs_by_pos[pos]]) for pos in uniq_hs_positions]
                        reg.add_record('Hotspots & Deleterious', '\n'.join(hs_breakable))
                        reg.add_record('Possible reasons', ', '.join(r.extra_fields))
                        reg.add_record('Samples', ',\n'.join(r.sample_name))
                        found_mutations = []
                        for sample in samples:
                            if sample.name in r.sample_name:
                                for mut in mutations[sample.name]:
                                    if mut.gene.name == r.gene_name and r.start <= mut.pos <= r.end:
                                        found_mutations.append(gray(Metric.format_value(mut.pos, human_readable=True, is_html=True)) + ':' + mut.ref + '>' + mut.alt + ' (' + sample.name +')')
                        reg.add_record('Found mutations', '\n'.join(found_mutations))
            flagged_regions_report.expandable = True
            flagged_regions_report.unique = True
            regions_dict[coverage_type] = create_section(flagged_regions_report, num_regions, regions_by_gene.keys(), region_type)
            total_regions += num_regions
        flagged_report_fpath = join(flagged_regions_report_dirpath, 'flagged_' + region_type + '.html')
        write_static_html_report(cnf, {
            'key_or_target': key_or_target_genes,
            'region_type': region_type,
            'genes_description': genes_description,
            'flagged_low': regions_dict['low'],
            'flagged_high': regions_dict['high'],
        }, flagged_report_fpath, tmpl_fpath=join(dirname(abspath(__file__)), 'template_flagged_regions.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'flagged_regions.js')],
           extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'flagged_regions.css')])
        #BaseReport.save_html(flagged_regions_report, cnf, flagged_report_fpath, caption='Flagged regions')
        info('')
        info('Flagged regions (total ' + str(total_regions) + ' ' + region_type + ') saved into:')
        info('  ' + flagged_report_fpath)


def create_section(report, num_regions, genes, region_type):
    flagged_dict = dict()
    if report.rows:
        # if cnf.debug:
        #     mutations_report.regions = mutations_report.regions[::20]
        flagged_dict['table'] = build_report_html(report, sortable=True)
        flagged_dict['total_regions'] = Metric.format_value(num_regions, is_html=True)
        flagged_dict['total_key_genes'] = Metric.format_value(len(genes), is_html=True)
        flagged_dict['region_type'] = Metric.format_value(region_type, is_html=True)
    return flagged_dict

# def _generate_flagged_regions_report(output_dir, sample, genes, depth_threshs):
#     report = PerRegionSampleReport(sample=sample, metric_storage=get_detailed_metric_storage(depth_threshs))
#     report.add_record('Sample', sample.name)
#
#     ''' 1. Detect depth threshold (ave sample coverage/2 but > 25x)
#         2. Select regions covered in less than 80% at threshold
#         3. Sort by % at threshold
#         4. Select those parts of those regions where % = 0, save to BED
#         5. Find HotSpots at those regions
#         6. Intersect HotSpots with tracks
#     '''
#
#     ave_coverages_per_sample = {
#         s.name: get_ave_coverage(cnf, s.targetcov_json_fpath)
#         for s in samples if verify_file(s.targetcov_json_fpath)}
#
#     regions = []
#     for gene in genes:
#         regions.extend(gene.get_exons())
#
#     depth_cutoff = max(ave_depth / 4, 25)
#     for thresh in depth_threshs[::-1]:
#         if thresh < depth_cutoff:
#             depth_cutoff = thresh
#             break
#
#     sorted_by_thresh = sorted(regions, key=lambda r: [r.rates_within_threshs[t] for t in depth_threshs])
#
#     low_cov_regions = [r for r in sorted_by_thresh if r.rates_within_threshs[depth_cutoff] < 1]
#
#     selected_regions = low_cov_regions
#
#     report = make_flat_region_report(sample, selected_regions, depth_threshs)
#
#     gene_report_basename = sample.name + '.' + source.targetseq_name + '.selected_regions'
#     txt_rep_fpath = report.save_txt(join(output_dir, gene_report_basename + '.txt'))
#     tsv_rep_fpath = report.save_tsv(join(output_dir, gene_report_basename + '.tsv'))
#     info('')
#     info('Selected regions (total ' + str(len(selected_regions)) + ') saved into:')
#     info('  ' + txt_rep_fpath)
#
#     return txt_rep_fpath


def _clip_vcf_by_bed(cnf, vcf_fpath, bed_fpath):
    info('Clipping VCF ' + vcf_fpath + ' using BED ' + bed_fpath)

    bedtools = get_system_path(cnf, 'bedtools')

    clipped_vcf_fpath = intermediate_fname(cnf, vcf_fpath, 'clip')
    cmdline = '{bedtools} intersect -header	-a {vcf_fpath} -b {bed_fpath}'.format(**locals())
    res = call(cnf, cmdline, output_fpath=clipped_vcf_fpath)

    clipped_gz_vcf_fpath = bgzip_and_tabix(cnf, clipped_vcf_fpath)

    return clipped_gz_vcf_fpath


def _get_depth_for_each_variant(cnf, samtools, bedtools,
        var_by_site, clipped_gz_vcf_fpath, bed_fpath,
        sample_name, bam_fpath):
    info()
    info('Processing sample ' + sample_name)

    # http://www.1000genomes.org/faq/what-depth-coverage-your-phase1-variants
    # bedtools intersect -a oncomine.vcf -b Exons.az_key.bed -header > oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/bgzip oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/tabix -h -p vcf oncomine.az_key.vcf.gz
    # samtools view -b TRF004223.sorted.bam -L Exons.az_key.bed | bedtools genomecov -ibam stdin -bg > coverage.bg
    # bedtools intersect -a oncomine.az_key.vcf.gz -b coverage.bg -wa | cut -f1,2,4,5,8,11,12,13,14 > oncomine.az_key.depth_numbers.vcf

    info()
    info('Depth of coverage for regions in BED ' + bed_fpath)
    cov_bg = join(cnf.work_dir, sample_name + '_coverage.bg')
    cmdline = '{samtools} view -L {bed_fpath} -b {bam_fpath} | {bedtools} genomecov -ibam stdin -bg'.format(**locals())
    call(cnf, cmdline, output_fpath=cov_bg)

    info()
    info('Intersection depth regions with VCF ' + clipped_gz_vcf_fpath)
    vcf_depth_numbers_fpath = join(cnf.work_dir, sample_name + '_vcf_bg.intersect')
    cmdline = '{bedtools} intersect -a {clipped_gz_vcf_fpath} -b {cov_bg} -wao'.format(**locals())
    res = call(cnf, cmdline, output_fpath=vcf_depth_numbers_fpath)
    # if res != oncomine_depth_numbers_fpath:
    #     info()
    #     info('Trying with uncompressed VCF')
    #     cmdline = 'gunzip {vcf_fpath} -c | {bedtools} intersect -a - -b {cov_bg} -wao | cut -f1,2,4,5,8,11,12,13,14,15'.format(**locals())
    #     call(cnf, cmdline, output_fpath=oncomine_depth_numbers_fpath)

    depths_per_var = defaultdict(list)
    with open(vcf_depth_numbers_fpath) as f:
        for l in f:
            # 1,2,4,5,8,11,12,13,14,15,16,17,18,19,20
            # c,p,r,a,f,ch,st,en,ge,ex,st,ft,bt,de,ov
            fs = l[:-1].split('\t')
            chrom, pos, _, ref, alt = fs[:5]
            depth, overlap = fs[-2:]
            var = var_by_site.get((chrom, pos, ref, alt))
            if var and depth != '.':
                depth, overlap = int(depth), int(overlap)
                for i in range(overlap):
                    depths_per_var[(chrom, pos, ref, alt)].append(depth)

    # Getting avarage depth of coverage of each variant (exactly for those parts that were in BED)
    depth_by_var = {var: (sum(depths) / len(depths)) if len(depths) != 0 else None
                    for var, depths in depths_per_var.iteritems()}

    return sample_name, depth_by_var


def _prep_comb_report(metric_storage, samples, shared_general_metrics, shared_metrics):
    comb_general_metrics = shared_general_metrics[:]
    comb_general_metrics.append(Metric('For each sample'))
    for s in samples:
        comb_general_metrics.append(Metric(s.name + ' ave depth'))

    comb_metrics = shared_metrics[:]
    for s in samples:
        comb_metrics.append(DepthsMetric(s.name + ' hotspots depths/norm depths', short_name=s.name))

    comb_report_metric_storage = MetricStorage(
        general_section=ReportSection('general_section', metrics=comb_general_metrics),
        sections=[ReportSection(metrics=comb_metrics)])

    report = PerRegionSampleReport(sample='Combined', metric_storage=comb_report_metric_storage)

    report.add_record('Sample', 'contains values from all samples: ' + ', '.join([s.name for s in samples]))
    report.add_record('For each sample', 'Depths and normalized depths for each hotspot.')

    m = metric_storage.find_metric('Average sample depth')
    for s in samples:
        val = BaseReport.find_record(s.report.records, m.name).value
        report.add_record(s.name + ' ave depth', val)

    return report


def _prep_best_report(metric_storage, samples):
    report = PerRegionSampleReport(sample='Best', metric_storage=metric_storage)

    report.add_record('Sample', 'contains best values from all samples: ' + ', '.join([s.name for s in samples]))

    m = metric_storage.find_metric('Average sample depth')
    ave_sample_depth = max(BaseReport.find_record(s.report.records, m.name).value for s in samples)
    report.add_record('Average sample depth', ave_sample_depth)

    return report


def _report_normalize_coverage_for_variant_sites(cnf, summary_threads, output_dir, samples, vcf_key, bed_fpath):
    step_greetings('Combined normalized coverage for ' + vcf_key + ' hotspots')

    # vcf_fpath = cnf.genome.get(vcf_key)
    # if not vcf_fpath:
    #     err('Error: no ' + vcf_key + ' for ' + cnf.genome.name + ' VCF fpath specified in ' + cnf.sys_cnf)
    #     return None
    #
    # ave_coverages_per_sample = {
    #     s.name: get_ave_coverage(cnf, s.targetcov_json_fpath)
    #     for s in samples if verify_file(s.targetcov_json_fpath)}
    #
    # clipped_gz_vcf_fpath, regions_in_order, vars_by_region, var_by_site = \
    #     _read_vars_per_region_and_clip_vcf(cnf, vcf_fpath, bed_fpath)
    #
    # samtools = get_system_path(cnf, 'samtools')
    # bedtools = get_system_path(cnf, 'bedtools')
    #
    # vars_by_region_per_sample = OrderedDict(Parallel(n_jobs=summary_threads)
    #    (delayed(_get_depth_for_each_variant)(
    #         CallCnf(cnf.__dict__), samtools, bedtools, var_by_site, clipped_gz_vcf_fpath, bed_fpath,
    #         s.name, s.bam)
    #     for s in samples))

    ############################ Combined ############################
    # info()
    # info('*' * 70)
    # info('Saving for all samples: combined reports.')
    # # best_report = _prep_best_report(single_report_metric_storage, samples)
    # comb_report = _prep_comb_report(single_report_metric_storage, samples, shared_general_metrics, shared_metrics)
    #
    # total_variants = 0
    # nth_regions_from_each_sample = [s.report.get_regions() for s in samples]
    # while True:
    #     nth_region_from_each_sample = [rs[total_variants] for rs in nth_regions_from_each_sample if total_variants < len(rs)]
    #     total_variants += 1
    #     if len(nth_region_from_each_sample) == 0:
    #         break
    #     assert len(nth_region_from_each_sample) == len(nth_regions_from_each_sample), 'Region files for samples are not euqal size'
    #
    #     # best_report_reg = best_report.add_region()
    #     comb_report_reg = comb_report.add_region()
    #     rand_line = nth_region_from_each_sample[0]
    #     for i in range(9):
    #         # best_report_reg.records.append(rand_line.records[i])
    #         comb_report_reg.records.append(rand_line.records[i])
    #
    #     # best_depth = select_best(r.records[10].value for r in nth_region_from_each_sample)
    #     # best_norm_depth = select_best(r.records[11].value for r in nth_region_from_each_sample)
    #     # best_report_reg.add_record('Depth', best_depth)
    #     # best_report_reg.add_record('Norm depth', best_norm_depth)
    #
    #     for s, r in zip(samples, nth_region_from_each_sample):
    #         comb_report_reg.add_record(s.name + ' hotspots depths/norm depths', r.records[9].value)
    #
    # best_report_basename = 'Best.' + source.targetseq_name  + '_' + vcf_key
    # comb_report_basename = 'Comb.' + source.targetseq_name  + '_' + vcf_key
    # # best_targetcov_norm_depth_vcf_txt = best_report.save_txt(output_dir, best_report_basename)
    # # best_targetcov_norm_depth_vcf_tsv = best_report.save_tsv(output_dir, best_report_basename)
    # comb_targetcov_norm_depth_vcf_txt = comb_report.save_txt(output_dir, comb_report_basename)
    # comb_targetcov_norm_depth_vcf_tsv = comb_report.save_tsv(output_dir, comb_report_basename)
    # info('')
    # info('Depths for Oncomine variants (total: {0:,} variants, {0:,} regions) saved into:'.format(total_variants))
    # # info('  Best:     ' + best_targetcov_norm_depth_vcf_txt)
    # info('  Combined: ' + comb_targetcov_norm_depth_vcf_txt)
    #
    return None, None  # comb_targetcov_norm_depth_vcf_txt


def _get_targqc_metric(metric, header_metric_storage, report_type='targetcov'):  # report type is in ['targetcov', 'qualimap', 'ngscat']
    qualimap_to_targetcov_dict = {
        'Number of reads': header_metric_storage.find_metric('Reads'),
        'Mapped reads': header_metric_storage.find_metric('Mapped reads'),
        'Unmapped reads': header_metric_storage.find_metric('Unmapped reads'),
        'Mapped reads (on target)': header_metric_storage.find_metric('Reads mapped on target'),
        'Coverage Mean': header_metric_storage.find_metric('Average target coverage depth'),
        'Coverage Standard Deviation': header_metric_storage.find_metric('Std. dev. of target coverage depth')
    }

    ngscat_to_targetcov_dict = {
        'Number reads': header_metric_storage.find_metric('Mapped reads'),
        # '% target bases with coverage >= 1x': cov.header_metric_storage.get_metric('Percentage of target covered by at least 1 read'),
        '% reads on target': header_metric_storage.find_metric('Reads mapped on target'),
        'mean coverage': header_metric_storage.find_metric('Average target coverage depth')
    }

    if report_type == 'targetcov':
        return metric
    elif report_type == 'qualimap':
        if metric.name in qualimap_to_targetcov_dict:
            return qualimap_to_targetcov_dict[metric.name]
        return metric
    elif report_type == 'ngscat':
        if metric.name in ngscat_to_targetcov_dict:
            return ngscat_to_targetcov_dict[metric.name]
        return metric
    # critical('Incorrect usage of get_targqc_metric(), report_type is %s but should be one of the following: %s' %
    #          (report_type, ", ".join(['targetcov', 'qualimap', 'ngscat'])))
    return metric


def _get_targqc_metric_storage(metric_storages_by_report_type):
    class SectionId:
        def __init__(self, name, title):
            self.name = name
            self.title = title

        def __hash__(self):
            #return hash((self.name, self.title))
            return hash(self.name)  # use title from the first metric_storage

        def __eq__(self, other):
            #return (self.name, self.title) == (other.name, other.title)
            return self.name == other.name  # use title from the first metric_storage

    metrics_by_sections = OrderedDict()
    general_section_id = None
    general_section_metric_list = []

    for report_type, metric_storage in metric_storages_by_report_type:
        for section in metric_storage.sections:
            section_id = SectionId(section.name, section.title)
            if section_id not in metrics_by_sections.keys():
                metrics_by_sections[section_id] = []

            metrics_by_sections[section_id] += [metric
                for metric in metric_storage.get_metrics(sections=[section], skip_general_section=True)
                if metric == _get_targqc_metric(metric, dict(metric_storages_by_report_type)['targetcov'], report_type)]

        # specific behaviour for general section
        general_section_metric_list += [metric
            for metric in metric_storage.general_section.metrics
            if metric == _get_targqc_metric(metric, dict(metric_storages_by_report_type)['targetcov'], report_type)]
        if not general_section_id:
            general_section_id = SectionId(metric_storage.general_section.name, metric_storage.general_section.title)

    sections = []
    for section_id, metric_list in metrics_by_sections.items():
        sections.append(ReportSection(section_id.name, section_id.title, metric_list))

    return MetricStorage(
        general_section=ReportSection(
            general_section_id.name, general_section_id.title, general_section_metric_list),
        sections=sections)


def _get_targqc_records(records_by_report_type, header_storage):
    targqc_records = []
    filled_metric_names = []
    for report_type, records in records_by_report_type:
        for record in records:
            new_metric = _get_targqc_metric(record.metric, header_storage, report_type)
            if not new_metric or new_metric.name not in filled_metric_names:
                filled_metric_names.append(record.metric.name)
                record.metric = new_metric
                targqc_records.append(record)
    return targqc_records


def _correct_qualimap_genome_results(cnf, samples):
    """ fixing java.lang.Double.parseDouble error on entries like "6,082.49"
    """
    for s in samples:
        if verify_file(s.qualimap_genome_results_fpath):
            correction_is_needed = False
            with open(s.qualimap_genome_results_fpath, 'r') as f:
                content = f.readlines()
                metrics_started = False
                for line in content:
                    if ">> Reference" in line:
                        metrics_started = True
                    if metrics_started:
                        if line.find(',') != -1:
                            correction_is_needed = True
                            break
            if correction_is_needed:
                with open(s.qualimap_genome_results_fpath, 'w') as f:
                    metrics_started = False
                    for line in content:
                        if ">> Reference" in line:
                            metrics_started = True
                        if metrics_started:
                            if line.find(',') != -1:
                                line = line.replace(',', '')
                        f.write(line)


def _correct_qualimap_insert_size_histogram(cnf, samples):
    """ replacing Qualimap insert size histogram with Picard one.
    """
    for s in samples:
        qualimap1_dirname = dirname(s.qualimap_ins_size_hist_fpath).replace('raw_data_qualimapReport', 'raw_data')
        qualimap2_dirname = dirname(s.qualimap_ins_size_hist_fpath)
        if exists(qualimap1_dirname):
            if not exists(qualimap2_dirname):
                shutil.move(qualimap1_dirname, qualimap2_dirname)
            else:
                shutil.rmtree(qualimap1_dirname)
        elif not exists(qualimap2_dirname):
            continue  # no data from both Qualimap v.1 and Qualimap v.2

        # if qualimap histogram exits and reuse_intermediate, skip
        if verify_file(s.qualimap_ins_size_hist_fpath, silent=True) and cnf.reuse_intermediate:
            pass
        else:
            if verify_file(s.picard_ins_size_hist_txt_fpath):
                with open(s.picard_ins_size_hist_txt_fpath, 'r') as picard_f:
                    one_line_to_stop = False
                    for line in picard_f:
                        if one_line_to_stop:
                            break
                        if line.startswith('## HISTOGRAM'):
                            one_line_to_stop = True

                    with file_transaction(cnf.work_dir, s.qualimap_ins_size_hist_fpath) as tx:
                        with open(tx, 'w') as qualimap_f:
                            for line in picard_f:
                                qualimap_f.write(line)


def select_best(values, fn=max):
    vs = [v for v in values if v is not None]
    return fn(vs) if len(vs) > 0 else None


def get_int_val(v):
    v = _get_num(v)
    return int(v) if v else None

def get_float_val(v):
    v = _get_num(v)
    return float(v) if v else None

def _get_num(v):
    v = get_val(v)
    return ''.join(c for c in v if c.isdigit() or c == '.') if v else None

def get_val(v):
    return v.strip() if v.strip() not in ['.', '-', ''] else None


def _save_best_details_for_each_gene(depth_threshs, samples, output_dir):
    metric_storage = get_detailed_metric_storage(depth_threshs)

    report = PerRegionSampleReport(sample='Best', metric_storage=metric_storage)
    report.add_record('Sample', 'contains best values from all samples: ' + ', '.join([s.name for s in samples]))

    total_regions = 0
    fpaths = [s.targetcov_detailed_tsv for s in samples if verify_file(s.targetcov_detailed_tsv)]
    if not fpaths:
        err('No targetcov detailed per-gene report was generated; skipping.')
        return None

    open_tsv_files = [open(fpath) for fpath in fpaths]

    first_col = 0
    while True:
        lines_for_each_sample = [next(f, None) for f in open_tsv_files]
        if not all(lines_for_each_sample):
            break
        l = lines_for_each_sample[0]
        if l.startswith('##'):
            continue
        elif l.startswith('#'):
            if l.startswith('#Sample'):
                first_col = 1
            break

    while True:
        lines_for_each_sample = [next(f, None) for f in open_tsv_files]
        if not all(lines_for_each_sample):
            break

        if all([not l.startswith('#') and ('Whole-Gene' in l or 'Gene-Exon' in l) for l in lines_for_each_sample]):
            shared_fields = lines_for_each_sample[0].split('\t')[first_col:first_col+8]
            reg = report.add_row()
            reg.add_record('Chr', get_val(shared_fields[0]))
            reg.add_record('Start', get_int_val(shared_fields[1]))
            reg.add_record('End', get_int_val(shared_fields[2]))
            reg.add_record('Size', get_int_val(shared_fields[3]))
            reg.add_record('Gene', get_val(shared_fields[4]))
            reg.add_record('Strand', get_val(shared_fields[5]))
            reg.add_record('Feature', get_val(shared_fields[6]))
            reg.add_record('Biotype', get_val(shared_fields[7]))

            min_depths, ave_depths, stddevs, withins = ([], [], [], [])
            percents_by_threshs = {t: [] for t in depth_threshs}

            for l in lines_for_each_sample:
                fs = l.split('\t')

                min_depths.append(get_int_val(fs[first_col+8]))
                ave_depths.append(get_float_val(fs[first_col+9]))
                stddevs.append(get_float_val(fs[first_col+10]))
                withins.append(get_float_val(fs[first_col+11]))
                for t, f in zip(depth_threshs, fs[first_col+12:]):
                    percents_by_threshs[t].append(get_float_val(f))

            # counting bests
            reg.add_record('Min depth', select_best(min_depths))
            reg.add_record('Ave depth', select_best(ave_depths))
            reg.add_record('Std dev', select_best(stddevs, max))
            reg.add_record('W/n 20% of ave depth', select_best(withins))
            for t in depth_threshs:
                reg.add_record('{}x'.format(t), select_best(percents_by_threshs[t]))

            total_regions += 1

    for f in open_tsv_files:
        f.close()

    gene_report_basename = 'Best.' + source.targetseq_name + source.detail_gene_report_baseending
    txt_rep_fpath = report.save_txt(join(output_dir, gene_report_basename + '.txt'))
    tsv_rep_fpath = report.save_tsv(join(output_dir, gene_report_basename + '.tsv'))
    info('')
    info('Best values for the regions (total ' + str(total_regions) + ') saved into:')
    info('  ' + txt_rep_fpath)

    return txt_rep_fpath


def get_bed_targqc_inputs(cnf, bed_fpath=None):
    if bed_fpath:
        bed_fpath = verify_bed(bed_fpath, description='input bed file', is_critical=True)

    exons_bed_fpath = adjust_path(cnf.exons if cnf.exons else cnf.genome.exons)
    if exons_bed_fpath:
        info('Exons: ' + exons_bed_fpath)

    if bed_fpath:
        info('Using amplicons/capture panel ' + bed_fpath)

    genes_fpath = None
    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)

    return bed_fpath, exons_bed_fpath, genes_fpath

