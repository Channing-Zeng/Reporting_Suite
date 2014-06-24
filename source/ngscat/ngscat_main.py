#!/usr/bin/python

import glob
import time
import sys
import os
import string
import multiprocessing
import shutil
import math

import numpy
import pysam
import coverageHisto
import coverage_target
import coveragecorr
import bam_file
import bed_file
import bedgraph_file
import target_coverage
import coverage_saturation
import exon_coverage_std
import gcbias
import config

from source.file_utils import intermediate_fname
from source.logger import step_greetings
from os.path import join


def launch_coveragebed(bamfilenames, bedfilename, legend, outdir, executiongranted):
    coveragefiles = []
    Pcoveragebeds = []
    cnf = config.cnf

    for i, bamfilename in enumerate(bamfilenames):
        coveragefile = join(cnf["work_dir"], os.path.basename(bamfilename)[:-len('bam')] + 'coverage')
        coveragebedgraph = join(outdir, 'data', legend[i][:-len('bam')] + 'bed')

        print 'Coveragefile = ' + coveragefile
        bam = bam_file.BamFile(bamfilename, 'rb')

        print 'Launching coverageBed...'
        Pcoveragebed = multiprocessing.Process(target=bam.myCoverageBed, args=(
        bedfilename, None, coveragefile, executiongranted, coveragebedgraph,))
        Pcoveragebed.start()

        print '    Done.'

        coveragefiles.append(coveragefile)
        Pcoveragebeds.append(Pcoveragebed)

    #    return [positions,coverage,chromosomes,processedbed]
    return [Pcoveragebeds, coveragefiles]


def launch_covered_positions(coveragefiles, coveragethresholds, outdir, legend, executiongranted):
    status = multiprocessing.Value('b', False)
    coveredbases = multiprocessing.Array('f', len(coveragefiles))

    #    print 'outdir = '+outdir
    #    print 'coveragefile = '+coveragefile
    Pcoveredpositions = multiprocessing.Process(target=target_coverage.target_coverage_lite,
                                                args=(coveragefiles, coveragethresholds, outdir, legend, None,
                                                      executiongranted, status, coveredbases, config.warnbasescovered,))

    print 'Launching covered positions calculation...'
    Pcoveredpositions.start()

    return Pcoveredpositions, status, coveredbases


def launch_coverage_saturation(bamfilenames, bedfilename, depthlist, legend, outdir, executiongranted):
    status = multiprocessing.Value('b', False)
    saturationslopes = multiprocessing.Array('f', len(bamfilenames))
    Psaturation = multiprocessing.Process(target=coverage_saturation.coverage_saturation_lite,
                                          args=(bamfilenames, [bedfilename for i in range(len(bamfilenames))],
                                                depthlist, 10, legend, outdir + '/coverage_saturation_10x.png',
                                                executiongranted, status, saturationslopes,
                                                config.warnsaturation,))
    print 'Launching coverage saturation calculation...'
    Psaturation.start()

    return Psaturation, status, saturationslopes


def launch_coverage_distribution(coveragefiles, outdir, legend, executiongranted):
    status = multiprocessing.Value('b', False)
    meancoverage = multiprocessing.Array('f', len(coveragefiles))

    Pcoveragedistribution = multiprocessing.Process(target=coverageHisto.histo_CV, args=(
    coveragefiles, outdir, legend, executiongranted, status, meancoverage,
    config.warnmeancoverage,))
    print 'Launching coverage distribution calculation...'
    Pcoveragedistribution.start()

    return Pcoveragedistribution, status, meancoverage


def launch_coveragecorr(coveragefiles, fileout, legend, executiongranted):
    status = multiprocessing.Value('b', False)
    corr = multiprocessing.Value('f', 0)

    Pcoveragecorr = multiprocessing.Process(target=coveragecorr.coveragecorr,
                                            args=(coveragefiles, fileout, legend, executiongranted, status, corr,
                                                  config.warncoveragecorrelation,))
    print 'Launching coverage correlation calculation...'
    Pcoveragecorr.start()

    return Pcoveragecorr, status, corr


def launch_onoff_reads(bamfilenames, bedfilename, legend, outdir, executiongranted):
    onoff_status = multiprocessing.Value('b', False)
    duplicates_status = multiprocessing.Value('b', False)
    enrichment = multiprocessing.Array('f', len(bamfilenames))
    percontarget = multiprocessing.Array('f', len(bamfilenames))
    onduplicates = multiprocessing.Array('f', len(bamfilenames))
    offduplicates = multiprocessing.Array('f', len(bamfilenames))

    bam = bam_file.BamFile(bamfilenames[0], 'rb')
    print 'Launching on/off target enrichment calculation...'
    Ponoff_reads = multiprocessing.Process(target=bam.reads_on_target, args=(
    bedfilename, outdir, [bam_file.BamFile(bamfilenames[i]) for i in range(1, len(bamfilenames))],
    legend, executiongranted, onoff_status, duplicates_status, onduplicates,
    offduplicates, enrichment, percontarget, config.warnontarget,))
    Ponoff_reads.start()
    bam.close()

    return Ponoff_reads, onoff_status, onduplicates, offduplicates, duplicates_status, enrichment, percontarget


def sequential_offclusters_call(offtargetoffset, offtargetthreshold, bedgraphfilenames, bedfilename, executiongranted):
    if (executiongranted <> None):
        executiongranted.acquire()

    for bedgraphfilename in bedgraphfilenames:
        bedgraph = bedgraph_file.bedgraph_file(bedgraphfilename)
        outfile = bedgraphfilename.replace('.bed', '.off.bed') if '.off' not in bedgraphfilename else bedgraphfilename
        bedgraph.getOffTarget(offtargetoffset, offtargetthreshold, bedfilename, outfile)

    if (executiongranted <> None):
        executiongranted.release()


def launch_offclusters(bedgraphfilenames, bedfilename, executiongranted):
    ### sorting all at once and reuse in the future
    bed = bed_file.BedFile(bedfilename)
    extendedBed = bed.extend(config.offtargetoffset)
    sortedBed = extendedBed.sort_bed()
    nonOverlappingBed = sortedBed.non_overlapping_exons(-1)  # Base 0, it is a standard BED
    nonOverlappingBed.sort_bed()
    ###

    Poffclusters = multiprocessing.Process(target=sequential_offclusters_call,
                                           args=(config.offtargetoffset, config.offtargetthreshold, bedgraphfilenames,
                                                 bedfilename, executiongranted,))
    print 'Launching off target clusters...'
    Poffclusters.start()

    return Poffclusters


def launch_coverage_through_target(coveragefiles, outdir, legend, executiongranted):
    status = multiprocessing.Value('b', False)
    lowcovbases = multiprocessing.Array('f', len(coveragefiles))

    if (len(coveragefiles) < 2):
        legend = None

    Pcoveragethroughtarget = multiprocessing.Process(target=coverage_target.print_coverage, args=(
    coveragefiles, 1000, outdir, legend, executiongranted, status,
    lowcovbases,
    config.warncoverageregion, config.warncoveragethreshold,))
    print 'Launching coverage throught target calculation...'
    Pcoveragethroughtarget.start()

    return Pcoveragethroughtarget, status, lowcovbases


def launch_coverage_std(coveragefiles, outdir, legend, executiongranted):
    status = multiprocessing.Value('b', False)
    coveragestd = multiprocessing.Array('f', len(coveragefiles))

    Pcoveragestd = multiprocessing.Process(target=exon_coverage_std.exon_coverage_std_lite,
                                           args=([[coveragefile] for coveragefile in coveragefiles], outdir,
                                                 legend, executiongranted, status, coveragestd, config.warnstd,))

    print 'Launching coverage std calculation...'
    Pcoveragestd.start()

    return Pcoveragestd, status, coveragestd


def launch_gcbias(coveragefile, bedfilename, reference, fileout, graphtitle, executiongranted):
    status = multiprocessing.Value('b', False)

    Pgcbias = multiprocessing.Process(target=gcbias.gcbias_lite, args=(
    coveragefile, bedfilename, reference, fileout, graphtitle, executiongranted, status,))

    print 'Launching coverage std calculation...'
    Pgcbias.start()

    return Pgcbias, status


def generate_report(cnf, bamfilenames, sortedbams, bedfilename, outdir, coveredpositions_status, coveredbases,
                    coverage_saturation_status, saturationslopes, onoff_status,
                    duplicates_status, onduplicates, offduplicates, coveragedistribution_status, meancoverage,
                    coveragecorr_status, corr, throughtarget_status, lowcovbases, coveragestd_status, coveragestd,
                    gcbias_status, enrichment, percontarget,
                    reference, nthreads,
                    depthlist,
                    coveragethresholds):
    dest = join(outdir, 'img')
    shutil.copy(join(config.IMGSRC, 'xls_icon.png'), dest)
    shutil.copy(join(config.IMGSRC, 'txt_icon.png'), dest)
    shutil.copy(join(config.IMGSRC, 'ok.jpg'), dest)
    shutil.copy(join(config.IMGSRC, 'warning.jpg'), dest)
    shutil.copy(join(config.IMGSRC, 'coverage_histogram_example.png'), dest)

    shutil.copy(join(config.DATASRC, 'styles.css'), outdir)

    # ********************************************************* INput parameters ******************************************************************
    if (coverage_saturation_status <> None):
        saturationcurve = 'Yes'
    else:
        saturationcurve = 'No'

    fd = open(join(config.DATASRC, 'captureQC.html'))
    reportcontent = string.join(fd.readlines(), sep='').replace('bamfilename',
                                                                string.join(bamfilenames, sep=', ')).replace(
        'bedfilename', bedfilename).replace('reportdate', time.ctime()).replace('reference', str(reference)).replace(
        'saturationcurve', saturationcurve).replace('nthreads', str(nthreads)).replace('tmpdir', cnf['work_dir'])
    fd.close()

    # ********************************************************* Result summary ******************************************************************

    jsonstr = ''
    for i, bam in enumerate(bamfilenames):
        jsonstr += '{"bamfile":"' + bam + '"'
        jsonstr += ',"nreads":' + str(bam_file.BamFile(sortedbams[i]).nreads())
        jsonstr += ',"coveredbases":' + str(coveredbases[i])

        if (coverage_saturation_status <> None):
            jsonstr += ',"saturationslope":' + str(saturationslopes[i])

        jsonstr += ',"percontarget":' + str(percontarget[i])
        jsonstr += ',"onduplicates":' + str(onduplicates[i])
        jsonstr += ',"offduplicates":' + str(offduplicates[i])
        jsonstr += ',"meancoverage":' + str(meancoverage[i])
        jsonstr += ',"lowcovbases":' + str(lowcovbases[i])

        if (not math.isnan(coveragestd[i])):
            jsonstr += ',"coveragestd":' + str(coveragestd[i]) + '}'
        else:
            jsonstr += '}'

    fd = file(outdir + '/data/summary.json', 'w')
    fd.write(jsonstr)
    fd.close()

    summaryrows = ''
    for i, bam in enumerate(bamfilenames):
        summaryrows += '<tr>\n'
        summaryrows += '<td class="table-cell"> ' + bam + '</td>'
        summaryrows += '<td class="table-cell"> ' + str(bam_file.BamFile(sortedbams[i]).nreads()) + ' </td>'
        summaryrows += '<td class="table-cell">%.1f' % (coveredbases[i]) + '% </td>'

        if (coverage_saturation_status <> None):
            summaryrows += '<td class="table-cell">%.1e</td>\n' % saturationslopes[i]

        summaryrows += '<td class="table-cell">%.1f' % (percontarget[i]) + '% </td>\n'
        summaryrows += ('<td class="table-cell">ON-%.1f%%' % onduplicates[i]) + '; OFF: %.1f' % (
        offduplicates[i]) + '% </td>'
        summaryrows += '<td class="table-cell">%.1fx' % meancoverage[i] + '</td>\n'
        summaryrows += '<td class="table-cell">%d consecutive bases<br>with coverage <= <WARNCOVERAGETHRESHOLD></td>\n' % (
        lowcovbases[i])

        if (coveragecorr_status <> None):
            summaryrows += '<td class="table-cell">%.2f</td>\n' % corr.value

        summaryrows += '<td class="table-cell">%.2f</td>\n' % coveragestd[i]
        summaryrows += '</tr>\n'

    summarystatus = '<td class="table-header">Overall status</td>\n'
    summarystatus += '<td class="table-header"></td>\n'
    summarystatus += '<td class="table-header"><a href="#targetbases"><img src="img/<TARGETBASESSTATUS>.jpg" height=23px /></a></td>\n'
    if (coverage_saturation_status <> None):
        summarystatus += '<td class="table-header"><a href="#coveragesaturation"><img src="img/<COVERAGESATURATIONSTATUS>.jpg" height=23px /></a></td>\n'
    summarystatus += '<td class="table-header"><a href="#onoff"><img src="img/<ONOFFSTATUS>.jpg" height=23px /></a></td>\n'
    summarystatus += '<td class="table-header"><a href="#dup"><img src="img/<DUPSTATUS>.jpg" height=23px /></a></td>\n'
    summarystatus += '<td class="table-header"><a href="#distribution"><img src="img/<DISTRIBUTIONSTATUS>.jpg" height=23px /></a></td>\n'
    summarystatus += '<td class="table-header"><a href="#coveragethroughtarget"><img src="img/<COVERAGETHROUGHTARGETSTATUS>.jpg" height=23px /></a></td>\n'
    if (coveragecorr_status <> None):
        summarystatus += '<td class="table-header"><a href="#coveragecorr"><img src="img/<COVERAGECORRSTATUS>.jpg" height=23px /></a></td>\n'
    summarystatus += '<td class="table-header"><a href="#coveragestd"><img src="img/<COVERAGESTDSTATUS>.jpg" height=23px /></a></td>\n'

    reportcontent = reportcontent.replace('<SUMMARYROWS>', summaryrows)
    reportcontent = reportcontent.replace('<SUMMARYSTATUS>', summarystatus)

    if (coverage_saturation_status <> None):
        reportcontent = reportcontent.replace('<SUMMARYSATURATION>',
                                              '<td class="table-header"><a href="#coveragesaturation">Coverage saturation<br>(slope at the end of the curve)</a></td>')
    else:
        reportcontent = reportcontent.replace('<SUMMARYSATURATION>', '')

    if (coveragecorr_status <> None):
        reportcontent = reportcontent.replace('<SUMMARYCOVCORRELATION>',
                                              '<td class="table-header"><a href="#coveragecorr">Coverage correlation<br>per ROI</a></td>')
    else:
        reportcontent = reportcontent.replace('<SUMMARYCOVCORRELATION>', '')

    reportcontent = reportcontent.replace('<SUMMARYCOVERAGETHRS>', str(coveragethresholds[0]))
    reportcontent = reportcontent.replace('<SUMMARYTARGETSIZE>', str(bed_file.BedFile(bedfilename).size()))


    # ********************************************************* Detailed results ******************************************************************
    chromosomeimages = ''
    ontarget_coverage_files = glob.glob(outdir + '/data/*_Ontarget_Coverage.png')
    ontarget_coverage_files.sort()
    for afile in ontarget_coverage_files:
        chromosomeimages += '<a href="data/' + os.path.basename(
            afile) + '"><img style="width: 33%; float: left;" src="data/' + os.path.basename(afile) + '" /></a>'
    reportcontent = reportcontent.replace('<CHROMOSOMEIMAGES>', chromosomeimages)

    if (coveredpositions_status.value):
        reportcontent = reportcontent.replace('<TARGETBASESSTATUS>', 'ok')
    else:
        reportcontent = reportcontent.replace('<TARGETBASESSTATUS>', 'warning')
    reportcontent = reportcontent.replace('<WARNBASESCOVERED>', str(config.warnbasescovered))

    percentagestr = '\n<ul>'
    enrichmentstr = '\n<ul>'
    for i, bamfilename in enumerate(bamfilenames):
        percentagestr += '<li>' + bamfilename + ': %.1f' % (percontarget[i]) + '%</li>\n'
        enrichmentstr += '<li>' + bamfilename + ': %.1f' % (enrichment[i]) + '</li>\n'
    percentagestr += '</ul>'
    enrichmentstr += '</ul>'
    reportcontent = reportcontent.replace('<PERCENTAGEONTARGET>', percentagestr)
    reportcontent = reportcontent.replace('<ENRICHMENT>', enrichmentstr)

    reportcontent = reportcontent.replace('<WARNONTARGET>', str(config.warnontarget))
    if (onoff_status.value):
        reportcontent = reportcontent.replace('<ONOFFSTATUS>', 'ok')
    else:
        reportcontent = reportcontent.replace('<ONOFFSTATUS>', 'warning')

    duplicates_files = glob.glob(outdir + '/data/duplicates*.png')
    duplicates_files.sort()
    dupimages = ''
    for afile in duplicates_files:
        dupimages += '<img style="width: 50%; float: left;" src="data/' + os.path.basename(afile) + '" /></a>'
    reportcontent = reportcontent.replace('<DUPIMAGES>', dupimages)

    if (duplicates_status.value):
        reportcontent = reportcontent.replace('<DUPSTATUS>', 'ok')
    else:
        reportcontent = reportcontent.replace('<DUPSTATUS>', 'warning')

    reportcontent = reportcontent.replace('<WARNMEANCOVERAGE>', str(config.warnmeancoverage))
    if (coveragedistribution_status.value):
        reportcontent = reportcontent.replace('<DISTRIBUTIONSTATUS>', 'ok')
    else:
        reportcontent = reportcontent.replace('<DISTRIBUTIONSTATUS>', 'warning')

    if (coveragecorr_status <> None):
        fd = open(join(config.DATASRC, 'coveragecorr_content.html'))
        coveragecorr_content = string.join(fd.readlines(), sep='')
        fd.close()
        reportcontent = reportcontent.replace('<COVERAGECORRCONTENT>', coveragecorr_content)

        reportcontent = reportcontent.replace('<WARNCOVERAGECORRELATION>', str(config.warncoveragecorrelation))
        if (coveragecorr_status.value):
            reportcontent = reportcontent.replace('<COVERAGECORRSTATUS>', 'ok')
        else:
            reportcontent = reportcontent.replace('<COVERAGECORRSTATUS>', 'warning')
    else:
        reportcontent = reportcontent.replace('<COVERAGECORRCONTENT>', '\n')

    reportcontent = reportcontent.replace('<WARNCOVERAGEREGION>', str(config.warncoverageregion))
    reportcontent = reportcontent.replace('<WARNCOVERAGETHRESHOLD>', str(config.warncoveragethreshold))
    if (throughtarget_status.value):
        reportcontent = reportcontent.replace('<COVERAGETHROUGHTARGETSTATUS>', 'ok')
    else:
        reportcontent = reportcontent.replace('<COVERAGETHROUGHTARGETSTATUS>', 'warning')

    reportcontent = reportcontent.replace('<WARNSTD>', str(config.warnstd))
    if (coveragestd_status.value):
        reportcontent = reportcontent.replace('<COVERAGESTDSTATUS>', 'ok')
    else:
        reportcontent = reportcontent.replace('<COVERAGESTDSTATUS>', 'warning')

    if (coverage_saturation_status <> None):
        fd = open(join(config.DATASRC, 'saturation_content.html'))
        saturation_content = string.join(fd.readlines(), sep='')
        fd.close()
        reportcontent = reportcontent.replace('<SATURATIONCONTENT>', saturation_content).replace('<DEPTHLIST>',
                                                                                                 string.join(map(str,
                                                                                                                 depthlist[
                                                                                                                 :-1]),
                                                                                                             sep='x10<sup>6</sup>, ') + 'x10<sup>6</sup> and ' + str(
                                                                                                     depthlist[
                                                                                                         -1]) + 'x10<sup>6</sup>').replace(
            'depthlist', str(depthlist)[1:-1])
        reportcontent = reportcontent.replace('<WARNSATURATION>', str(config.warnsaturation))

        if (coverage_saturation_status.value):
            reportcontent = reportcontent.replace('<COVERAGESATURATIONSTATUS>', 'ok')
        else:
            reportcontent = reportcontent.replace('<COVERAGESATURATIONSTATUS>', 'warning')
    else:
        reportcontent = reportcontent.replace('<SATURATIONCONTENT>', '\n').replace('depthlist', 'None')

    reportcontent = reportcontent.replace('coveragethrs', string.join(map(str, coveragethresholds), sep=', '))

    if (gcbias_status <> None):
        fd = open(join(config.DATASRC, 'gcbias_content.html'))
        gcbias_content = string.join(fd.readlines(), sep='')
        fd.close()
        reportcontent = reportcontent.replace('<GCBIASCONTENT>', gcbias_content)

        gcbiasimages = ''
        for afile in glob.glob(outdir + '/data/gcbias*.png'):
            gcbiasimages += '<img style="width:40%" src="data/' + os.path.basename(afile) + '" />'
        reportcontent = reportcontent.replace('<GCBIASIMAGES>', gcbiasimages)

    else:
        reportcontent = reportcontent.replace('<GCBIASCONTENT>', '\n')

    report_fpath = join(outdir, 'captureQC.html')
    fd = open(report_fpath, 'w')
    fd.write(reportcontent)
    fd.close()
    return report_fpath


def ngscat(cnf, bamfilenames, originalbedfilename, outdir, reference=None, saturation=False, nthreads=2, extend=None,
           depthlist='auto', coveragethresholds=[1, 5, 10, 20, 30],
           onefeature=None):

    step_greetings('ngsCAT starting! Assesses capture performance in terms of sensibility, '
                   'specificity and uniformity of the coverage.\n')

    if (not (os.path.isdir(outdir) or os.path.islink(outdir))):
        print 'WARNING: ' + outdir + ' does not exist. Creating directory.'
        os.mkdir(outdir)

    if (not (os.path.isdir(outdir + '/data') or os.path.islink(outdir + '/data'))):
        print 'Creating ' + outdir + '/data'
        os.mkdir(outdir + '/data')

    if (not (os.path.isdir(outdir + '/img') or os.path.islink(outdir + '/img'))):
        print 'Creating ' + outdir + '/img'
        os.mkdir(outdir + '/img')

    sortedbams = []
    for bamfilename in bamfilenames:
        filelink = join(cnf["work_dir"], os.path.basename(bamfilename))
        try:
            os.symlink(bamfilename, filelink)
        except OSError:
            print 'WARNING: when trying to create a symbolic link at the temporary directory pointing to ' + bamfilename + ', a file named ' + filelink + ' was already found.'
            print '    Probably the temporary and origin directories are the same. The only problem this could cause is that the new index overwrites an existing one.'
            #print '    Continue (y/n)?'

            #goahead = raw_input()
            goahead = 'y'
            if (goahead == 'n' or goahead == 'N'):
                print 'Exiting...'
                sys.exit(1)
            elif (goahead <> 'y' and goahead <> 'Y'):
                print 'Unknown choice ' + goahead
                print 'Exiting...'
                sys.exit(1)

            if os.path.dirname(bamfilename) != os.path.dirname(filelink):
                os.remove(filelink)
                os.symlink(bamfilename, filelink)

        print 'Indexing...'
        pysam.index(filelink)
        print '    Done.'

        if (not bam_file.BamFile(filelink).issorted()):
            print 'WARNING: ' + bamfilename + ' is not sorted'
            print 'Sorting...'
            newsortedbam_fpath = intermediate_fname(cnf, bamfilename, 'sorted')
            sortedbams.append(newsortedbam_fpath)
            pysam.sort(filelink, newsortedbam_fpath[:-len('.bam')])
            print 'Indexing...'
            pysam.index(sortedbams[-1])
            print '    Done.'
        else:
            sortedbams.append(filelink)

    if (saturation and depthlist == 'auto'):
        maxdepth = max([bam_file.BamFile(bamfilename).nreads() for bamfilename in sortedbams])
        depthlist = numpy.arange(maxdepth / 5.0, maxdepth + (maxdepth / 5.0) - 1, maxdepth / 5.0)
        depthlist = depthlist / 1000000.0

    legend = [os.path.basename(bamfilename) for bamfilename in bamfilenames]
    executiongranted = multiprocessing.Semaphore(nthreads)

    if (extend <> None):
        bedfilename = bed_file.BedFile(originalbedfilename).extend(extend).filename
    else:
        bedfilename = originalbedfilename

    if onefeature == None or onefeature == 'specificity' or onefeature == 'gcbias':
        ### sorting all at once and reuse in the future
        bed = bed_file.BedFile(bedfilename)
        sortedBed = bed.sort_bed()
        nonOverlappingBed = sortedBed.non_overlapping_exons(1)  # This generates a BED file in base 1 (Non-standard BED)
        nonOverlappingBed.sort_bed()
        ###

    if (onefeature == None or onefeature <> 'saturation' or onefeature <> 'specificity'):
        Pcoveragebeds, coveragefiles = launch_coveragebed(sortedbams, bedfilename, legend, outdir, executiongranted)

    if ((saturation and onefeature == None) or onefeature == 'saturation'):
        Psaturation, coverage_saturation_status, saturationslopes = launch_coverage_saturation(sortedbams, bedfilename,
                                                                                               depthlist, legend,
                                                                                               outdir + '/data/',
                                                                                               executiongranted)
    else:
        coverage_saturation_status = None
        saturationslopes = None

    if (onefeature == None or onefeature == 'specificity'):
        Ponoff_reads, onoff_status, onduplicates, offduplicates, duplicates_status, enrichment, percontarget = launch_onoff_reads(
            sortedbams, bedfilename, legend, outdir + '/data/', executiongranted)

    for i in range(len(Pcoveragebeds)):
        Pcoveragebeds[i].join()
        Pcoveragebeds[i].terminate()

    if (onefeature == None or onefeature == 'specificity'):
        Poffclusters = launch_offclusters(glob.glob(outdir + '/data/*.bed'), bedfilename, executiongranted)

    if (onefeature == None or onefeature == 'coveragefreq'):
        Pcoveragedistribution, coveragedistribution_status, meancoverage = launch_coverage_distribution(coveragefiles,
                                                                                                        outdir + '/data/',
                                                                                                        legend,
                                                                                                        executiongranted)

    if (onefeature == None or onefeature == 'percbases'):
        Pcoveredpositions, coveredpositions_status, coveredbases = launch_covered_positions(coveragefiles,
                                                                                            coveragethresholds,
                                                                                            outdir + '/data/', legend,
                                                                                            executiongranted)

    if (onefeature == None or onefeature == 'coveragedistr'):
        Pcoveragethroughtarget, throughtarget_status, lowcovbases = launch_coverage_through_target(coveragefiles,
                                                                                                   outdir + '/data/',
                                                                                                   legend,
                                                                                                   executiongranted)

    if (len(coveragefiles) > 1 and (onefeature == None or onefeature == 'coveragecorr')):
        Pcoveragecorr, coveragecorr_status, corr = launch_coveragecorr(coveragefiles, outdir + '/data/coveragecorr.png',
                                                                       legend, executiongranted)
    else:
        coveragecorr_status = None
        corr = None

    if (onefeature == None or onefeature == 'coveragestd'):
        Pcoveragestd, coveragestd_status, coveragestd = launch_coverage_std(coveragefiles, outdir + '/data/', legend,
                                                                            executiongranted)

    if ((reference <> None and onefeature == None) or onefeature == 'gcbias'):
        Pgcbias = []
        for i, coveragefile in enumerate(coveragefiles):
            onePgcbias, gcbias_status = launch_gcbias(coveragefile, bedfilename, reference,
                                                      outdir + '/data/gcbias' + str(i) + '.png', legend[i],
                                                      executiongranted)
            Pgcbias.append(onePgcbias)
        for onePgcbias in Pgcbias:
            onePgcbias.join()
            onePgcbias.terminate()
    else:
        gcbias_status = None

    # LAUNCH BASIC STATS

    if ((saturation and onefeature == None) or onefeature == 'saturation'):
        Psaturation.join()
        Psaturation.terminate()

    if (onefeature == None or onefeature == 'coveragefreq'):
        Pcoveragedistribution.join()
        Pcoveragedistribution.terminate()

    if (onefeature == None or onefeature == 'percbases'):
        Pcoveredpositions.join()
        Pcoveredpositions.terminate()

    if (onefeature == None or onefeature == 'coveragedistr'):
        Pcoveragethroughtarget.join()
        Pcoveragethroughtarget.terminate()

    if (len(coveragefiles) > 1 and (onefeature == None or onefeature == 'coveragecorr')):
        Pcoveragecorr.join()
        Pcoveragecorr.terminate()

    if (onefeature == None or onefeature == 'coveragestd'):
        Pcoveragestd.join()
        Pcoveragestd.terminate()

    if (onefeature == None or onefeature == 'specificity'):
        Ponoff_reads.join()
        Ponoff_reads.terminate()

        Poffclusters.join()
        Poffclusters.terminate()

    #    if(onefeature==None or onefeature<>'saturation'):
    #        for coveragefile in coveragefiles:
    #            os.remove(coveragefile)

    if (onefeature == None):
        return generate_report(cnf, bamfilenames, sortedbams, originalbedfilename, outdir, coveredpositions_status, coveredbases,
                        coverage_saturation_status, saturationslopes,
                        onoff_status,
                        duplicates_status, onduplicates, offduplicates, coveragedistribution_status, meancoverage,
                        coveragecorr_status, corr, throughtarget_status, lowcovbases, coveragestd_status, coveragestd,
                        gcbias_status, enrichment, percontarget,
                        reference, nthreads, depthlist,
                        coveragethresholds)
    return None
