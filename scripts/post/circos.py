#!/usr/bin/env python
from genericpath import exists

import bcbio_postproc
from optparse import OptionParser
import vcf
from os.path import join

from source.calling_process import call
from source.clinical_reporting.clinical_parser import get_key_or_target_bed_genes
from source.config import Config
from source.file_utils import file_transaction
from source.logger import critical
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, determine_run_cnf, \
    determine_sys_cnf
from source.tools_from_cnf import get_system_path
from source.utils import is_us, is_uk

_circos_R_script = """#!/usr/bin/env Rscript

library(circlize) # call required R package - installed on the US (R 3.2.2) and UK (R 3.2.0) HPC

# Read into R
# seq2c
# look up these fields
seq2cfields = c("Chr","Start","End","Log2ratio","Amp_Del","KeyGene")
# read in header
seq2chead = strsplit(readLines(pipe("head -n1 {seq2c}")), '\\t')
# figure out indices to required fields
cutfields = paste(which(seq2chead[[1]] %in% seq2cfields), collapse = ",")
# read in the values for the given sample
sample.seq2c.bed = read.table(pipe(sprintf('cut -f %s -d "\\t" {seq2c} | sort -V | uniq', cutfields)) , sep = "\\t" )
sample.seq2c.bed = sample.seq2c.bed[-1,]
colnames(sample.seq2c.bed) = seq2cfields
sample.seq2c.bed$Log2ratio <- as.numeric(as.character(sample.seq2c.bed$Log2ratio))
sample.seq2c.bed$Start <- as.numeric(as.character(sample.seq2c.bed$Start))
sample.seq2c.bed$End <- as.numeric(as.character(sample.seq2c.bed$End))

# dbsnp
dbsnpfields = c("Chr","Start","AlleleFreq")
# read in header
dbsnphead = strsplit(readLines(pipe("head -n1 {vardict_all}")), '\\t')
# figure out indices to required fields
cutfieldsdbsnp = paste(which(dbsnphead[[1]] %in% dbsnpfields), collapse = ",")

sample.dbSNPs.bed = read.table(pipe(sprintf('grep -Fw {mysample} {vardict_all} | cut -f %s | sort -V | uniq ', cutfieldsdbsnp)) , sep = "\\t" )

sample.dbSNPs.bed = cbind(sample.dbSNPs.bed[,1],sample.dbSNPs.bed[,2]-1, sample.dbSNPs.bed[,c(2,3)])
colnames(sample.dbSNPs.bed) = c("Chr","Start","End","AlleleFreq")

# SVS
svloci = read.table(pipe('grep -vE "^(GL)|_|\\\\*" {svsbed}'), sep = "\\t") # only use std chromosomes
bed1 = data.frame(chr=svloci[,1], start=svloci[,2]-1, end=svloci[,2]-1, value=rep(1,times=nrow(svloci)))
bed2 = data.frame(chr=svloci[,3], start=svloci[,4]-1, end=svloci[,4]-1, value=rep(1,times=nrow(svloci)))
# Set up the plot

circos.clear() # best to always call this so as to clear the palette
png("{outdir}/{mysample}.png", res=300,units='cm',width=42, height=42) # produces high-res png which can be zoomed into
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
cytoband = "{cytoband}"
if (cytoband == "hg19"){{
    circos.initializeWithIdeogram() # uses hg19 as default cytoband
}}else{{
    circos.initializeWithIdeogram(cytoband = "{cytoband}") # uses hg19 as default cytoband
}}

om = circos.par("track.margin")
oc = circos.par("cell.padding")
circos.par(track.margin = om, cell.padding = oc)

# plot copy number
circos.genomicTrackPlotRegion(track.height=.2, sample.seq2c.bed, numeric.column = 4, panel.fun = function(region, value, ...)
{{
    xlim = get.cell.meta.data("xlim")
    circos.lines(xlim, c(0, 0), col = "#00000040")
    ampdel = value[[2]]
    keygene = value[[3]]
    value = value[[1]]

    circos.genomicPoints(subset(region, ampdel == 1 & keygene == 1), subset(value, ampdel == 1 & keygene == 1), pch = 16, cex = 0.5, col = "red")
    circos.genomicPoints(subset(region, ampdel == 1 & keygene == 0), subset(value, ampdel == 1 & keygene == 0), pch = 16, cex = 0.3, col = "#FF8080")
    circos.genomicPoints(subset(region, ampdel == -1 & keygene == 1), subset(value, ampdel == -1 & keygene == 1), pch = 16, cex = 0.5, col = "blue")
    circos.genomicPoints(subset(region, ampdel == -1 & keygene == 0), subset(value, ampdel == -1 & keygene == 0), pch = 16, cex = 0.3, col = "#67ddff")
    circos.genomicPoints(subset(region, ampdel == 0 & keygene == 1), subset(value, ampdel == 0 & keygene == 1), pch = 16, cex = 0.4, col = "black")
    circos.genomicPoints(subset(region, ampdel == 0 & keygene == 0), subset(value, ampdel == 0 & keygene == 0), pch = 16, cex = 0.2, col = "gray")
}}
)

# plot SNPs
circos.genomicTrackPlotRegion(track.height=.2, sample.dbSNPs.bed, ylim=c(0,1), panel.fun = function(region, value, ...)
{{
    cell.xlim = get.cell.meta.data("cell.xlim")
    for(h in c(0, 0.5)) {{
       circos.lines(cell.xlim, c(h, h), col = "#00000040")
    }}
    circos.genomicPoints(region, value, pch = 16, cex = 0.1)
}}
)

# plot SVs
circos.genomicLink(bed1, bed2, col = sample(1:5, nrow(bed1), replace = TRUE), border = NA)

text(0,0,"{mysample}", font=2, cex=4) # prints name of sample in the centre of the figure
dev.off()

"""

def get_args():
    description = (
        'Plots a Circos plot given vardict variant file (with all dbSNP SNPs, not the PASS one), '
        'Seq2C CNV calls and Manta SVs.')
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)
    parser.add_option('--bed', dest='bed_fpath', help='Path to BED file')
    parser.add_option('-v', '--mutations', dest='mutations_fpath', help='Path to VarDict.txt file')
    parser.add_option('-c', '--seq2c', dest='seq2c_tsv_fpath', help='Path to seq2c copy number file')
    parser.add_option('--sv', dest='sv_fpath', help='Path to Manta SV call vcf.gz file')
    parser.add_option('-s', '--sample', dest='sample', help='Identifier of sample in VarDict and Seq2c files')
    parser.add_option('-o', '--output-dir', dest='output_dir', default="./",
                        help='Output directory. Defaults to ./')
    (opts, args) = parser.parse_args()
    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    return cnf


def modify_seq2c(cnf, key_genes_chrom, seq2c_file, out_seq2c_fpath):
    key_genes = set(g for (g,c) in key_genes_chrom)
    with open(seq2c_file) as seq2c:
        with file_transaction(cnf.work_dir, out_seq2c_fpath) as tx:
            with open(tx, 'w') as out:
                for i, l in enumerate(seq2c):
                    if i == 0:
                        out.write('\t'.join(['Chr', 'Start', 'End', 'Log2ratio', 'Amp_Del', 'KeyGene']) + '\n')
                        continue
                    fs = l.replace('\n', '').split('\t')
                    sname, gname = fs[0], fs[1]
                    if sname != cnf.sample: continue
                    if 'not_a_gene' in gname: continue

                    sname, gname, chrom, start, end, length, log2r, sig, fragment, amp_del, ab_seg, total_seg, \
                        ab_log2r, log2r_diff, ab_seg_loc, ab_samples, ab_samples_pcnt = fs[:17]
                    is_key_gene = 0
                    if gname in key_genes:
                        is_key_gene = 1
                    is_amp_del = 0
                    if amp_del == 'Amp':
                        is_amp_del = 1
                    elif amp_del == 'Del':
                        is_amp_del = -1
                    out.write('\t'.join([chrom, start, end, str(ab_log2r) if ab_log2r else str(log2r),
                                         str(is_amp_del), str(is_key_gene)]) + '\n')



def parse_svs(cnf, sv_file, out_bed_fpath):
    """
    Parse sv vcf into a bed file
    """
    bp_dict = {}
    vcf_reader = vcf.Reader(filename=sv_file)
    with file_transaction(cnf.work_dir, out_bed_fpath) as tx:
        with open(tx, 'w') as out:
            for record in vcf_reader:
                if record.FILTER is None:
                    record.FILTER = []
                try: # if there is no SVTYPE or MATEID, ignore for now
                    if record.INFO['SVTYPE'] == "BND":
                        if record.INFO['MATEID'][0] not in bp_dict:
                            #Store the record in the dict until its pair is found
                            bp_dict[record.ID] = record
                        else:
                            #If the other BND is in the dict, annotate
                            record2 = bp_dict[record.INFO['MATEID'][0]]
                            try:
                                if record.samples[0]["PR"][1] + record.samples[0]["SR"][1] >= 5:
                                    out.write('\t'.join([str(record.CHROM), str(record.POS), str(record2.CHROM), str(record2.POS) + '\n']))
                            except AttributeError:
                                pass
                            #remove used record from the dict
                            del bp_dict[record.INFO['MATEID'][0]]
                    else:
                        #first check if 'END' is specified
                        if 'END' in record.INFO:
                            try:
                                # require 1Mb difference in coordinates and evidence from 10 reads or more
                                if record.samples[0]["PR"][1] + record.samples[0]["SR"][1] >= 10 and abs(record.INFO['END'] - record.POS) > 1000000:
                                    out.write('\t'.join([str(record.CHROM), str(record.POS), str(record.CHROM), str(record.INFO['END']) + '\n']))
                            except AttributeError:
                                pass
                except KeyError:
                    pass
                    #record.FILTER.append("REJECT")
                    #vcf_writer.write_record(record)


def main():
    cnf = get_args()

    vardict_res_fpath = cnf.mutations_fpath
    seq2c_tsv_fpath = cnf.seq2c_tsv_fpath
    sv_fpath = cnf.sv_fpath
    output_dir = cnf.output_dir
    sample_name = cnf.sample

    cytoband = None
    if 'hg38' in cnf.genome.name:
        cytoband = cnf.genome.circos_cytoband
    elif cnf.genome.name == 'hg19':
        cytoband = 'hg19'
    if not cytoband:
        critical('Circos plot does not support ' + cnf.genome + ' genome')

    svs_bed_fpath = join(output_dir, sample_name + '_svs.bed')
    parse_svs(cnf, sv_fpath, svs_bed_fpath)
    if not exists(svs_bed_fpath):
        return None

    modified_seq2c_fpath = join(output_dir, sample_name + '_seq2c.tsv')

    key_genes_chrom, _ = get_key_or_target_bed_genes(cnf.bed_fpath, cnf.key_genes)
    modify_seq2c(cnf, key_genes_chrom, seq2c_tsv_fpath, modified_seq2c_fpath)

    out_r_script = join(output_dir, sample_name + '.R')

    with open(out_r_script, 'w') as out_r_script_handle:
        out_r_script_handle.write(_circos_R_script.format(vardict_all=vardict_res_fpath,
                                  seq2c=modified_seq2c_fpath, outdir=output_dir,
                                  cytoband=cytoband, mysample=sample_name,
                                  svsbed=svs_bed_fpath))

    r_script = get_system_path(cnf, 'Rscript')
    if is_us():
        r_script = '/opt/az/local/R/R-3.2.2/installdir/bin/Rscript'
    elif is_uk():
        r_script = '/apps/R/3.2.0/rhel6-x64/bin/Rscript'
    cmdline = '{r_script} {out_r_script}'.format(**locals())
    res = call(cnf, cmdline)
    #os.remove(modified_seq2c_fpath)
    #os.remove(svs_bed_fpath)
    #os.remove(out_r_script)


if __name__ == "__main__":
    main()