
# This file was produced by plot-vcfstats, the command line was:
#	plot-vcfstats -s /Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_bcftools.report -p /Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/ --no-PDF
#
# Edit as necessary and recreate the plots by running
#	python /Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/plot.py 
#
# Title abbreviations:
# 	 0 .. sampl .. /Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1.vcf.gz
#

# Set to 1 to plot in PDF instead of PNG
pdf_plots = 1

# Set to 1 to use sample names for xticks instead of numeric sequential IDs
#   and adjust margins and font properties if necessary
sample_names   = 1
sample_margins = {'right':0.98, 'left':0.07, 'bottom':0.2}
sample_font    = {'rotation':45, 'ha':'right', 'fontsize':8}

if sample_names==0: sample_margins=(); sample_font=();

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

import numpy
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1: raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.get_size < window_len: return x
	if window_len<3: return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']: raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	s = numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
	if window == 'flat': # moving average
		w = numpy.ones(window_len,'d')
	else:
		w = eval('numpy.'+window+'(window_len)')
	y = numpy.convolve(w/w.sum(),s,mode='valid')
	return y[(window_len/2-1):-(window_len/2)]


dat0 = [
[ 0.000000, 10 ],
]

fig = plt.figure(figsize=(2*4.33070866141732,3.93700787401575*0.7))
ax1 = fig.add_subplot(111)
ax1.set_ylabel('Number of sites')
ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y') 
ax1.set_yscale('log')
ax1.set_xlabel('Non-reference allele frequency')
ax1.set_xlim(-0.05,1.05)

dat0 = [
[ 0.000000, 1 ],
]

fig = plt.figure(figsize=(2*4.33070866141732,3.93700787401575*0.7))
ax1 = fig.add_subplot(111)
ax1.set_ylabel('Number of sites')
ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y') 
ax1.set_yscale('log')
ax1.set_xlabel('Non-reference allele frequency')
ax1.set_xlim(-0.05,1.05)

dat = [
]
dat = [
	[ 0.000000, 10, 4.000000 ],
]
dat = []
with open('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/tstv_by_qual.0.dat', 'rb') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
ax1 = fig.add_subplot(111)
ax1.plot([row[1] for row in dat], [row[2] for row in dat], '^-', ms=3, mec='orange', color='orange')
ax1.set_ylabel('Ts/Tv',fontsize=10)
ax1.set_xlabel('Number of sites\n(sorted by QUAL, descending)',fontsize=10)
ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')
ax1.set_ylim(min(2,min(row[2] for row in dat))-0.3,0.3+max(2.2,max(row[2] for row in dat)))

plt.subplots_adjust(right=0.88,left=0.15,bottom=0.15)
plt.title('sampl')
plt.savefig('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/tstv_by_qual.0.png')
if pdf_plots: plt.savefig('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/tstv_by_qual.0.pdf')
plt.close()

dat = [
	[2,1],
]

fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
ax1 = fig.add_subplot(111)
ax1.bar([row[0]-0.5 for row in dat], [row[1] for row in dat], color='orange')# , edgecolor='orange')
ax1.set_xlabel('InDel Length')
ax1.set_ylabel('Count')
ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
ax1.set_xlim(-20,20)
plt.subplots_adjust(bottom=0.17)
plt.title('sampl')
plt.savefig('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/indels.0.png')
if pdf_plots: plt.savefig('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/indels.0.pdf')
plt.close()

dat = [
	[0,'A>C',1],
	[1,'A>G',1],
	[2,'A>T',1],
	[3,'C>A',0],
	[4,'C>G',0],
	[5,'C>T',3],
	[6,'G>A',1],
	[7,'G>C',0],
	[8,'G>T',0],
	[9,'T>A',0],
	[10,'T>C',3],
	[11,'T>G',0],
]

fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
cm  = mpl.cm.get_cmap('autumn')
n = 12
col = range(n)
for i in range(n): col[i] = cm(1.*i/n)
ax1 = fig.add_subplot(111)
ax1.bar([row[0] for row in dat], [row[2] for row in dat], color=col)
ax1.set_ylabel('Count')
ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
ax1.set_xlim(-0.5,n+0.5)
plt.xticks([row[0] for row in dat],[row[1] for row in dat],rotation=45)
plt.title('sampl')
plt.savefig('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/substitutions.0.png')
if pdf_plots: plt.savefig('/Users/vladsaveliev/vagrant/ngs_analysis/test/results_qc/work/sample1_qc_plots/substitutions.0.pdf')
plt.close()
