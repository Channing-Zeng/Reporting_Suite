import sys
import subprocess

cmdline = 'cat ' + sys.argv[1] + ' | ' \
          'perl /group/ngs/src/snpEff/snpEff3.5/scripts/vcfEffOnePerLine.pl | ' \
          'java -jar /group/ngs/src/snpEff/snpEff3.5/SnpSift.jar extractFields - ' \
          'CHROM POS ID CNT GMAF REF ALT QUAL FILTER TYPE '

subprocess.call(cmdline.split(), stdout=open('aaa.tsv', 'w'))
