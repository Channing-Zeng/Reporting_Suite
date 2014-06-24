from os.path import join, dirname, realpath

warnbasescovered = 90
warnsaturation = 1e-5
warnontarget = 80
warnmeancoverage = 40
warncoverageregion = 100
warncoveragethreshold = 6
warncoveragecorrelation = 0.90
warnstd = 0.3

offtargetoffset = 1000
offtargetthreshold = 15

chr_lengths = None

DATASRC = join(dirname(realpath(__file__)), 'html')
IMGSRC = join(dirname(realpath(__file__)), 'img')

cnf = None