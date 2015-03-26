with open('__comparison_summary.txt') as f:
	txt=f.read().strip()

report = [l.split('\t')[1:] for l in txt.split('\n')]

def getval(s):
	return ''.join(c for c in s if c.isdigit() or c == '.')

samples = []
for sn in report[4]:
	if '-' in sn:
		sn = '-'.join(sn.split('-')[:-1])
	if '.' in sn:
		sn = '-'.join(sn.split('.')[:-1])
	samples.append(sn)

dup_rate = report[10]
reads_mapped = report[5]
reads_mapped_in_target = report[17]
depth1x = report[26]
depth10x = report[28]
depth25x = report[31]

print  '\t'.join(samples)
print  '\t'.join([getval(v) for v in dup_rate])
print  '\t'.join([getval(v) for v in reads_mapped])
print  '\t'.join([getval(v) for v in reads_mapped_in_target])
print ''
print  '\t'.join([getval(v) for v in depth1x])
print  '\t'.join([getval(v) for v in depth10x])
print  '\t'.join([getval(v) for v in depth25x])








