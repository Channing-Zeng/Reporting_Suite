with open('__comparison_summary.txt') as f:
    txt = f.read().strip()

report = [l.split('\t') for l in txt.split('\n') if l]

def getval(s):
    return ''.join(c for c in s if c.isdigit() or c == '.')

samples = []
for sn in report[4][1:]:
    if '-' in sn:
        sn = '-'.join(sn.split('-')[:-1])
    if '.' in sn:
        sn = '-'.join(sn.split('.')[:-1])
    samples.append(sn)

dup_rate = next(l[1:] for l in report if l[0] == 'Dup rate')
reads_off_target = next(l[1:] for l in report if l[0] == '% off trg')
depth1x = next(l[1:] for l in report if l[0] == '1x')
depth10x = next(l[1:] for l in report if l[0] == '10x')
depth25x = next(l[1:] for l in report if l[0] == '25x')
depth100x = next(l[1:] for l in report if l[0] == '100x')

print '\t'.join(samples)
print '\t'.join([getval(v) for v in dup_rate])
print ''
print '\t'.join([getval(v) for v in reads_off_target])
print '\t'.join([getval(v) for v in depth1x])
print '\t'.join([getval(v) for v in depth10x])
print '\t'.join([getval(v) for v in depth25x])
print '\t'.join([getval(v) for v in depth100x])








