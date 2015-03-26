with open('__comparison_oncomine_in.txt') as f:
	lines = f.readlines()

with open('__comparison_oncomine_out.txt', 'w') as f:
	i = 2
	header = []
	for l in lines[i:]:
		if not l.startswith('##'):
			break
		i += 1
		##PBMC-Resistant-4-ready ave depth=8.55288701409
		sn = l.split()[0][2:]
		if '-' in sn:
			sn = '-'.join(sn.split('-')[:-1])	
		if '.' in sn:
			sn = '-'.join(sn.split('.')[:-1])
		dp = '{:.2f}'.format(float(l.split()[-1].split('=')[-1]))
		print sn
		print dp
		header.append([sn, dp])

	ts = lines[i][1:].split('\t')
	f.write('\t'.join(ts[0:3] + ts[4:9]))
	for sn, dp in header:
		f.write('\t' + sn + ', ' + dp + ' ave dp')
	f.write('\n')

	for l in lines[i+1:]:
		ts = l.split('\t')
		f.write('\t'.join(ts[0:3] + ts[4:]))
