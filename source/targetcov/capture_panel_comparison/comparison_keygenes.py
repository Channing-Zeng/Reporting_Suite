with open('__comparison_keygenes_in.txt') as f:
    lines = f.readlines()

lns = []
with open('__comparison_keygenes_out.txt', 'w') as f:
    for l in lines:
        if l.startswith('##'):
            fl = 'Best values from all samples: '
            sns = l[48:].split(', ')
            for sn in sns:
                if '-' in sn:
                    sn = '-'.join(sn.split('-')[:-1])
                if '.' in sn:
                    sn = '-'.join(sn.split('.')[:-1])
                fl += sn + ', '
            fl += '\n'
            f.write(fl)
        else:
            c, s, e, sz, g, _, ft, bt, md, ad, sd, wn, x1, x5, x10, x25, x50, x100, x500, x1000, x5000, x10000, x50000 = l.strip().split('\t')
            if not l.startswith('#'):
                [wn, x1, x25, x100] = [str(float(v)) for v in [wn, x1, x25, x100]]
                lns.append([g, c, s, e, sz, bt, md, ad, sd, wn, x1, x25, x100])

    lns.sort(key=lambda l: l[0])
    for l in lns:
        f.write('\t'.join(l[6:]) + '\n')



#cat __comparison_tmp | cut -f9,10,11,12,13,16,18 > __comparison_tmp2

