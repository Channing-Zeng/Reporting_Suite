# with open('__comparison_keygenes_in.txt') as f:
#     lines = f.readlines()
#
# with open('__comparison_keygenes_out.txt', 'w') as f:
#     for l in lines:
#         if l.startswith('##'):
#             f.write('Best values from all samples: ')
#             sns = l[48:].split(', ')
#             for sn in sns:
#                 if '-' in sn:
#                     sn = '-'.join(sn.split('-')[:-1])
#                 if '.' in sn:
#                     sn = '-'.join(sn.split('.')[:-1])
#                 f.write(sn + ', ')
#             f.write('\n')
#         else:
#             c, s, e, sz, g, _, ft, bt, md, ad, sd, wn, x1, x5, x10, x25, x50, x100, x500, x1000, x5000, x10000, x50000 = l.strip().split('\t')
#             if not l.startswith('#'):
#                 [wn, x1, x25, x100] = [str(float(v) * 100) for v in [wn, x1, x25, x100]]
#             f.write('\t'.join([md, ad, sd, wn, x1, x25, x100]) + '\n')


#cat __comparison_tmp | cut -f9,10,11,12,13,16,18 > __comparison_tmp2



def key_to_sort(name):
    parts = []

    cur_part = []
    prev_was_num = False

    for c in name:
        if prev_was_num == c.isdigit():
            cur_part.append(c)
        else:
            part = ''.join(cur_part)
            if prev_was_num:
                part = int(part)
            parts.append(part)
            cur_part = []

    return parts


names = ['DNA_10-ready', 'DNA_1-ready', 'DNA_2-ready']
names.sort(key=lambda n: key_to_sort(n))
print names