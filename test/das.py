contigs = '[chrM, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY, chr1_gl000191_random, chr1_gl000192_random, chr4_ctg9_hap1, chr4_gl000193_random, chr4_gl000194_random, chr6_apd_hap1, chr6_cox_hap2, chr6_dbb_hap3, chr6_mann_hap4, chr6_mcf_hap5, chr6_qbl_hap6, chr6_ssto_hap7, chr7_gl000195_random, chr8_gl000196_random, chr8_gl000197_random, chr9_gl000198_random, chr9_gl000199_random, chr9_gl000200_random, chr9_gl000201_random, chr11_gl000202_random, chr17_ctg5_hap1, chr17_gl000203_random, chr17_gl000204_random, chr17_gl000205_random, chr17_gl000206_random, chr18_gl000207_random, chr19_gl000208_random, chr19_gl000209_random, chr21_gl000210_random, chrUn_gl000211, chrUn_gl000212, chrUn_gl000213, chrUn_gl000214, chrUn_gl000215, chrUn_gl000216, chrUn_gl000217, chrUn_gl000218, chrUn_gl000219, chrUn_gl000220, chrUn_gl000221, chrUn_gl000222, chrUn_gl000223, chrUn_gl000224, chrUn_gl000225, chrUn_gl000226, chrUn_gl000227, chrUn_gl000228, chrUn_gl000229, chrUn_gl000230, chrUn_gl000231, chrUn_gl000232, chrUn_gl000233, chrUn_gl000234, chrUn_gl000235, chrUn_gl000236, chrUn_gl000237, chrUn_gl000238, chrUn_gl000239, chrUn_gl000240, chrUn_gl000241, chrUn_gl000242, chrUn_gl000243, chrUn_gl000244, chrUn_gl000245, chrUn_gl000246, chrUn_gl000247, chrUn_gl000248, chrUn_gl000249]'[1:-1].split(', ')

asdad1 = contigs
for c in contigs:
    print c,

print

c2 = '''##contig=<ID=chrM,length=16571,assembly=hg19>
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chr1_gl000191_random,length=106433,assembly=hg19>
##contig=<ID=chr1_gl000192_random,length=547496,assembly=hg19>
##contig=<ID=chr4_ctg9_hap1,length=590426,assembly=hg19>
##contig=<ID=chr4_gl000193_random,length=189789,assembly=hg19>
##contig=<ID=chr4_gl000194_random,length=191469,assembly=hg19>
##contig=<ID=chr6_apd_hap1,length=4622290,assembly=hg19>
##contig=<ID=chr6_cox_hap2,length=4795371,assembly=hg19>
##contig=<ID=chr6_dbb_hap3,length=4610396,assembly=hg19>
##contig=<ID=chr6_mann_hap4,length=4683263,assembly=hg19>
##contig=<ID=chr6_mcf_hap5,length=4833398,assembly=hg19>
##contig=<ID=chr6_qbl_hap6,length=4611984,assembly=hg19>
##contig=<ID=chr6_ssto_hap7,length=4928567,assembly=hg19>
##contig=<ID=chr7_gl000195_random,length=182896,assembly=hg19>
##contig=<ID=chr8_gl000196_random,length=38914,assembly=hg19>
##contig=<ID=chr8_gl000197_random,length=37175,assembly=hg19>
##contig=<ID=chr9_gl000198_random,length=90085,assembly=hg19>
##contig=<ID=chr9_gl000199_random,length=169874,assembly=hg19>
##contig=<ID=chr9_gl000200_random,length=187035,assembly=hg19>
##contig=<ID=chr9_gl000201_random,length=36148,assembly=hg19>
##contig=<ID=chr11_gl000202_random,length=40103,assembly=hg19>
##contig=<ID=chr17_ctg5_hap1,length=1680828,assembly=hg19>
##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>
##contig=<ID=chr17_gl000204_random,length=81310,assembly=hg19>
##contig=<ID=chr17_gl000205_random,length=174588,assembly=hg19>
##contig=<ID=chr17_gl000206_random,length=41001,assembly=hg19>
##contig=<ID=chr18_gl000207_random,length=4262,assembly=hg19>
##contig=<ID=chr19_gl000208_random,length=92689,assembly=hg19>
##contig=<ID=chr19_gl000209_random,length=159169,assembly=hg19>
##contig=<ID=chr21_gl000210_random,length=27682,assembly=hg19>
##contig=<ID=chrUn_gl000211,length=166566,assembly=hg19>
##contig=<ID=chrUn_gl000212,length=186858,assembly=hg19>
##contig=<ID=chrUn_gl000213,length=164239,assembly=hg19>
##contig=<ID=chrUn_gl000214,length=137718,assembly=hg19>
##contig=<ID=chrUn_gl000215,length=172545,assembly=hg19>
##contig=<ID=chrUn_gl000216,length=172294,assembly=hg19>
##contig=<ID=chrUn_gl000217,length=172149,assembly=hg19>
##contig=<ID=chrUn_gl000218,length=161147,assembly=hg19>
##contig=<ID=chrUn_gl000219,length=179198,assembly=hg19>
##contig=<ID=chrUn_gl000220,length=161802,assembly=hg19>
##contig=<ID=chrUn_gl000221,length=155397,assembly=hg19>
##contig=<ID=chrUn_gl000222,length=186861,assembly=hg19>
##contig=<ID=chrUn_gl000223,length=180455,assembly=hg19>
##contig=<ID=chrUn_gl000224,length=179693,assembly=hg19>
##contig=<ID=chrUn_gl000225,length=211173,assembly=hg19>
##contig=<ID=chrUn_gl000226,length=15008,assembly=hg19>
##contig=<ID=chrUn_gl000227,length=128374,assembly=hg19>
##contig=<ID=chrUn_gl000228,length=129120,assembly=hg19>
##contig=<ID=chrUn_gl000229,length=19913,assembly=hg19>
##contig=<ID=chrUn_gl000230,length=43691,assembly=hg19>
##contig=<ID=chrUn_gl000231,length=27386,assembly=hg19>
##contig=<ID=chrUn_gl000232,length=40652,assembly=hg19>
##contig=<ID=chrUn_gl000233,length=45941,assembly=hg19>
##contig=<ID=chrUn_gl000234,length=40531,assembly=hg19>
##contig=<ID=chrUn_gl000235,length=34474,assembly=hg19>
##contig=<ID=chrUn_gl000236,length=41934,assembly=hg19>
##contig=<ID=chrUn_gl000237,length=45867,assembly=hg19>
##contig=<ID=chrUn_gl000238,length=39939,assembly=hg19>
##contig=<ID=chrUn_gl000239,length=33824,assembly=hg19>
##contig=<ID=chrUn_gl000240,length=41933,assembly=hg19>
##contig=<ID=chrUn_gl000241,length=42152,assembly=hg19>
##contig=<ID=chrUn_gl000242,length=43523,assembly=hg19>
##contig=<ID=chrUn_gl000243,length=43341,assembly=hg19>
##contig=<ID=chrUn_gl000244,length=39929,assembly=hg19>
##contig=<ID=chrUn_gl000245,length=36651,assembly=hg19>
##contig=<ID=chrUn_gl000246,length=38154,assembly=hg19>
##contig=<ID=chrUn_gl000247,length=36422,assembly=hg19>
##contig=<ID=chrUn_gl000248,length=39786,assembly=hg19>
##contig=<ID=chrUn_gl000249,length=38502,assembly=hg19>'''.split()

asdad2=[]
for c in c2:
    c_n = c.split(',')[0][13:]
    asdad2.append(c_n)

for c1, c2 in zip(asdad1, asdad2):
    print c1, c2, c1==c2

