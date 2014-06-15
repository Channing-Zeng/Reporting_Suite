#!/usr/bin/env python
from __future__ import print_function
import fileinput
from scipy.stats import fisher_exact
from sys import argv, exit

def main():
    """Read vcf from first argument or stdin.
    Filter SomaticIndelDetector based InDels and pass everything else.
    """
    if len(argv) > 1 and (argv[1] == "-h" or argv[1] == "--help"):
        print(main.__doc__)
        exit(1)
    for line in fileinput.input():
        line = line.strip()
        if line.startswith("#") or not is_indel(line):
            print(line)
            # add a few header lines
            if "fileformat=VCF" in line:
                # ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
                print("##INFO=<ID=AFDIFF,Number=1,Type=Float,Description=\"Allele frequency difference between tumor and normal\">")
                print("##INFO=<ID=SIDFILTERS,Number=.,Type=String,Description=\"Custom filters for SomaticIndelDetector InDels {AFDIFF,SB,OFFSET}\">")
            continue
        lineSplit = line.split("\t")
        if len(lineSplit) < 10:
            # some weird line, print anyway
            print(line)
            continue
        myFilters = []
        add_af_diff_gt_filter(lineSplit, myFilters)
        strand_bias_filter(lineSplit, myFilters)
        offset_filter(lineSplit, myFilters)
        if lineSplit[6].strip() == ".":
            lineSplit[6] = "PASS"
        if len(myFilters) > 0:
            infoString = "SIDFILTERS=" + ",".join(myFilters)
            add_INFO(lineSplit, infoString)
        print("\t".join(lineSplit))
    
def is_indel(line):
    # split variant by tab
    myParts = line.split("\t")
    if len(myParts) < 5:
        return False
    # consider multiallelic sites
    myRefs = myParts[3].split(",")
    myAlts = myParts[4].split(",")
    # if any of the multiallelic variants are of length 2 or higher -> indel
    for myRef in myRefs:
        if len(myRef.strip()) > 1:
            return True
    for myAlt in myAlts:
        if len(myAlt.strip()) > 1:
            return True
    # must be a SNP
    return False        

def add_af_diff_gt_filter(linesplit, myFilters):
    """ Add allele frequency difference (as per maximum observer ALT) into INFO
    Further estimate Fisher's exact test based p-value for allele frequency
    difference (contingency table).
    """
    AD = get_format(linesplit, "AD")
    # a: sample1 allele count of the major allele 
    # b: sample1 count for the minor allele
    # c: sample2 allele count of the major allele
    # d: sample2 allele count for the minor allele
    if not AD == []:
        sampleAFs = []
        sampleCounts = []
        for sample in AD:
            maxNonRefAD = max(sample[1:])
            sampleAFs.append(float(maxNonRefAD)/(float(maxNonRefAD)+float(sample[0])))
            sampleCounts += [[sample[0], maxNonRefAD]]
        # test for significant difference in AF
        if len(sampleCounts) > 1:
            AFDiffSignificant = False
            for myCounts in sampleCounts[1:]:
                oddsratio, p_value = fisher_exact([sampleCounts[0],
                                               myCounts])
                if p_value < 0.05:
                    AFDiffSignificant = True
            if not AFDiffSignificant:
                myFilters.append("AFDIFF")
                add_FILTER(linesplit, "REJECT")
                    
        # add maximum allele frequency difference into INFO field
        if len(sampleAFs) == 1:
            sampleAFs += [float(0)]
        AFDiffString = "AFDIFF=" + str(max(sampleAFs)-min(sampleAFs))
        add_INFO(linesplit  , AFDiffString)

def strand_bias_filter(linesplit, myFilters):
    # SC: (FwdRef,RevRef,FwdIndel,RevIndel)
    # a: represents the forward strand allele count of the major allele 
    # b: represents the forward strand allele count for the minor allele
    # c: represents the reverse strand allele count of the major allele
    # d: represents the reverse strand allele count for the minor allele
    abcd = get_format(linesplit, "SC")
    genotypes = get_format(linesplit, "GT")
    SB = False
    if not abcd == []:
        for sample in abcd:
            a, c, b, d = int(sample[0]), int(sample[1]), int(sample[2]), int(sample[3])
            if b == 0 and d == 0 and ["0/0"] in genotypes: # not an indication of strand bias but no evidence for alternate allele
                continue 
            oddsratio, p_value = fisher_exact([[a,b],
                                               [c,d]])
            if p_value < 0.05 or b == 0 or d == 0:
                SB = True
    if SB == True:
        myFilters.append("SB")
        add_FILTER(linesplit, "REJECT")

def offset_filter(linesplit, myFilters):
    REnd = get_format(linesplit, "REnd")
    RStart = get_format(linesplit, "RStart")
    for offsetPair in REnd:
        if int(offsetPair[0]) <= 5 and int(offsetPair[1]) != 0:
            myFilters.append("OFFSET")
            add_FILTER(linesplit, "REJECT")
    for offsetPair in RStart:
        if int(offsetPair[0]) <= 5 and int(offsetPair[1]) != 0:
            myFilters.append("OFFSET")
            add_FILTER(linesplit, "REJECT")

def add_FILTER(linesplit, filterstring):
    if linesplit[6].strip() == ".":
        linesplit[6] = filterstring
    else:
        if filterstring not in linesplit[6]:
            linesplit[6] += ";" + filterstring

def add_INFO(linesplit, infostring):
    linesplit[7] += ";" + infostring
    linesplit[7] = ";".join(sorted(linesplit[7].split(";")))
    
def get_format(linesplit, key):
    """ For given FORMAT key (e.g. "AD") return list of values for all samples.
    Returned as a list of lists (samples by values)
    """
    myKeys = linesplit[8].split(":")
    #import pdb; pdb.set_trace()
    if key in myKeys:
        myIndex = myKeys.index(key)
        output = []
        for sample in linesplit[9:]:
            samplesplit = sample.split(":")
            if len(samplesplit) == len(myKeys) and sample != ".":
                output.append(samplesplit[myIndex].split(","))
        return output
    else:
        return None

if __name__ == "__main__":
    main()

"""
##FORMAT=<ID=MM,Number=2,Type=Float,Description="Average # of mismatches per consensus indel-supporting read/per reference-supporting read">
##FORMAT=<ID=MQS,Number=2,Type=Float,Description="Average mapping qualities of consensus indel-supporting reads/reference-supporting reads">
##FORMAT=<ID=NQSBQ,Number=2,Type=Float,Description="Within NQS window: average quality of bases from consensus indel-supporting reads/from reference-supporting reads">
##FORMAT=<ID=NQSMM,Number=2,Type=Float,Description="Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
##FORMAT=<ID=REnd,Number=2,Type=Integer,Description="Median/mad of indel offsets from the ends of the reads">
##FORMAT=<ID=RStart,Number=2,Type=Integer,Description="Median/mad of indel offsets from the starts of the reads">
##FORMAT=<ID=SC,Number=4,Type=Integer,Description="Strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  normal     tumor
chr1    9034815 .       C       CAGA    0.0     .       
AC=2;AF=0.500;AN=4;EFF=DOWNSTREAM(MODIFIER||313||313|CA6||CODING|NM_001270500.1||1),UTR_3_PRIME(MODIFIER||||180|CA6||CODING|NM_001270502.1|6|1),UTR_3_PRIME(MODIFIER||||248|CA6||CODING|NM_001270501.1|7|1),UTR_3_PRIME(MODIFIER||||308|CA6||CODING|NM_001215.3|8|1)    
GT:AD:DP:MM:MQS:NQSBQ:NQSMM:REnd:RStart:SC      
0/1:20,112:133:0.3,0.16071428:60.0,60.0:31.158228,34.844643:0.0,0.0:41,18:56,18:5,15,35,77      
0/1:12,69:81:0.0,0.1884058:60.0,60.0:33.919193,35.169567:0.0,0.0:39,15:58,15:7,5,22,47
GT     0/1
AD     12,69
DP     81
MM     0.0,0.1884058
MQS    60.0,60.0
NQSBQ  33.919193,35.169567
NQSMM  0.0,0.0
REnd   39,15
RStart 58,15
SC     7,5,22,47 


"""
