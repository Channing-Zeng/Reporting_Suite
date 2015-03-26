import sys
import os
import shutil


def __remove_quotes(str):
    if str and str[0] == '"':
        str = str[1:]
    if str and str[-1] == '"':
        str = str[:-1]
    return str


def __correct_cosmic_dbsnp(sample_fpath):
    result_fpath = sample_fpath + '_tmp'

    with open(sample_fpath) as vcf, open(result_fpath, 'w') as out:
        for i, line in enumerate(vcf):
            clean_line = line.strip()
            if not clean_line or clean_line[0] == '#':
                out.write(line)
            else:
                tokens = line.split()
                chr_field = tokens[0]
                if not chr_field.startswith('chr'):
                    chr_field = 'chr' + chr_field
                alls = __remove_quotes(tokens[4])
                line = '\t'.join([chr_field] + [tokens[1]] + ['.'] + tokens[3:4] + [alls] + tokens[5:]) + '\n'
                out.write(line)

    os.remove(sample_fpath)
    shutil.move(result_fpath, sample_fpath)
    return sample_fpath


if __name__ == '__main__':
    __correct_cosmic_dbsnp(sys.argv[1])