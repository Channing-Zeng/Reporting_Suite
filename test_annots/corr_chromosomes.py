import sys
import os


def __correct_cromosomes(sample_fpath):
    result_fpath = sample_fpath + '_tmp'

    with open(sample_fpath) as vcf, open(result_fpath, 'w') as out:
        for i, line in enumerate(vcf):
            clean_line = line.strip()
            if not clean_line or clean_line[0] == '#':
                out.write(line)
            else:
                tokens = line.split()
                chr_field = tokens[0]
                line = '\t'.join(['chr' + chr_field] + tokens[1:]) + '\n'
                out.write(line)

    os.remove(sample_fpath)
    os.rename(result_fpath, sample_fpath)
    return sample_fpath


if __name__ == '__main__':
    __correct_cromosomes(sys.argv[1])