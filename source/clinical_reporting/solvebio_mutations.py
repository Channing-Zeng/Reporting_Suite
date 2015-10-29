from os.path import join, expanduser
from source.file_utils import verify_module, verify_file
from source.logger import err, info, warn


class SolvebioRecord:
    def __init__(self, clinsig=None, url=None):
        self.clinsig = clinsig
        self.url = url


def parse_response(res, mut):
    ok = True
    for f in ['allele_origin',
              'clinical_significance',
              'genomic_coordinates']:
        if f not in res:
            warn('No ' + f + ' in SolveBio for mutation ' + str(mut))
            ok = False
    if not ok:
        return None

    rec = SolvebioRecord()

    rec.clinsig = res['clinical_significance']
    if rec.clinsig.lower() == 'other':
        rec.clinsig = 'Uncertain'

    coords = res['genomic_coordinates']
    rec.url = 'https://www.solvebio.com/variant/GRCH37-{chrom}-{start}-{stop}-{alt}'.format(
        chrom=coords['chromosome'],
        start=coords['start'],
        stop=coords['stop'],
        alt=res['allele'])

    return rec

def query_mutations(cnf, mutations):
    info('')
    info('SolveBio')
    if not verify_module('solvebio'):
        err('Cannot import solvebio')
        return None

    saves_fpath = join(cnf.work_dir, 'solvebio.txt')
    if cnf.debug:
        if verify_file(saves_fpath, silent=True):
            info('Debug, reading hits from ' + saves_fpath)
            with open(saves_fpath) as f:
                for mut, l in zip(mutations, f):
                    if l:
                        l = l.strip().split('\t')
                        mut.solvebio = SolvebioRecord(clinsig=l[0], url=l[1])
            info('Done, read ' + str(sum(1 for m in mutations if m.solvebio)) + ' hits')
            return mutations

    from solvebio import login, Depository, Dataset, BatchQuery
    login()
    ds = Dataset.retrieve('ClinVar/Variants')

    queries = [
        ds.query().filter(gene_symbol=m.gene.name) \
            .position(m.chrom, position=m.pos) \
            .filter(allele=m.alt, reference_allele=m.ref)
        for m in mutations]

    info('Querying mutations against SolveBio ClinVar depository')
    for mut, res in zip(mutations, BatchQuery(queries).execute()):
        if res['results']:
            solvbio_record = parse_response(res['results'][0], mut)
            if solvbio_record:
                mut.solvebio = solvbio_record
    info('Done, found ' + str(sum(1 for m in mutations if m.solvebio)) + ' hits')

    if cnf.debug:
        with open(saves_fpath, 'w') as f:
            for mut in mutations:
                if mut.solvebio:
                    f.write(mut.solvebio.clinsig)
                    f.write('\t')
                    f.write(mut.solvebio.url)
                    f.write('\n')
        info('Saved to ' + saves_fpath)

    info('-' * 70)
    return mutations

