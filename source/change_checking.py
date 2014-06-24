import hashlib
from os.path import basename

from source.logger import info
from source.utils_from_bcbio import file_exists, open_gzipsafe


def md5_for_file(f, block_size=2**20):
    md5 = hashlib.md5()

    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)

    return md5.hexdigest()


def check_file_changed(cnf, new, in_work):
    if not file_exists(in_work):
        cnf['reuse_intermediate'] = False

    if cnf.get('reuse_intermediate'):
        if (basename(in_work) != basename(new) or
            md5_for_file(open(in_work, 'rb')) !=
            md5_for_file(open_gzipsafe(new, 'rb'))):

            info('Input file %s changed, setting "reuse_intermediate" '
                'to False.' % str(new))
            cnf['reuse_intermediate'] = False


# def check_inputs_changed(cnf, new_inputs):
#     prev_input_fpath = join(cnf['work_dir'], 'prev_inputs.txt')
#
#     new_inp_hashes = {realpath(fn): md5_for_file(fn) for fn in new_inputs}
#
#     if cnf.get('reuse_intermediate'):
#         if not file_exists(prev_input_fpath):
#             info('File %s does not exist, setting "reuse_intermediate" to '
#                       'False.' % str(prev_input_fpath))
#             cnf['reuse_intermediate'] = False
#
#         else:
#             prev_inp_hashes = dict()
#
#             with open(prev_input_fpath) as f:
#                 for l in f:
#                     input_fname, md5 = l.strip().split('\t')
#                     prev_inp_hashes[input_fname] = md5
#
#             if len(new_inp_hashes) != len(prev_inp_hashes):
#                 info('Number of input files changed, setting "reuse_intermediate" to False.')
#                 cnf['reuse_intermediate'] = False
#
#             for inp_fpath, inp_hash in new_inp_hashes.items():
#                 if inp_fpath not in prev_inp_hashes:
#                     info('Input changed, setting "reuse_intermediate" to False.')
#                     cnf['reuse_intermediate'] = False
#
#                 if inp_hash != prev_inp_hashes[inp_fpath]:
#                     info('Input %s changed, setting "reuse_intermediate" '
#                               'to False.' % str(inp_fpath))
#                     cnf['reuse_intermediate'] = False
#
#     with open(prev_input_fpath, 'w') as f:
#         for inp_fpath, inp_hash in new_inp_hashes.items():
#             f.write(inp_fpath + '\t' + inp_hash + '\n')
