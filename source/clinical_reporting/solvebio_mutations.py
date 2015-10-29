from source.file_utils import verify_module
from source.logger import err


def search_mutations(mutations):
    if not verify_module('solvebio'):
        err('Cannot import solvebio')
        return None

    from solvebio import Depository

    Depository.all()