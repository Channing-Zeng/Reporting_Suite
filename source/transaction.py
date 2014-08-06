"""Handle file based transactions allowing safe restarts at any point.

To handle interrupts,this defines output files written to temporary
locations during processing and copied to the final location when finished.
This ensures output files will be complete independent of method of
interruption.
"""
import os
import shutil
import contextlib
from os.path import exists, join
from source.file_utils import add_suffix


