# -*- coding: windows-1252 -*-
from source.ngscat.xlwt import Workbook, Column, Worksheet, Row

__VERSION__ = '0.7.2'

import sys
if sys.version_info[:2] < (2, 3):
    print >> sys.stderr, "Sorry, xlwt requires Python 2.3 or later"
    sys.exit(1)

from source.ngscat.xlwt.Formatting import Font, Alignment, Borders, Pattern, Protection
from source.ngscat.xlwt.Style import XFStyle, easyxf
from source.ngscat.xlwt.ExcelFormula import *
