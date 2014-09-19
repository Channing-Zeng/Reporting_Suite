#!/usr/bin/env python
import sys
import os
from bs4 import BeautifulSoup
import collections

header = ["Basic Statistics",
          "Per base sequence quality",
          "Per tile sequence quality",
          "Per sequence quality scores",
          "Per base sequence content",
          "Per sequence GC content",
          "Per base N content",
          "Sequence Length Distribution",
          "Sequence Duplication Levels",
          "Overrepresented sequences"
          "Adapter Content",
          "Kmer Content"]


def get_graphs(input_files):
    m = collections.OrderedDict()
    for input_file in input_files:
        if os.path.exists(input_file):
            with open(input_file) as source_file_obj:
                broken_html = source_file_obj.read()
                soup = BeautifulSoup(broken_html)
                bam_file_name = soup.title.string.strip(" FastQC Report")
                #report_name = input_file
                divs_with_class = (div for div in soup.find_all('div') if
                                   div.get('class') is not None and 'module' in div.get('class'))

                sort_graph_by_type(m, divs_with_class, bam_file_name)
    return m


def sort_graph_by_type(m, divs, bam_file_name):
    i = 0
    for div in divs:
        if header[i] not in m:
            m[header[i]] = collections.OrderedDict()

        m[header[i]][bam_file_name] = div
        i += 1


if __name__ == "__main__":
    m = get_graphs(sys.argv[1:])
    #print m
    for keys, values in m.items():
        print keys
        for a, b in values.items():
            print a

