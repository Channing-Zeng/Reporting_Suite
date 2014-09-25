#!/usr/bin/env python
import sys
import os
from bs4 import BeautifulSoup
import collections
import itertools
import re

header = ["Basic Statistics",
          "Per base sequence quality",
          "Per tile sequence quality",
          "Per sequence quality scores",
          "Per base sequence content",
          "Per sequence GC content",
          "Per base N content",
          "Sequence Length Distribution",
          "Sequence Duplication Levels",
          "Overrepresented sequences",
          "Adapter Content",
          "Kmer Content"]


def get_graphs(input_files):
    m = collections.OrderedDict()
    for input_file in input_files:
        if os.path.exists(input_file):
            with open(input_file) as source_file_obj:
                broken_html = source_file_obj.read()
                # print broken_html
                soup = BeautifulSoup(broken_html)
                #divs_module = soup.select("div > div")
                bam_file_name = soup.title.string.strip(" FastQC Report")

                divs_module = soup.find_all("div", class_="module")

                sort_graph_by_type(m, divs_module, bam_file_name)
                soup.decompose()

    return m


def sort_graph_by_type(m, divs, bam_file_name):
    i = 0

    for div in divs:

        table = ""
        ok_img = ""
        graph = ""

        for content in div.contents:

            if content.name == 'table':
                table = content
            if content.name == 'h2':
                ok_img = content.next_element
            if content.name == 'p':
                graph = content.next_element

        if header[i] not in m:
            m[header[i]] = list()

        m[header[i]].append([bam_file_name, ok_img, graph, table])
        i += 1


def group2(iterator, count):
    if count > len(iterator): count = len(iterator)
    return itertools.imap(None, *([iter(iterator)] * count))


if __name__ == "__main__":
    m = get_graphs(sys.argv[1:])
    # print m
    #for keys, values in m.items():
    #print keys
    #for a, b in values:
    #print a
    #for a in group2(values, 3):
    #print  a[0][0]





