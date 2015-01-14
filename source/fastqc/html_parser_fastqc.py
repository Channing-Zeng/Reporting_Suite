#!/usr/bin/env python
import sys
import os
from bs4 import BeautifulSoup
import collections
import itertools
from source import verify_file


_header = ["Basic Statistics",
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


def get_graphs(samples):
    parsed_data = collections.OrderedDict()
    for s in samples:
        if verify_file(s.fastqc_html_fpath):
            with open(s.fastqc_html_fpath) as source_file_obj:
                html = source_file_obj.read()
                soup = BeautifulSoup(html)
                #bam_file_name = soup.title.string.strip(" FastQC Report")

                module_divs = soup.find_all("div", class_="module")
                _sort_graph_by_type(parsed_data, module_divs, s.name)
                soup.decompose()

    return parsed_data


def _sort_graph_by_type(parsed_data, divs, bam_file_name):
    i = 0
    for div in divs:
        table = ""
        ok_img = ""
        graph = ""
        for content in div.contents:
            if 'table' == content.name: table = content
            if 'h2' == content.name: ok_img = content.next_element
            if 'p' == content.name: graph = content.next_element

        if _header[i] not in parsed_data:
            parsed_data[_header[i]] = list()
        parsed_data[_header[i]].append([bam_file_name, ok_img, graph, table])
        i += 1


def _group2(iterator, count):
    if count > len(iterator): count = len(iterator)
    return itertools.imap(None, *([iter(iterator)] * count))


if __name__ == "__main__":
    k = 0
    d = collections.OrderedDict()
    for i in sys.argv[1:]:
        k += 1
        d[k]=i
    m = get_graphs(d)

    #print m
    #for keys, values in m.items():
        #print keys
        #for a, b in values:
            #print a
    #for a in group2(values, 3):
    #print  a[0][0]





