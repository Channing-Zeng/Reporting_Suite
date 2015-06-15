#!/usr/bin/env python
import sys
import os
# from bs4 import BeautifulSoup
import collections
import itertools
from source import verify_file
from source.logger import err, info


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
    parsed_data = collections.OrderedDict((h, list()) for h in _header)

    for s in samples:
        if verify_file(s.fastqc_html_fpath):
            with open(s.fastqc_html_fpath) as source_file_obj:
                html = source_file_obj.read()
                parts = [p.split('</div>')[0] for p in html.split('<div class="module">')[1:]]
                # <h2><img/></h2><table></table></div>  OR  <h2><img/></h2><p><img/></p></div>
                for i, part in enumerate(parts):
                    # info('Parsing ' + _header[i])
                    info(str(part))
                    table, graph = None, None
                    ok_img = '<img ' + part.split('"><img ')[1].split('"/>')[0] + '"/>'
                    if '<table>' in part:
                        table = '<table>' + part.split('<table>')[1]
                    if '<p><img ' in part:
                        graph = '<img ' + part.split('<p><img ')[1].split('"/>')[0] + '"/>'
                    parsed_data[_header[i]].append([s.name, ok_img, graph, table])

                # module_divs = soup.find_all("div", class_="module")
                # _sort_graph_by_type(parsed_data, module_divs, s.name)
                # soup.decompose()
        else:
            err('Could not find fastqc html fpath for sample ' + s.name + ': ' + str(s.fastqc_html_fpath))

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
        print 'Graph for', bam_file_name
        print graph
        print ''
        print 'Table'
        print table
        print ''
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
