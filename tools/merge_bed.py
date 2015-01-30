def _merge_bed(cnf, bed_fpath, collapse_gene_names=True):
    output_fpath = intermediate_fname(cnf, bed_fpath, 'merge')

    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = ('{bedtools} merge ' + ('-c 4,5,6 -o distinct' if collapse_gene_names
               else '') + ' -i {bed_fpath}').format(**locals())
    res = call(cnf, cmdline, output_fpath, exit_on_error=False)

    if res is not None:
        return res
    else:
        warn('Old version of bedtools. Trying different arguments...')

        cmdline = ('{bedtools} merge ' + ('-nms -scores collapse' if collapse_gene_names
                   else '') + ' -i {bed_fpath}').format(**locals())
        call(cnf, cmdline, output_fpath)

        def fn(l, i):
            ts = l.split('\t')
            if len(ts) < 4:
                return l
            gns = ts[3].split(';')
            l = l.replace(ts[3], ','.join(set(gns)))
            return l

        return iterate_file(cnf, output_fpath, fn, 'dstnct')
