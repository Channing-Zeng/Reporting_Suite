AZ pipeline for postprocessing BCBio (bcbio-nextgen.readthedocs.org) WGS and WES analysis results.

Includes:
- variant annotation,
- variant quality control,
- variant filtration,
- amplicon and exon coverage statistics using built-in means and also ngsCAT (http://www.bioinfomgp.org/ngscat),
- integrated Qualimap for alignment quality control (http://qualimap.bioinfo.cipf.es),
- CNV caller Seq2C.

Main running script: 
az-reporting.py [path to BCBio final directory] [--run-cnf run_info.yaml] [--sys-cnf system_info.yaml]
