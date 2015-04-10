find . -mindepth 2 -maxdepth 4 -follow -path '*_targetSeq/*' -not -name 'picard_*' -type f -delete
find . -mindepth 2 -follow -path 'work/*' not -name 'picard_*' -type f -delete
find . -mindepth 2 -maxdepth 3 -follow -name '*.targetSeq*' -type f -delete
find . -mindepth 2 -maxdepth 3 -follow -name '*.targetSeq*' -type l -delete
find . -mindepth 2 -maxdepth 3 -follow -name 'targQC.*' -type f -delete

find . -maxdepth 2 -type f -path '*work/*' -delete


#####################################
# find . -mindepth 2 -maxdepth 3 -follow -name work_* -type d -exec rm -r {} \;
# find . -mindepth 2 -maxdepth 3 -follow -name work* -type d -exec rm -r {} \;
# find . -mindepth 2 -maxdepth 3 -follow -name *_targetSeq -type d -exec rm -r {} \;

# find . -mindepth 2 -maxdepth 3 -follow -name *targetSeq.details.gene.txt -type l -exec rm -r {} \;
# find . -mindepth 2 -maxdepth 3 -follow -name targQC.* -type f -exec rm -r {} \;
# find . -mindepth 2 -maxdepth 3 -follow -name Best.targetSeq.details.gene.tsv -type f -exec rm -r {} \;
