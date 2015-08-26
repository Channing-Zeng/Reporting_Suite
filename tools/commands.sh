find . -name varqc -exec bash -c 'mv "{}" "$(dirname "{}")"/qc/varQC' \;
find . -name varqc_after -exec bash -c 'mv "{}" "$(dirname "{}")"/qc/varQC_postVarFilter' \;
find . -name varannotate -exec bash -c 'mv "{}" "$(dirname "{}")"/varAnnotate' \;
find . -name varfilter -exec bash -c 'mv "{}" "$(dirname "{}")"/varFilter' \;
find . -name ngscat -exec bash -c 'mv "{}" "$(dirname "{}")"/qc/ngscat' \;
find . -name qualimap -exec bash -c 'mv "{}" "$(dirname "{}")"/qc/qualimap' \;
find . -name targetcov -exec bash -c 'mv "{}" "$(dirname "{}")"/targetSeq' \;
for DIR in ./*; do mkdir $DIR/var; done


~/symlink_filt_vcf.sh .

find . -wholename "*/*-vardict.vcf*" -exec bash -c 'mv {} $(dirname {})/var' \;
find . -wholename "*/*-freebayes.vcf*" -exec bash -c 'mv {} $(dirname {})/var' \;
find . -wholename "*/*-mutect.vcf*" -exec bash -c 'mv {} $(dirname {})/var' \;


find . -wholename "*/varFilter/*.filt.vcf" -exec bash -c 'ln -s {} $(dirname {})/../' \;
find . -wholename "*/varFilter/*.filt.maf" -exec bash -c 'ln -s {} $(dirname {})/../' \;
find . -wholename "*/varFilter/*.filt.tsv" -exec bash -c 'ln -s {} $(dirname {})/../' \;


find . -name "*mble.bed" -exec bash -c 'mkdir $(dirname {})/cnv; mv {} $(dirname {})/cnv' \;
find . -name "*mops.bed" -exec bash -c 'mkdir $(dirname {})/cnv; mv {} $(dirname {})/cnv' \;

rsync -tavz --include '*/' --include '*mble.bed' --exclude '*' klpf990@ukapdlnx115.ukapd.astrazeneca.net:/ngs/oncology/Analysis/translation/TS_0025_Platinum_Resistant_Ovarian/NGS-99/final/* .

rsync -tavz --include '*/' --include '*mble.bed' --exclude '*' klpf990@ukapdlnx115.ukapd.astrazeneca.net:/ngs/oncology/analysis/translation/TS_0023_TGEN_DLBCL/exomes/final/* .
rsync -tavz --include '*/' --include '*mops.bed' --exclude '*' --prune-empty-dirs klpf990@ukapdlnx115.ukapd.astrazeneca.net:/ngs/oncology/analysis/translation/TS_0023_TGEN_DLBCL/exomes/final/* .

find . -type d -empty -delete


rsync -tavz --include '*/' --include 'varQC' --exclude '*' --prune-empty-dirs . ~/NGS-99-reports



rsync -tavz --include '*/' --include '*varFilter/*vardict.anno.filt.*' --exclude '*' --prune-empty-dirs . /gpfs/ngs/oncology/Analysis/translation/TS_0025_Platinum-Resistant-Ovarian/bcbio/NGS-99/final


klpf990@ukapdlnx115.ukapd.astrazeneca.net:/ngs/oncology/Analysis/translation/TS_0025_Platinum_Resistant_Ovarian/NGS-99/final/* .
