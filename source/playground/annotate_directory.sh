#! /bin/bash
# Warning: globstar excludes hidden directories.
# Turn on recursive globbing (in this script) or exit if the option is not supported:
shopt -s globstar || exit

source /etc/profile.d/modules.sh
module load bcbio-nextgen/0.7.6

for f in **
do
  if [[ "$f" =~ \.txt$ ]] ; then
    echo "$f"
    qsub -j y -cwd -q batch.q -l h=chara -l huge_ram=1 /gpfs/users/klpf990/az_scripts/test_annots/filter_and_annotate.sh $f
  fi
done