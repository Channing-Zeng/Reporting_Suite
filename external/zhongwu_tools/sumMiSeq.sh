#!/bin/bash

mapping_stat_summary.pl -s read_counts.txt */stat* > mapping_summary_human.txt
cd PhiX
mapping_stat_summary.pl -s ../read_counts.txt */stat* > ../mapping_summary_PhiX.txt
cd ..
joinLine mapping_summary_PhiX.txt mapping_summary_human.txt  | getcols -c 1..2:5:8:11:10:12 | perl -ne '{ if ( /^Sample/ ) { print join("\t", qw(Sample Mapped_human %Mapped_human Mapped_PhiX %Mapped_PhiX Total %Total)), "\n"; } else { print; } }' > mapping_summary.txt
mapping_cov_summary.pl -d 1:5:10:25:50:100:500:1000:5000:10000:50000 */*cov* > ROISample.txt
cat */*cov.txt | grep -v Whole-G | grep -v Sample | pivot.pl -I "Gene\tChr\tStart\tEnd" -i 2..5 -c 1 -v 8 > DepthTargets.txt
cat */*cov.txt | grep -v Whole | grep -v Sample | myGroupByFunc -g "2:3..5" -c 8 -f mean:min:max -I "Gene:Chr:Start:End" > mean_amplicon_cov.txt

cat */*cosmic.vcf | grep -v raw_ | vcf2txt.pl -u -p 5 -q 25 -n 6 -F 0.01 | pickLine -c 35 -l PASS:TRUE > variants_PASS.txt

cat variants_PASS.txt | perl -ne '{if ( /^Sample/ ) { print; next; } @a = split(/\t/); next if ( $a[20] < 8 || $a[21] < 0.025 || $a[20]*$a[21] < 4.95); next if ( ($a[22] eq "2;1" || $a[22] eq "2;0") && $a[35] =~ /Novel/ && $a[21] < 0.3); next unless( $a[24] > 0 ); print;}' > variants_PASS.filtered.txt
grep -vw dbSNP variants_PASS.filtered.txt | grep -v UTR_ | grep -vw SILENT | grep -v INTRON | grep -v UPSTREAM | grep -v DOWNSTREAM | grep -v INTERGENIC | grep -v INTRAGENIC | grep -v NON_CODING | pickLine -v -i 12:3 -c 13:11 /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/snpeffect_export_Polymorphic.txt > variants_PASS.filtered.mut.txt

