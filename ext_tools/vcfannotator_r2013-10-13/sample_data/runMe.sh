#!/bin/bash -e

../VCF_annotator.pl --gff3 Eco_TY_2482.gff3 --genome Eco_TY_2482.genome --vcf 04-8351_TY-2482.vcf.snps_only  > 04-8351_TY-2482.vcf.snps_only.annotated

echo Done. See file: 04-8351_TY-2482.vcf.snps_only.annotated

