grep -v "\|^H\|^G" 
grep -v "^H"


Ensembl.gtf
1       miRNA   gene    5875931 5876019 .       -       .       gene_id "ENSG00000216045"; gene_name "AL356693.1"; gene_source "ensembl"; gene_biotype "miRNA";
1       miRNA   transcript      5875931 5876019 .       -       .       gene_id "ENSG00000216045"; transcript_id "ENST00000401226"; gene_name "AL356693.1"; gene_source "ensembl"
1       miRNA   exon    5875931 5876019 .       -       .       gene_id "ENSG00000216045"; transcript_id "ENST00000401226"; exon_number "1"; gene_name "AL356693.1"; gene_source

subtruction
chr1    5875992 5876019 AL356693.1      chr1    5875931 5876019 miRNA   gene    27
chr1    5875992 5876019 AL356693.1      chr1    5875931 5876019 miRNA   transcript      27
chr1    5875992 5876019 AL356693.1      chr1    5875931 5876019 miRNA   exon    27


# get the union.bed
grep -v '^#' /ngs/reference_data/genomes/Hsapiens/hg19/bed.20150528/Exons/Ensembl.gtf | grep -w 'transcript' > Ensembl.transcripts.gtf
grep 'protein_coding\|nonsense_mediated_decay\|miRNA\|IG_C_gene\|IG_D_gene\|IG_J_gene\|IG_V_gene\|TR_C_gene\|TR_D_gene\|TR_J_gene\|TR_V_gene' Ensembl.transcrits.gtf > Ensembl.transcrits.protein_coding.gtf
cut -f1,4,5 Ensembl.transcrits.protein_coding.gtf > transcripts.txt
bedtools slop -l 1 -r 0 -i transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > transcripts.bed
grch37_to_hg19.py transcripts.bed | grep -v "^chrH" | grep -v "^chrG" > transcripts.fixedchr.bed
cd ..
cat padded_panels/*.bed transcripts.bed > joined.bed
sort -k1,1 -k2,2n joined.bed > joined.sorted.bed
bedtools merge -i joined.sorted.bed > union.bed

# subtruct
bedtools subtract -b union.bed -a ../AZ-Exome.bed > subtraction.bed

# get the bed for annotation
awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed

# intersect
bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
cut -f8 subtraction.annotated.bed | sort -u

bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
cut -f7 subtraction.back.annotated.bed | sort -u





# get the union.exons.bed
grep -v '^#' /ngs/reference_data/genomes/Hsapiens/hg19/bed.20150528/Exons/Ensembl.gtf | grep -w 'CDS\|UTR\|stop_codon' | grep 'protein_coding\|nonsense_mediated_decay\|miRNA\|IG_C_gene\|IG_D_gene\|IG_J_gene\|IG_V_gene\|TR_C_gene\|TR_D_gene\|TR_J_gene\|TR_V_gene' > Ensembl.cds_utr.gtf
grep -w 'CDS' Ensembl.cds_utr.gtf | cut -f1,4,5 > Ensembl.csd.prebed
grep -vw 'CDS' Ensembl.cds_utr.gtf | cut -f1,4,5 > Ensembl.utr.prebed
bedtools slop -l 1 -r 0 -i Ensembl.csd.prebed -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.csd.bed
bedtools slop -l 1 -r 0 -i Ensembl.utr.prebed -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.utr.bed
grch37_to_hg19.py Ensembl.csd.bed | grep -v "^chrH" | grep -v "^chrG" > Ensembl.csd.fixedchr.bed
grch37_to_hg19.py Ensembl.utr.bed | grep -v "^chrH" | grep -v "^chrG" > Ensembl.utr.fixedchr.bed
bedtools slop -b 50 -i Ensembl.csd.fixedchr.bed -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.csd.fixedchr.slop50.bed
cat Ensembl.utr.bed Ensembl.csd.fixedchr.slop50.bed > Ensembl.cds_utr.bed
cd ..
cat padded_panels/*.bed transcripts/Ensembl.cds_utr.bed > joined.cds_utr.bed
sort -k1,1 -k2,2n joined.cds_utr.bed > joined.cds_utr.sorted.bed
bedtools merge -i joined.cds_utr.sorted.bed > union.cds_utr.bed

# subtruct
bedtools subtract -b union.bed -a ../AZ-Exome.bed > subtraction.bed

# get the bed for annotation
awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed

# intersect
bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
cut -f8 subtraction.annotated.bed | sort -u

bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
cut -f7 subtraction.back.annotated.bed | sort -u
