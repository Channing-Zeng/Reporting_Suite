## Mutation prioritization in [AZ post-processing pipeline](https://github.com/AstraZeneca-NGS/Reporting_Suite) for BCBio-nextgen

### Variant calling
Variant calling is performed by [VarDict](https://github.com/AstraZeneca-NGS/VarDict) ([Lai Z, 2016](http://www.ncbi.nlm.nih.gov/pubmed/27060149)) that calls SNV, MNV, indels, complex and structural variants, and differences in somatic and LOH variants when processed in paired mode. 

##### Allele frequency threshold
The allowed allele frequency is 1% for exomes, hybrid capture, and whole genome sequencing, and .1% for deep sequencing. 

##### Target regions
In exomes, the variant are called in a special target called AZ exome, that combines all Ensembl CDS regions, UTR, plus regions from commonly used panels padded by 50pb:
```
ExomeSNP_ID.bed
FM_T5.bed
IDT_Exome.bed
IDT-PanCancer_AZ1-IDT_orig.bed
IDT_PanCancer_Exons.bed
Illumina_Nextera_Exome.bed
Illumina_TruSeq_Exome.bed
Illumina_TruSight_Cancer.bed
Personalis.bed
SeqCap_EZ_Exome_v3_capture.bed
SureSelect_Human_AllExon_V4.bed
SureSelect_Human_AllExon_V5.bed
Xgen-PanCancer.bed
```
For WGS projects, BCBio-nextgen removes tricky regions like centromers from the target.

### Annotation
Variants in form of VCF are annotated using [SnpEff](http://snpeff.sourceforge.net/) tool that predicts effect of variants on proteins according to RefSeq gene model. It considered canonical (longest) transcripts only, except for the following genes:
```
FANCL   NM_018062.3
MET     NM_000245.2
CDKN2A  NM_000077.4
BRCA1   NM_007294.3
MYD88   NM_002468.4
PPP2R2A NM_002717.3
RAD51D  NM_002878.3
RAD54L  NM_003579.3
ESR1    NM_000125.3
AKT1    NM_005163.2
FGFR3   NM_000142.4
CD79B   NM_000626.2
CHEK2   NM_007194.3
CHEK1   NM_001274.5
```
The full list is located under `/ngs/reference_data/genomes/Hsapiens/hg19/canonical_transcripts.txt` (hg19) and `/ngs/reference_data/genomes/Hsapiens/hg38/canonical_transcripts.txt` (hg38).

SnpEff assigns:
- gene and transcript IDs
- gene and transcript biotypes: coding/ncRNA/pseudogene
- region: coding/splice/intron/upstream/downstream/non-coding
- functional class: silent/missence/nonsence
- codon change
- aminoacid change

Variants are also searched against the following mutation databases:
- [COSMIC](http://cancer.sanger.ac.uk/cosmic) - cancer somatic mutations database, assigns ID and hits count
- [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) - assigns rsID and CAF (global allele frequencies)
- [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - assigns CLNSIG (clinical significance)

### Raw filtering: vardict.txt
The first filtering step is performed using [vcf2txt.pl](https://github.com/AstraZeneca-NGS/VarDict/blob/master/vcf2txt.pl) script from the VarDict package. Its input is a VarDict VCF file annotated as described above, its output is `vardict.txt` sample-level file. The program (1) performs hard and soft filtering for low quality variants, (2) assigns _variant class_ ("novelty"), (3) assigns _variant type_ (CNV, MNV, deletion, insertion, compelex)

Hard filtering (where the variants are discarded) is performed using the following parameters:
- Locus total depth (≥ 5x for exomes, hybrid capture, and targeted)
- Mean position in reads (≥ 5)
- Mean base quality phred score (≥ 25)

Soft filtering (the variants are reproted into `vardict.txt` with a reject reason in the `PASS` column) is done based on the following:
- Variant depth (≥ 2x for WGS; ≥ 3x for exomes, hybrid capture, and targeted)
- Alelle frequency (≥ 7.5% for whole genome sequencing; ≥ 5% for exomes and hybrid capture; ≥ .5% for targeted)
- Mean mapping quality (≥ 10)

The thresholds are specified in the run_info.yaml configuration file. Depending on the analysis, it can be either of:
- [run_info_ExomeSeq.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/run_info_ExomeSeq.yaml)
- [run_info_WGS.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/run_info_WGS.yaml)
- [run_info_DeepSeq.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/run_info_DeepSeq.yaml)

The defaults are pulled from [RUNINFO_DEFAULTS.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/tree/master/configs/RUNINFO_DEFAULTS.yaml)

Mutation class (`Var_Class`) is assigned in the following order:
- _COSMIC_ - any mutation found in COSMIC
- _ClnSNP_known_ - any other mutation, labeled significant in ClinVar (3 < `CLNSIG` < 7)
- _dbSNP_del_ - deletion found in dbSNP
- _dbSNP_ - the remaining dbSNP variants
- _Novel_ - all remaining variants

The results of the script are saved under `final/YYYY-MM-DD_projectname/var/vardict.txt`

### Cancer mutation filtering: vardict.PASS.txt
This step is mostly based on [vardict2mut.pl](https://github.com/AstraZeneca-NGS/VarDict/blob/master/vardict2mut.pl) script from the VarDict package, and partly on Rob McEven, Hedley Carr, and Brian Dougherty knowledge. Its input is a sample-level `vardict.txt`, its output is the sample-level `vardict.PASS.txt`.

This script does hard-filtering, and classified the remaining variants into the following classes:
- known (actionable; tier 1; clinically significant)
- likely (not actionable but high-impact; suppressors LOF; tier 2; Compendia and COSMIC hotspots; dbSNP deletions)
- unknown (remaining variants that were not hard-filtered)
- silent (any variant that would fall into known/likely/unknown if it weren't silent)
- incidentalome (found in blacklisted low-complexity genes and regions)

##### Raw hard filtering
The following records are removed:
- all soft-filtered variants in the previous step
- all variants reported as `protein_protein_contact` according to SnpEff

##### Actionability
The variant is checked against a set of rules that defined "actionable" (known driver) variants. A rules may specify any specific variant feature out of the following:
- gene
- exon 
- genomic position and change 
- protein position and change
- genomic region
- protein region
- indel type (deletion, frameshift deletion, insertion, frameshift insertion, indel, etc.)

For position-based rules, there exist rules for somatic and germline actionable variants.

For actionable variants, a special allele frequency (AF) threshold is defined. For somatic variants, by default it corresponds to the AF theshold used for variant calling (1% for exomes, hybrid capture, and whole genome sequencing, and .1% for deep sequencing); for germline, it is set to 15%. If the variant does not pass the threshold, it is not reported as actionable.

##### Filter out rules

##### Known
Robust evidence base linking aberration with likely sensitization or high likelihood of response in human in proposed trial and disease setting – good enough to go into a PhIII pivotal trial as a selection marker (registration intent). Well characterized pathogenic germline variants may be included. 
Similar to: Class 5 – pathogenic.

/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics/actionable.txt

##### Likely
frameshifts, stop gained, start lost, splice site, compendia ms7 hotspots, dbSNP deletions, or tier 2 according to your rules, or COSMIC hotspots with count more than 5

##### Rules

##### GMAF

##### MSI

##### Incidentalome

### Cohort sequencing artefacts filtering
If the samples are not homogeneous, but come from one sequencer's run, we expect recurring variants to be caused by sequencing artefacts. For all not-known and not-actionable variants, we calcualte the number and the percentage of samples harbouring this variant. If it's more than 40% of all samples and at the same time at least 5 samples, such variants are filtered out. Cohort filtering is done on the stage of merging all sample-level `vardict.PASS.txt` into the project level `vardict.PASS.txt` located in `final/<datestamp>` directory.

### Reporting

### Conversion to VCF


Germline variants – germline SNPs occur at approximately 100%, 50%, or 0% frequency and every effort is made to filter out germline variants (some exceptions like BRCA1/2).

Sample types - in general, focus on somatic variants above ~30% allele frequency for cell line and explant models, above ~5% and above for tumors, and above ~0.7% for ctDNA. Formalin treated tumors have the added complication of C>T and G>A de-amination artefacts.

Somatic variants – somatic variants occur at nearly any allele frequency but as allele frequencies approach 20% and lower, data become noisier and false positive variants increase.

False positive mutations – advance algorithms were employed to generate the highest quality variant calls, but false positive will still occur and are likely present in the candidate mutation output below. False positive variants are more found in homopolymer tracts, in genes with repetitive regions (within the gene or elsewhere in the genome), in very large genes, and in genes residing in repetitive or unstable regions of the genome.


