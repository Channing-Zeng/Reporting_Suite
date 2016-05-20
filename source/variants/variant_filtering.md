## Mutation prioritization in [AZ post-processing pipeline](https://github.com/AstraZeneca-NGS/Reporting_Suite) for BCBio-nextgen

### Variant calling
Variant calling is performed by [VarDict](https://github.com/AstraZeneca-NGS/VarDict) ([Lai Z, 2016](http://www.ncbi.nlm.nih.gov/pubmed/27060149)) in BCBio-nextgen, using the allele frequency threshold of 1% for exomes, hybrid capture, and whole genome sequencing, or .1% for targeted deep sequencing. Samples analysed either in paired (tumor vs. normal) or single sample mode.

### Annotation
Variants are annotated using [SnpEff](http://snpeff.sourceforge.net/) tool that predicts effect of variants on proteins with the reference to the RefSeq gene model. It considered canonical (longest) transcripts only, except for the following genes:
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

Variants are also searched against:
- [COSMIC](http://cancer.sanger.ac.uk/cosmic) - cancer somatic mutations database, assigns ID and hits count
- [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) - assigns rsID and CAF (global allele frequencies)
- [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - assigns CLNSIG (clinical significance)

### First step filtering: vardict.txt
First step is performed using [vcf2txt.pl](https://github.com/AstraZeneca-NGS/VarDict/blob/master/vcf2txt.pl) script from the VarDict package. The program (1) performs hard and soft filtering for low quality variants, (2) assigns _variant class_ ("novelty"), (3) assigns _variant type_ (CNV, MNV, deletion, insertion, compelex)

Hard filtering - the variants are discarded - is performed based on the following parameters:
- Locus total depth (5x for exomes, hybrid capture, and targeted)
- Mean position in reads (5)
- Mean base quality phred score (25)

Soft filtering - the variants are kept in vardict.txt, but with PASS!=True - is done based on the following:
- Variant depth (2x for WGS; 3x for exomes, hybrid capture, and targeted)
- Alelle frequency (7.5% for whole genome sequencing; 5% for exomes and hybrid capture; .5% for targeted)
- Mean mapping quality (10)

The thresholds are specified in the run_info.yaml configuration file. Depending on the analysis, it can be either of [run_info_ExomeSeq.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/run_info_ExomeSeq.yaml), [run_info_WGS.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/run_info_WGS.yaml), or [run_info_DeepSeq.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/run_info_DeepSeq.yaml); the defaults are pulled from [RUNINFO_DEFAULTS.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/tree/master/configs/RUNINFO_DEFAULTS.yaml).

Mutation class (Var_Class) is assigned in the following order:
- _COSMIC_ - any mutation found in COSMIC
- _ClnSNP_known_ - any other mutation, significant according to ClinVar (3 < CLNSIG < 7)
- _dbSNP_del_ - any other dbSNP deletion
- _dbSNP_ - the remaining dbSNP variants
- _Novel_ - all remaining variants

The results of the script are saved under `final/YYYY-MM-DD_projectname/var/vardict.txt`

### Second step filtering: vardict.PASS.txt
This step is partly based on [vardict2mut.pl](https://github.com/AstraZeneca-NGS/VarDict/blob/master/vardict2mut.pl) script from the VarDict package, and partly on Rob McEven, Hedley Carr, and Brian Dougherty knowledge.

#### Classification

##### Classes
- incidentalome - blacklisted genes and regions
- unknown - 
- known - actionable, tier 1, ClinVar known.
- likely - not known as actinable, but potentially important for cancer

frameshifts, stop gained, start lost, splice site, compendia ms7 hotspots, dbSNP deletions, or tier 2 according to your rules, or COSMIC hotspots with count more than 5

##### Actionable
/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics/actionable.txt

#### Rules

#### GMAF

#### MSI

##### Incidentalome

#### Reporting

