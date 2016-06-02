#### v1.4.1, 3 Jun 2016
Filtering
- Gene blacklist is split into 2 parts: hard-filtered genes, and soft-filtered genes/regions (faded in the reports)
- Filter more homopolymer artefacts

NGS reports 
- Reporting mutations in regions and genes with poor mapability
- Add COSMIC and dbSNP links back
- Transcript, AA len and exon number moved to gene tooltip

Added variant filtering documentation

#### v1.4.0, 26 May 2016
Filtering
- Gene and region blacklists introduced to filter variants
- Checking if ClinVar gene matches SnpEff gene
- Cohort frequency default threshold 1.0 (no filtering) changed to 0.4
- Germline default AF is 15%
- Defailt AF is 7.5% (2.5% for actionable) for exomes, 7.5% (5%) for WGS, 0.5% (0.2%) for deep-seq
- Other default filtering parameters changed (see [configs/RUNINFO_DEFAULT.yaml])
- MSI check constants changed
- If low ClinVar significance, filter even if COSMIC
- Artefacts and actionable filtering rules updated

NGS Reports
- Unique color scheme for Seq2C report and plot
- Column for HP (homopolymer length if > 3) added
- Indicendalome tab for mutations in blacklisted genes and regions
- SolveBio links for all mutations
- Actioanble mutations shown disregarding of the AF slider
- Buttons to report mutation and send an email
- Reporting for 820 key genes instead of 300
- Showing version and link to changes
- Showing link to documentation

#### v1.3.0, 21 Apr 2016
Features
- RNA-seq bcbio workflow support. QC reporting, counts heatmaps, RNA-seq Qualimap
- Support Sweden HPC
- Updated Vardict2mut filtering:
  - Keeping silent mutations
  - Filtering by MSI
  - Filtering by GMAF, AF, cohort frequency, depth after actionable filtering
  - Ignoring variants occurring after last known critical amino acid
- Using virtualenv to handle python dependencies
- More advanced BED annotation based on RefSeq
- NGS reports:
  - Separate tab for silent mutations
  - Key genes coverage theshold for NGS reports is set exactly as half mean coverage
  - SolveBio integration
  - Dynamic AF filtration slider
- Preproc: samblaster to deduplicate
- Add script to prepare canonical transcript list from SnpEff output
- Update script for making RefSeq BED reference data
- VarQC reports: showing cosmic and dbsnp version in tooltips
- Support normal sample only configurations

#### v1.2.0, 28 Mar 2016
Features
- Circos plots integrated with all mutations visualized by frequency accorss the genome
- SV in NGS oncology reports: more consise, highlighting known fusion, linking to jBrowse.

Optimization
- New filtration features (vardict2mut):
  - more effective cohort filtering
  - cohort filtering performed only in the final step, so known and actionable mutations are not affected
  - filtering by MSI
  - more comprehensive splice site mutation detection
  - filter_artifacts hits AF threshold reduced 50% -> 20%


#### v1.1.0, 24 Mar 2016
Features
- Oncoptints integration
- jBrowse integration with UK
- New sex determining functionality
- TargQC.py can be fed with fastqc files, and accepts the   -downsample-to option

Optimization
- Seq2C optimized, more time and memory efficient. Using "sambamba depth" to call coverage

Other
- More clear Download button in NGS oncology report
- Vardict2mut filtering updated, reporting more comprehensive rationale


#### v1.0.1, 2 Mar 2016
- Stanalone variants.py - handling incorrect input
- BED processing scripts updates: annotate_bed.py, sort_bed.py, standardize_bed.py
- Using all transcripts for bed annotation, canonical only for reports and CNV. make_exons.py udpated to make both canon and all

Optimization
- TargQC memory optimization: Do not keep regions in memory. Run Picard before coverage
- Bam to bigwig converion using bedgraph as intermediate (faster, reduce memory usage)

Fixes
- More memory allocated to qsub for TargQC - fixes some crashes
- TargQC standalone fix - thanks Sally Luke (https://github.com/AstraZeneca-NGS/Reporting_Suite/issues/19)


#### v1.0.0, 29 Feb 2016
New features:
- New mutation filtering: using rules and manual curations by Robert McEwen and Hedley Carr. Reporting clasification assignment rationale.
- Standalone script for variant annotation and filtration - variants.py
- TargQC detailed coverage stats for unmerged regions.
- Using canonical RefSeq transcripts (coding and miRNA) for TargQC detailed reports. The same RefSeq canonical transcripts are used to report mutation effects, and RefSeq CDS are used for Seq2C cnv calls.
- Extended summary HTML: phenotype info, project and codebase datestamps, and versions.
- New snpEff 4.2.
- NGS oncology report has many new features: Allele frequency plot, flagged regions, exportable data, Seq2C plos for whole profile, detailed gene-level coverage plot.
- Integration and exposing mutations, regions and alignment data to jBrowse. Includes convertion to BigWig.
- Integration with LIMS to pull sample data.
- Support bcbio 0.9.6.

Optimization
- Memory efficient vcf2txt.pl (first round filtering script) - saving 80% of memory.
- Using BCFTools, which is faster and memory efficient, but also which addresses problem of compressed VCF reference data in hg19/hg38 notation.
- Using tweaked sambamba for region coverage calculation - faster and memory effecient, and also works for unmerged regions.

Fixes
- Fix issue when a intergenic/intragenic change is annotated by COSMIC as a missence change from another gene.
- Fix occasionally missing dbSNP annotations in chr1.
- Internally referring to genes by symbol+chromosome rather than by symbol, which fixes some region annotation.
