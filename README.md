# AZ reporting suite

Contents:
- Post-processing and reporting for [BCBio-nextgen](https://github.com/chapmanb/bcbio-nextgen)
- Pre-alignment processing and reporting


## Post-processing and reporting for BCBio-nextgen

- Read and alignment quality via Qualimap, AZ TargQC and FastQC at capture, exon and CDS level
- Variant (SNP, Indel and SV) annotation, QC, filtration and prioritization
- CNV caller Seq2C
- RNAseq QC, expression reporting
- Variant, SV, CNV, coverage reporting via single, easy to interpret HTML Report

### Installation
Requirements:
- python 2.7
- bcbio
- virtualenv
- git
```
$ virtualenv virtualenv -p $BCBIO/0.9.7/rhel6-x64/anaconda/bin
$ git clone https://github.com/AstraZeneca-NGS/Reporting_Suite.git $AZ_REPORTING
$ export LD_LIBRARY_PATH=$BCBIO/0.9.7/rhel6-x64/anaconda/lib:$LD_LIBRARY_PATH
$ source activate $AZ_REPORTING/virtualenv/bin/activate
$ pip install -r $AZ_REPORTING/python_requirements.txt
$ R
> install.packages("DESeq2")
> install.packages("RColorBrewer")
> install.packages("gplots")
> install.packages("amap")
> install.packages("ggplot2")
> install.packages("pheatmap")
> install.packages("reshape2")
```

In Waltham, UK and Sweden AZ HPC, it is already installed under "module load bcbio". The locations are the following:

Waltham
```
BCBIO=/group/ngs/src/bcbio-nextgen
AZ_REPORTING=/group/ngs/src/az.reporting
```
UK
```
BCBIO=/apps/bcbio-nextgen
AZ_REPORTING=/ngs/RDI/SCRIPTS/az.reporting
```
Sweden
```
BCBIO=/ngs/apps/bcbio-nextgen
AZ_REPORTING=/ngs/apps/az.reporting
```
China
```
AZ_REPORTING=/ngs/RDI/SCRIPTS/az.reporting
```

### Usage
Run BCBio in a standard manner, make sure the folder contains at least these non-empty directories:
```
config/
final/
```
Then load the script path:
```
module load bcbio-nextgen
az-reporting.py [/path/to/bcbio/final/] [--jira https://jira.rd.astrazeneca.net/browse/NGSG-351]
```

The tool reads the bcbio-nextgen YAML configuration file inside the `config` directory in order to extract the information on samples, batches, variant callers.

#### BED file
As a default BED file for QC reports and Seq2C, the tool takes `cnv_regions` BED file from the bcbio-nextgen config. If `cnv_regions` is not provided, the tools takes RefSeq CDS located at ```/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/RefSeq_CDS.hg19-chr20.bed```. Optionally, the script can be fed with a custom BED file using the --bed option, e.g.:
```
bcbio_postproc.py --bed SureSelect_Human_AllExon_V5.bed
```
The BED file is automatically annotated with gene names based on RefSeq features, first based on intersections with CDS, then for everything remaining - with UTRs and ncRNA exons, then with introns, and finally with other genes: http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Annotating+a+BED+file+with+HUGO+genes+symbols

If the BED files contains 8 and more columns, it is assumed to to contain amplicon information for Seq2C.

#### SGE (decoupling in progress over to iPython)
The pipeline uses qsub to submit jobs. The qsub command line template is akin to (however does not need to be run manually by the user):
```
qsub -pe smp [threads] -S /bin/bash -q batch.q -j n -o [OUT] -e [LOG] -N [NAME] runner_Waltham.sh "[cmdline]"
```

#### Configuration files
The tool uses a project specific configuration file to tune up the post-processing steps and parameters on variant annotation, filtration, etc. The default parameters are stored in the `RUNINFO_DEFAULTS.yaml` file (located in the `configs` directory in the source root. eg. `/group/ngs/src/az.reporting/configs/run_info_ExomeSeq.yaml`). After the first run, `run_info.yaml` is put into the `config` directory. It is almost empty, meaning that all the parameters used are default. We recommend modifying this file if needed, using `RUN_INFO_DEFAULTS.yaml` as a reference - just copy lines and edit. The file in `config/run.yaml` will be used as a configuration for the subsequent re-runs. Also, as an alternative, a config file can be explicitly fed to the tool using the `--run-cnf` option.

The are 3 pre-set configurations (stored in the ```configs``` directory in the source root):

- `run_info_ExomeSeq.yaml` (default, also used for hybrid capture)
- `run_info_WGS.yaml` (used automatically if no `variant_regions` are specified in the bcbio.yaml, or when the `--wgs` flag is set)
- `run_info_DeepSeq.yaml` (used automatically if the `--deep-seq` flag is set)

For different kinds of analysis, different allele frequencies thresholds are recommended. The default values are the following.
- Exome and hybrid capture: call call with 0.01, filter with 0.05 (tissue) or 0.01 (plasma)
- WGS: call with 0.01, filter with 0.075 (tissue) or 0.01 (plasma)
- Targeted deep sequencing: call with 0.001, filter with 0.005

But if otherwise is requested, use the `-f` parameter to forse another value. You can also modify it in the run_info yaml.

The tool also uses system configuration file with the paths to external tools and reference data. For Waltham, UK, Sweden, China and Cloud, it has presets that are set automatically. For other systems, you can make you own system config based on one of the predefined AstraZeneca YAMLs ([system_info_Waltham.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/system_info_Waltham.yaml), [system_info_AP.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/system_info_AP.yaml) or [system_info_Sweden.yaml](https://github.com/AstraZeneca-NGS/Reporting_Suite/blob/master/configs/system_info_Sweden.yaml)), and provide it with the `--sys-cnf` option.

#### Altering the step list
In `run_info.yaml`, the pipeline step list is specified, and by default looks the following:
```
steps:
- Variants
- TargQC
- Seq2C
```

The step list is customizable, for example,
```
steps:
#- Variants
#- TargQC
- Seq2C
```
will run only Seq2C.
```
steps:
#- Variants
- VarFilter
#- TargQC
#- Seq2C
```
Will run only re-filter annotated VCF files and re-generate the `*.txt` and `*.PASS.txt` files in the `final/YYYY-MM-DD_projectname/` folder, as well as `*.anno.filt.vcf.gz` for each sample, and afterwards re-calculate the VarQC stats for post-filtered variants, and re-build the project-level and sample reports.

If you want only re-generate only summary project-level reports, just comment out all steps:
```
steps:
#- Variants
#- TargQC
#- Seq2C
```

#### Email notification
When all jobs are finished, an email will be sent with the path to the project HTML report, and list of all errors appeared during the processing. You can provide your email-address with the `--email` option, e.g.:
```
az-reporting.py --email Vlad.Saveliev@astrazeneca.com
```
The default SMTP server address is `localhost` and can be modified in the bottom of the system_info config file.

#### Exposing to the NGS webserver (Waltham and UK only)
When the processing is finished, a record with basic information about project is added to `/ngs/oncology/NGS.Projects.csv` (US) or `/ngs/oncology/reports/NGS.Project.csv` (UK). The information includes:
- Project name 
- Path to project-level HTML report,
- Jira url (if `--jira` is provided),
- Information automatically extracted from the Jira case: reporter username, assignee username, project type.
This makes the project shown on the [NGS website](ngs.usbod.astrazeneca.net) (see [SOP](http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Project+Dashboard) for details).

By default, the project name for reports, email notifications, and the NGS webserver list is extracted from the project path, e.g. `/ngs/oncology/Analysis/dev/Dev_0123_WGS_MiSeq_13_BN/bcbio/final` becomes `Dev_0123_WGS_MiSeq_13_BN_bcbio`. You can forse a different name using the --project option:
```
bcbio_postproc.py --project Dev_0123_WGS_MiSeq_13_BN
```
#### Output
###### Mutations
In `final/YYYY-MM-DD_projectname/`, you can find *.PASS.txt files with all prioritized mutations passed the filtering.
The annotated VCF variant files can be found in the individual sample folders `final/<sample>/` as `<sample>-vardict.anno.filt.vcf.gz`

You can find detailed information on mutation filtering/prioritization/classification in [this article](https://github.com/AstraZeneca-NGS/Reporting_Suite/tree/master/source/variants/variant_filtering.md).

###### CNV
CNV
Seq2C CNV calls are located in the `final/YYYY-MM-DD_projectname/cnv/` folder.

###### Coverage
The per-gene, per-exons and per-amplicon coverage reports can be found in the individual sample folders under `final/<sample>/targetSeq/<sample>.targetSeq.details.gene.tsv`

###### QC
The folder `final/YYYY-MM-DD_projectname/` contains the project-level summary HTML landing page `<projectname>.html` that contains links to the following reports:
  - NGS oncology reports with mutations, CNV, SV, and coverage stats and visualizations
  - VarQC (variant calling stats)
  - TargQC (target coverage stats and plots)
  - FastQC (reads stats)
  - RNAseq QC and expression heatmaps (for RNAseq projects only)

#### Troubleshooting
Logs are stored in `final/<datestamp>/log/reporting/`. The main log is called `log.txt`
- Issues with TargQC reports? May be caused by BAM files located under `final/<sample>/<sample>-ready.bam`, or a BED file, or memory issues. Look at TargQC reports logs: `final/<datestamp>/log/reporting/<sample>/targQC.log`
- Issues with variants? Look at the VCF files located under `final/<datestamp>/var/raw/`, also check annotation logs `final/<datestamp>/log/reporting/<sample>/varAnnotate-<caller>.log` and filtering logs `final/<datestamp>/log/reporting/<sample>/varFilter-<caller>.log`.
- Issues with Seq2C? Look at the BAMs, BEDs, logs `final/<datestamp>/log/reporting/Seq2C.log`, and also take a look at the coverage file located in `final/<datestamp>/cnv`.

[bcbio-nextgen]:bcbio-nextgen.readthedocs.org
[Qualimap]:qualimap.bioinfo.cipf.es
[VarAnnotate]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+Annotation
[VarQC]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+QC
[VarFilter]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+Filtration
[TargQC]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Targeted+Reseq+Reports
[Seq2C]:http://wiki.rd.astrazeneca.net/display/caninfra/Seq2C+for+copy+number+analysis


## Pre-alignment processing and reporting
This is to be run after the basecalling step. Make sure the bcl2fastq software was ran for the project (see BaseCalling / BaseCalling for HiSeq4000).

The following script detects the sequencing projects and samples from the SampleSheet.tsv file and the project structure (automatically detecting the type by base directory name: HiSeq / MiSeq / HiSeq4000). Then it does the following steps:

- Get the project name from the JIRA case, and create a symlink in the datasets folder as datasets/<division>/<project_name>
- If the fastq files are not yet merged by lanes, it merges them and saves into:
  - `Unalign/<project>/fastq` (for HiSeq4000)
  - `Unalign/Project_<project>/fastq` (for HiSeq)
  - `Data/Intensities/BaseCalls/fastq` (for MiSeq)
- Runs FastQC for each merged fastq file (separately for left and right reads), then combines them into one per project and saves into fastq/FastQC, and exposes both combined and separate htmls to the NGS webserver.
- Downsamples FastQC files to random 500.000 read pairs, aligns them to the reference, then runs and exposes the TargQC reports on resulting BAM files.
- Exposes the bcl2fastq reports:
  - `Unalign/Reports/html/index.html` (for MiSeq and HiSeq4000)
  - `Unalign/Basecall_Stats_<id>/All.htm`, `Unalign/Basecall_Stats_<id>/IVC.htm`, `Unalign/Basecall_Stats_<id>/Demultiplex_Stats.htm`, (for HiSeq).

### Usage
```
module load bcbio-nextgen
preproc \
   /ngs/oncology/datasets/HiSeq/150612_D00443_0168_AHMNFGADXX \
   [--jira https://jira.rd.astrazeneca.net/browse/NGSG-313] \
   [--bed target.bed] \
   [--project-name Dev_0104_HiSeq_DS] \
   [--expose-only]
```

Arguments:
- First argument - path to the dataset (which contains 'Unalign' folder).
- `--jira`: JIRA case URL.
- `--bed`: optional, used for downsampled TargQC reporting.
- `--project-name`: optional, by default it is retrieved from JIRA.
- `--expose-only`: Skip merging, FastQC, downsampling, TargQC. Only expose available reports to the NGS webserver.

### Multiple-project run
```
preproc \
   /ngs/oncology/Datasets/HiSeq4000/151222_K00172_0042_BH5GTFBBXX \
   --conf conf.csv
```

conf.csv
```
#project_name,desired_project_name,jira_url,bed
AZ300.Crown_Rerun,Dev_0138_HiSeq4000_CrownModels_Rerun,https://jira.rd.astrazeneca.net/i#browse/NGSG-440,target.bed
HyperPlus.Alpha_WGS,Dev_0165_HiSeq4000_AKT.WGS3,https://jira.rd.astrazeneca.net/i#browse/NGSG-436,
```

Optionally, you can process a subproject separately by providing the exact subproject path in Unalign:
```
preproc \
   /ngs/oncology/datasets/HiSeq/150612_D00443_0168_AHMNFGADXX/Unalign/HyperPlus.Alpha_WGS \
   [--jira https://jira.rd.astrazeneca.net/i#browse/NGSG-436] \
   [--project-name Dev_0165_HiSeq4000_AKT.WGS3]
```

### UK specific notes
Run bcl2fastq as 'sbsuser' first: http://wiki.rd.astrazeneca.net/display/NG/NextSeq+500+-+Basecalling+UK

As 'sbsuser', run
```
screen
module load bcbio-nextgen
/ngs/RDI/SCRIPTS/az.reporting/preproc.py /ngs/oncology/datasets/NextSeq500/160318_NB501188_0011_AH3TCVAFXX --jira https://jira.rd.astrazeneca.net/browse/NGSG-536 --genome hg38 --project-name TS_UK_0010_Bladder_2_preqc
```
where the date stamped folder, Jira URL and project-name must be changed accordingly.

--------------------------------

The following scripts are used inside the pre-processing, but can be used outside as well.

1. Running FastQC on qsub:

```
/group/ngs/src/az.reporting/scripts/pre/fastqc.py \
   -1 Colo-829_R1.fastq.gz -2 Colo-829_R2.fastq.gz \
   -o output/dir/path [--sample Colo-829] [--downsample-to 1e7]
```

2. Downsampled TargQC on qsub (fastq paired will be automatically detected based on _1 _2 _R1 _R2 suffixes:
```
/group/ngs/src/az.reporting/scripts/pre/downsampled_targqc_all.py \
   *.fq.gz -o output/dir/path --genome hg19 [--bed target.bed] [--downsample-to 1e7]
```

3. Combining FastQC reports:
```
/group/ngs/src/az.reporting/scripts/pre/fastqc_combine.py \
   /ngs/scratch/vlad/HiSeq/150612_D00443_0168_AHMNFGADXX/Unalign/fastq/fastqc/*.html \
   -o /ngs/scratch/vlad/HiSeq/150612_D00443_0168_AHMNFGADXX/Unalign/fastq/fastqc/FastQC.html
```

