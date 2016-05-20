AZ post-processing pipeline for [link](https://github.com/chapmanb/bcbio-nextgen "BCBio-nextgen"): coverage QC, variant annotation, variant filtration, clinical reporting.

Includes:
- Read and alignment quality via Qualimap, AZ TargQC and FastQC at capture, exon and CDS level
- Variant (SNP, Indel and SV) annotation
- Variant (SNP, Indel and SV) quality control
- Variant (SNP, Indel and SV) filtration and prioritization
- CNV caller Seq2C
- Variant, CNV, Coverage reporting via single, easy to interpret HTML Report

<br>
####Installation
```
python 2.7 required
bcbio required
virtualenv virtualenv -p $BCBIO/0.9.7/rhel6-x64/anaconda/bin
Clone repo https://github.com/AstraZeneca-NGS/Reporting_Suite.git into $AZ_REPORTING
export LD_LIBRARY_PATH=$BCBIO/0.9.7/rhel6-x64/anaconda/lib:$LD_LIBRARY_PATH
cd $AZ_REPORTING
source activate virtualenv/bin/activate
pip install -r python_requirements.txt
R packages:
install.packages("DESeq2")
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("amap")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("reshape2")
```

#US:
BCBIO=/group/ngs/src/bcbio-nextgen
AZ_REPORTING=/group/ngs/src/az.reporting
#UK:
BCBIO=/apps/bcbio-nextgen
AZ_REPORTING=/ngs/RDI/SCRIPTS/az.reporting
#Sweden:
BCBIO=/ngs/apps/bcbio-nextgen
AZ_REPORTING=/ngs/apps/az.reporting

<br>
####Usage (bcbio decoupling in progress)
```
module load bcbio-nextgen
az-reporting.py [/path/to/a/bcbio/directory] [--run-cnf run_info.yaml] [--bed target.bed]
```

Instead of specifying the first argument, you can change to bcbio-nextgen project directory:

```
cd /path/to/a/bcbio/directory
az-reporting.py
```

The tool reads the bcbio-nextgen YAML configuration file inside the ```config``` directory in order to extract the information on samples, batches, variant callers, and final directory name.

<br>
####SGE (SGE decoupling in progress over to iPython)

The pipeline uses qsub to submit jobs. The qsub command line template is the following:

```
qsub -pe smp [threads] -S /bin/bash -q batch.q -j n -o [OUT] -e [LOG] -N [NAME] runner_Waltham.sh "[cmdline]"
```

<br>
####Configuration files
The tool also can be provided with its own configuration file to tune up the post-processing steps and parameters on variant annotation, filtration, etc. Default run config is ```run_info.yaml```, and it gets copied to the ```config``` directory with the first launch. You can edit it later and restart the pipeline with different parameters; or you provide your own configuration file with ```--run-cnf``` option.

The tool also uses system configuration file with the paths to external tools, reference data. You can make you own system config based on one of the predefined AstraZeneca YAMLs (```system_info_Waltham.yaml``` or ```system_info_AP.yaml```) and provide it with the ```--sys-cnf``` option.

<br>
####Email notification

You can provide your email-address with the ```--email``` option, e.g.:

```
az-reporting.py --email Vlad.Saveliev@astrazeneca.com
```

After all the jobs to GRID are submitted, the tool sends an email with the lsit of log files to track the jobs. When all jobs are finished, one more email is sent.
The default SMTP server address is ```localhost``` and can be modified in the bottom of the system_info config file.

<br>
####Altering the step list
For each sample in the list, the tools runs the following steps:
- [Variants]
- [TargQC]
- [Seq2C] (copy number estimates)

Afterwards, the tool generates project-level summary reports located in a ```final/YYYY-MM-DD/qc``` directory:
- VarQC summary (based on VarQC),
- TargQC summary (based on TargetCov, ngsCAT and QualiMap)
- Fastqc summary (based on individual Fastqc reports from the bcbio run)
- CombinedReport (single project-level summary for all samples and all summary reports)

The step list is customisable at ```run_info.yaml```. For example,
```
steps:
#- Variants
- VarFilter
#- TargQC
#- Seq2C
```

will run only filtering for each sample, and afterwards generate summary.

<br>
####Custom BED file

You can specify a custom BED file for target QC (TargetCov, ngsCAT, QualiMap), and Seq2C:

````
az-reporting.py --bed /ngs/reference_data/genomes/Hsapiens/hg19/bed/Agilent/SureSelect_Human_AllExon_V5.bed
````

<br>
####Workflow


[bcbio-nextgen]:bcbio-nextgen.readthedocs.org
[Qualimap]:qualimap.bioinfo.cipf.es
[VarAnnotate]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+Annotation
[VarQC]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+QC
[VarFilter]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+Filtration
[TargQC]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Targeted+Reseq+Reports
[Seq2C]:http://wiki.rd.astrazeneca.net/display/caninfra/Seq2C+for+copy+number+analysis


