Target and WGS sequencing coverage analysis, variant annotation and quality control after alignment and variant calling.

Includes:
- Read and alignment quality via Qualimap, NGSCat, AZ TargQC and FASTQC at probe, exon and CDS level
- variant (SNP, Indel and SV) annotation
- variant (SNP, Indel and SV) quality control
- variant (SNP, Indel and SV) filtration and prioritization
- CNV caller Seq2C
- Variant, CNV, Coverage reporting via single, easy to interpret HTML Report

<br>
####Installation
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
#
#
git clone https://github.com/AstraZeneca-NGS/Reporting_Suite.git
cd $AZ_REPORTING
export LD_LIBRARY_PATH=$BCBIO/0.9.7/rhel6-x64/anaconda/lib:$LD_LIBRARY_PATH
module load virtualenv
virtualenv virtualenv -p $BCBIO/0.9.7/rhel6-x64/anaconda/bin
pip install -r python_requirements.txt
```

$R
install.packages("DESeq2")
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("amap")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("reshape2")

<br>
####Usage (bcbio decoupling in progress)
```
export LD_LIBRARY_PATH=$BCBIO/0.9.7/rhel6-x64/anaconda/lib:$LD_LIBRARY_PATH
source $AZ_REPORTING/virtualenv/bin/activate
az-reporting.py [/path/to/a/bcbio/directory] [--run-cnf run_info.yaml] [--sys-cnf system_info.yaml] [--bed target.bed]
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
- [VarAnnotate]
- [VarQC]
- [VarFilter]
- VarQC_postVarFilter (VarQC for the variants passed filters)
- [TargetCov]
- [Seq2C] (copy number estimates)
- [ngsCAT]
- [QualiMap]
- MongoLoader (disabled by default, turned on with the \--load-mongo option)
- AbnormalCovReport (disabled by default)

Afterwards, the tool generates project-level summary reports located in a ```final/YYYY-MM-DD/qc``` directory:
- VarQC summary (based on VarQC),
- TargQC summary (based on TargetCov, ngsCAT and QualiMap)
- Fastqc summary (based on individual Fastqc reports from the bcbio run)
- CombinedReport (single project-level summary for all samples and all summary reports)

The step list is customisable at ```run_info.yaml```. For example,
```
steps:
#- VarAnnotate
#- VarQC
#- VarFilter
#- VarQC_postVarFilter
#- TargetCov
#- Seq2C
- ngsCAT
#- QualiMap
#- MongoLoader
#- AbnormalCovReport
```

will run only ngsCAT for each sample, and afterwards generate summary.

<br>
The ```--load-mongo``` option adds the MongoLoader step - in Waltham, it loads final filtered variants into Mongo database:

```
az-reporting.py --load-mongo
```

<br>
####Custom BED file

You can specify a custom BED file for target QC (TargetCov, ngsCAT, QualiMap), and Seq2C:

````
az-reporting.py --bed /ngs/reference_data/genomes/Hsapiens/hg19/bed/Agilent/SureSelect_Human_AllExon_V5.bed
````

[bcbio-nextgen]:bcbio-nextgen.readthedocs.org
[ngsCAT]:www.bioinfomgp.org/ngscat
[Qualimap]:qualimap.bioinfo.cipf.es
[VarAnnotate]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+Annotation
[VarQC]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+QC
[VarFilter]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Variant+Filtration
[TargetCov]:http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Targeted+Reseq+Reports
[Seq2C]:http://wiki.rd.astrazeneca.net/display/caninfra/Seq2C+for+copy+number+analysis
