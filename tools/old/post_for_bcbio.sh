#!/bin/sh

usage()
{
cat << EOF
usage: $0 -s <SAMPLES> -d <BCBIO FINAL DIR> -v <VCF SUFFIX> -b <BED FILE>

This script runs reporting suite on the bcbio final directory.

options:
   -d    Path to bcbio-nextgen final directory (default is pwd)
   -s    List of samples (default is samples.txt in bcbio final directory)
   -b    BED file
   -v    Suffix to file VCF file (mutect, ensembl, freebayes, etc)
   -q    Run QualiMap
   -z    Verbose (prints all command lines executed)
EOF
}

DIR=`pwd`
SAMPLES=
vcf_suf=
BED=
QUALIMAP=0
VERBOSE=0
while getopts "qzd:s:b:v" OPTION
do
     case ${OPTION} in
         q)
             QUALIMAP=1
             ;;
         d)
             DIR=$OPTARG
             ;;
         s)
             SAMPLES=$OPTARG
             ;;
         b)
             BED=$OPTARG
             ;;
         v)
             vcf_suf="-$OPTARG"
             ;;
         z)
             VERBOSE=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z ${BED} ]] || [[ -z ${vcf_suf} ]]
then
     usage
     exit 1
fi

if [[ -z ${SAMPLES} ]]
then
    SAMPLES=${DIR}/samples.txt
fi

echo "Using BED file: "${BED}
echo "Running on bcbio-nextgen final diretory: "${DIR}
echo "Using the list of samples: "${SAMPLES}
echo "Searching VCFs with the suffix "${vcf_suf}
if [ ${QUALIMAP} = 1 ] ; then
    echo 'Running QualiMap for each sample.'
else
    echo 'Skipping QualiMap.'
fi

function realpath {
    echo `python -c "import os,sys; print os.path.realpath(sys.argv[1])" "${1}"`
}

function cur_path {
    echo "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
}

BED=$(realpath ${BED})
DIR=$(realpath ${DIR})
SAMPLES=$(realpath ${SAMPLES})

function run_on_grid {
    cmdline=$1
    name=$2
    output_dir=$3
    output_file=$4
    threads=$5
    hold_jid=$6

    echo ${name}

    if [ ${VERBOSE} = 1 ] ; then echo ${cmdline}; fi

    if [ ! -z "${threads}" ]; then threads=4; fi

    if [ -f "${output_file}" ]; then rm "${output_file}"; fi
    output_file=$(realpath "${output_file}")

    runner_script=$(realpath $(cur_path)"/basic_runner.sh")

    if [ ! -z "${hold_jid}" ]; then hold_jid="-hold_jid "${hold_jid}; fi

    qsub_command="qsub ${hold_jid} -pe smp ${threads} -S /bin/bash -q batch.q -b y -o \"${output_file}\" -e ${name}.log -N ${name} bash \"${runner_script}\""
    if [ ${VERBOSE} = 1 ] ; then echo ${qsub_command}; echo ""; fi
    eval "${qsub_command}"
}

qc_jobids=""
targetcov_jobids=""


if [ ${QUALIMAP} = 1 ] ; then
    qualimap_bed=${BED}

    bed_col_num=`awk '{print NF}' ${BED} | sort -nu | tail -n 1`
    if [ ${bed_col_num} -lt 6 ]; then
        tmp_bed="/ngs/tmp/bed_fixed_for_qualimap.bed"

        if [ ${bed_col_num} -lt 5 ]; then
            awk 'NR==1 {v="\t0\t+"}{print $0v}' ${BED} > ${tmp_bed}
        else
            awk 'NR==1 {v="\t+"}{print $0v}' ${BED} > ${tmp_bed}
        fi
        qualimap_bed=${tmp_bed}
    fi
fi

indel_filter_name="indel_filter"
varannotate_name="varannotate"
varqc_name="varqc"
targetcov_name="targetcov"
ngscat_name="ngscat"
qualimap_name="qualimap"
varqc_summary_name="varqc_summary"
targetcov_summary_name="targetcov_summary"


filtered_vcf_suf=${vcf_suf}.filtered

for sample in `cat ${SAMPLES}`
do
    echo "${sample}"
    echo ""

    cd "${DIR}/${sample}"

    # rm -rf "${sample}${filtered_vcf_suf}.vcf" annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*

    if [ ! -f "${sample}${vcf_suf}.vcf" ]; then
        gunzip -c "${sample}${vcf_suf}.vcf.gz" > "${sample}${vcf_suf}.vcf"
    fi

    ### InDelFilter ###
    name=${indel_filter_name}
    cmdline="python /group/ngs/bin/InDelFilter.py \"${sample}${vcf_suf}.vcf\" > \"${sample}${filtered_vcf_suf}.vcf\""
    run_on_grid "${cmdline}" ${name}_${sample} . "${sample}${filtered_vcf_suf}.vcf" 1

    ### VarAnn ###
    name=${varannotate_name}
    cmdline="python /group/ngs/src/varannotate.py --var \"${sample}${filtered_vcf_suf}.vcf\" --bam \"${sample}-ready.bam\" -o ${name}"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 1 ${indel_filter_name}_${sample}

    ### VarQC ###
    name=${varqc_name}
    cmdline="python /group/ngs/src/varqc.py --var \"${sample}${filtered_vcf_suf}.vcf\" -o ${name}"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 1 ${indel_filter_name}_${sample}

    if [ ! -z "${qc_jobids}" ]; then
        qc_jobids=${qc_jobids},${varqc_name}_${sample}
    else
        qc_jobids=${varqc_name}_${sample}
    fi

    ### targetCov ###
    name=${targetcov_name}
    cmdline="python /group/ngs/src/targetcov.py --bam \"${sample}-ready.bam\" --bed \"${bed}\" --nt 4 -o ${name}"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 4

    if [ ! -z "${targetcov_jobids}" ]; then
        targetcov_jobids=${targetcov_jobids},${targetcov_name}_${sample}
    else
        targetcov_jobids=${targetcov_name}_${sample}
    fi

    ### NGSCat ###
    name=${ngscat_name}
    cmdline="python /group/ngs/src/ngscat/ngscat.py --bams \"${sample}-ready.bam\" --bed \"${bed}\" --out ${name} --saturation y"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 4

    ## QualiMap ##
    name=${qualimap_name}
    cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam \"${sample}-ready.bam\" -outdir ${name} -gff \"${qualimap_bed}\" -c -gd HUMAN"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 4

    cd ..
done

## VarQC summary ##
name=${varqc_summary_name}
cmdline="python /gpfs/group/ngs/src/ngs_reporting/varqc_summary.py ${DIR} ${SAMPLES} ${varqc_name} ${filtered_vcf_suf}"
run_on_grid "${cmdline}" ${name} . log_${name} 1 ${qc_jobids}

## Target coverage summary ##
name=${targetcov_summary_name}
cmdline="python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py ${DIR} ${SAMPLES} ${targetcov_name}"
run_on_grid "${cmdline}" ${name} . log_${name} 1 ${targetcov_jobids}
