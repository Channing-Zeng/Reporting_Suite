##!/bin/sh
#
#bed=$1
#bcbio_final_dir=$2
#samples=$3
#vcf_suffix=$4
#
#if [ -z "${samples}" ]; then
#    samples="./samples.txt"
#fi
#
#if [ -z "${vcf_suffix}" ]; then
#    vcf_suffix="-mutect"
#fi
#
#echo "Using BED file: "${bed}
#echo "Running on BCBIO final dir: "${bcbio_final_dir}
#echo "Using the list of samples: "${samples}
#echo "Using VCFs with the suffix "${vcf_suffix}
#echo ""
#
##filter_indels=false
##if [[ $* == *--filter-indels* ]]; then
##    filter_indels=true
##fi
#
#function run_on_grid {
#    cmdline=$1
#    name=$2
#    output_dir=$3
#    output_file=$4
#    threads=$5
#    hold_jid=$6
#
#    if [ ! -z "${threads}" ]; then
#             threads=4
#    fi
#
#    if [ -f "${output_file}" ]; then
#        rm ${output_file}
#    fi
#
#    cwd=`pwd`
#    echo $cwd
#    cd ${output_dir}
#    runner_script=${name}.sh
#    cd ${cwd}
#    runner_script=${output_dir}/${runner_script}
#    rm ${runner_script}
#}
#
#qc_jobids=""
#targetcov_jobids=""
#
#bed=`python -c "import os,sys; print os.path.realpath(sys.argv[1])" ${bed}`
#
#
#bed_col_num=`awk '{print NF}' ${bed} | sort -nu | tail -n 1`
#if [ ${bed_col_num} -lt 5 ]; then
#    tmp_bed="/ngs/tmp/bed_fixed_for_qualimap.bed"
#    rm ${tmp_bed}
#fi
#
#for sample in `cat ${samples}`
#do
#    echo "cd to ${bcbio_final_dir}/${sample}"
#    cd ${bcbio_final_dir}/${sample}
#    echo ""
#
#    cmdline="rm -rf ${sample}${vcf_suffix}.filtered.vcf annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*"
#    echo ${cmdline}
#    eval "${cmdline}"
#
#    ### InDelFilter ###
#    cmdline="python /group/ngs/bin/InDelFilter.py ${sample}${vcf_suffix}.vcf > ${sample}${vcf_suffix}.filtered.vcf"
#    run_on_grid "${cmdline}" InDelFilter_${sample} . ${sample}${vcf_suffix}.filtered.vcf 1
#    vcf_suffix=${vcf_suffix}.filtered
#
#    cd ..
#done
#
### VarQC summary ##
#cmdline="python /gpfs/group/ngs/src/ngs_reporting/varqc_summary.py $bcbio_final_dir $samples varQC ${vcf_suffix}"
#run_on_grid "${cmdline}" VarQCSummary . ../work/log_varqc_summary 1 ${qc_jobids}
#
### Target coverage summary ##
#cmdline="python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py $bcbio_final_dir $samples targetSeq"
#run_on_grid "${cmdline}" targetSeqSummary . ../work/log_targetcov_summary 1 ${targetcov_jobids}
