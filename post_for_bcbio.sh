#!/bin/sh

bed=$1
bcbio_final_dir=$2
samples=$3
vcf_suffix=$4

if [ -z "${samples}" ]; then
    samples=${bcbio_final_dir}"/samples.txt"
fi

if [ -z "${vcf_suffix}" ]; then
    vcf_suffix="-mutect"
fi

echo "Using BED file: "${bed}
echo "Running on BCBIO final dir: "${bcbio_final_dir}
echo "Using the list of samples: "${samples}
echo "Using VCFs with the suffix "${vcf_suffix}
echo ""

#filter_indels=false
#if [[ $* == *--filter-indels* ]]; then
#    filter_indels=true
#fi

function run_on_grid {
    cmdline=$1
    name=$2
    output_dir=$3
    output_file=$4
    threads=$5
    hold_jid=$6

    if [ ! -z "${threads}" ]; then
         threads=4
    fi

    if [ -f "${output_file}" ]; then
        rm "${output_file}"
    fi
    output_file=`python -c "import os,sys; print os.path.realpath(sys.argv[1])" "${output_file}"`

    echo ${cmdline}

    cwd=`pwd`
    cd "${output_dir}"
    runner_name="${name}.sh"
    runner_script=`python -c "import os,sys; print os.path.realpath(sys.argv[1])" ${runner_name}`
    echo "#!/bin/bash" > "${runner_script}"
#    echo "date" >> "${runner_script}"
## TODO: module loads could print to output
    echo "source /etc/profile.d/modules.sh" >> "${runner_script}"
    echo "module load python/64_2.7.3 java bedtools samtools" >> "${runner_script}"
    echo ${cmdline} >> "${runner_script}"
#    echo "date" >> "${runner_script}"
    echo "rm -- \"${runner_script}\"" >> "${runner_script}"
    chmod +x "${runner_script}"
    cd ${cwd}

    if [ ! -z "${hold_jid}" ]; then
        hold_jid="-hold_jid "${hold_jid}
    fi

    qsub_command="qsub ${hold_jid} -pe smp ${threads} -S /bin/bash -cwd -q batch.q -b y -o \"${output_file}\" -e log -N ${name} bash \"${runner_script}\""
    echo "${qsub_command}"
    eval "${qsub_command}"
    echo ""

    # rm ${runner_script}
}

qc_jobids=""
targetcov_jobids=""

bed=`python -c "import os,sys; print os.path.realpath(sys.argv[1])" ${bed}`

qualimap_bed=${bed}
bed_col_num=`awk '{print NF}' ${bed} | sort -nu | tail -n 1`
if [ ${bed_col_num} -lt 6 ]; then
    tmp_bed="/ngs/tmp/bed_fixed_for_qualimap.bed"
    echo ${tmp_bed}

    if [ ${bed_col_num} -lt 5 ]; then
        awk 'NR==1 {v="\t0\t+"}{print $0v}' ${bed} > ${tmp_bed}
    else
        awk 'NR==1 {v="\t+"}{print $0v}' ${bed} > ${tmp_bed}
    fi
    qualimap_bed=${tmp_bed}
fi


filtered_vcf_suffix=${vcf_suffix}.filtered

for sample in `cat ${samples}`
do
    echo "cd to ${bcbio_final_dir}/${sample}"
    cd "${bcbio_final_dir}/${sample}"
    echo ""

    rm -rf "${sample}${filtered_vcf_suffix}.vcf" annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*

    if [ ! -f "${sample}${vcf_suffix}.vcf" ]; then
        gunzip -c "${sample}${vcf_suffix}.vcf.gz" > "${sample}${vcf_suffix}.vcf"
    fi

    ### InDelFilter ###
    name="indel_filter"
    cmdline="python /group/ngs/bin/InDelFilter.py \"${sample}${vcf_suffix}.vcf\" > \"${sample}${filtered_vcf_suffix}.vcf\""
    run_on_grid "${cmdline}" ${name}_${sample} . "${sample}${filtered_vcf_suffix}.vcf" 1
    indel_filter_name=${name}

    ### VarAnn ###
    name="varannotate"
    mkdir ${name}
    cmdline="python /group/ngs/src/varannotate.py --var \"${sample}${filtered_vcf_suffix}.vcf\" --bam \"${sample}-ready.bam\" -o ${name}"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 1 ${indel_filter_name}_${sample}
    varannotate_name=${name}

    ### VarQC ###
    name="varqc"
    mkdir ${name}
    cmdline="python /group/ngs/src/varqc.py --var \"${sample}${filtered_vcf_suffix}.vcf\" -o ${name}"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 1 ${indel_filter_name}_${sample}
    varqc_name=${name}

    if [ ! -z "${qc_jobids}" ]; then
        qc_jobids=${qc_jobids},${varqc_name}_${sample}
    else
        qc_jobids=${varqc_name}_${sample}
    fi

    ### targetCov ###
    name="targetcov"
    mkdir ${name}
    cmdline="python /group/ngs/src/targetcov.py --bam \"${sample}-ready.bam\" --bed \"${bed}\" --nt=4 -o ${name}"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 4
    targetcov_name=${name}

    if [ ! -z "${targetcov_jobids}" ]; then
        targetcov_jobids=${targetcov_jobids},${targetcov_name}_${sample}
    else
        targetcov_jobids=${targetcov_name}_${sample}
    fi

    ### NGSCat ###
    name="ngscat"
    mkdir ${name}
    cmdline="python /group/ngs/src/ngscat/ngscat.py --bams \"${sample}-ready.bam\" --bed \"${bed}\" --out ${name} --saturation y"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 4

    ## QualiMap ##
    name="qualimap"
    mkdir ${name}
    cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam \"${sample}-ready.bam\" -outdir ${name} -gff \"${qualimap_bed}\" -c -gd HUMAN"
    run_on_grid "${cmdline}" ${name}_${sample} ${name} ${name}/log 8

    cd ..
done

## VarQC summary ##
name="varqc_summary"
cmdline="python /gpfs/group/ngs/src/ngs_reporting/varqc_summary.py ${bcbio_final_dir} ${samples} ${varqc_name} ${filtered_vcf_suffix}"
run_on_grid "${cmdline}" ${name} . log_${name} 1 ${qc_jobids}

## Target coverage summary ##
name="targetcov_summary"
cmdline="python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py ${bcbio_final_dir} ${samples} ${targetcov_name}"
run_on_grid "${cmdline}" ${name} . log_${name} 1 ${targetcov_jobids}
