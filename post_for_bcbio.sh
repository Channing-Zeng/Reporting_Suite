#!/bin/sh

bed=$1
bcbio_final_dir=$2
samples=$3
vcf_suffix=$4

if [ -z "${samples}" ]; then
    samples="./samples.txt"
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
    echo "source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;" >> "${runner_script}"
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


bed_col_num=`awk '{print NF}' ${bed} | sort -nu | tail -n 1`
if [ ${bed_col_num} -lt 6 ]; then
    tmp_bed="/ngs/tmp/bed_fixed_for_qualimap.bed"
    echo ${tmp_bed}

    if [ ${bed_col_num} -lt 5 ]; then
        awk 'NR==1 {v="\t0\t+"}{print $0v}' ${bed} > ${tmp_bed}
    else
        awk 'NR==1 {v="\t+"}{print $0v}' ${bed} > ${tmp_bed}
    fi
    qualimap_cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam \"${sample}-ready.bam\" -outdir QualiMap -gff \"${tmp_bed}\" -c -gd HUMAN"
else
    qualimap_cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam \"${sample}-ready.bam\" -outdir QualiMap -gff \"${bed}\" -c -gd HUMAN"
fi


filtered_vcf_suffix=${vcf_suffix}.filtered

for sample in `cat ${samples}`
do
    echo "cd to ${bcbio_final_dir}/${sample}"
    cd "${bcbio_final_dir}/${sample}"
    echo ""

#    rm -rf "${sample}${filtered_vcf_suffix}.vcf" annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*
#
#    if [ ! -f "${sample}${vcf_suffix}.vcf" ]; then
#        gunzip -c "${sample}${vcf_suffix}.vcf.gz" > "${sample}${vcf_suffix}.vcf"
#    fi
#
#    ### InDelFilter ###
#    cmdline="python /group/ngs/bin/InDelFilter.py \"${sample}${vcf_suffix}.vcf\" > \"${sample}${filtered_vcf_suffix}.vcf\""
#    run_on_grid "${cmdline}" InDelFilter_${sample} . "${sample}${filtered_vcf_suffix}.vcf" 1
#
#    ### VarAnn ###
#    mkdir annotation
#    cmdline="python /group/ngs/src/varannotate.py --var \"${sample}${filtered_vcf_suffix}.vcf\" --bam \"${sample}-ready.bam\" -o annotation"
#    run_on_grid "${cmdline}" VarAnn_${sample} annotation annotation/log 1 InDelFilter_${sample}
#
#    ### VarQC ###
#    mkdir varQC
#    cmdline="python /group/ngs/src/varqc.py --var \"${sample}${filtered_vcf_suffix}.vcf\" -o varQC"
#    run_on_grid "${cmdline}" VarQC_${sample} varQC varQC/log 1 InDelFilter_${sample}
#
#    if [ ! -z "${qc_jobids}" ]; then
#        qc_jobids=${qc_jobids},VarQC_${sample}
#    else
#        qc_jobids=VarQC_${sample}
#    fi
#
#    ### targetCov ###
#    mkdir targetSeq
#    cmdline="python /group/ngs/src/targetcov.py --bam \"${sample}-ready.bam\" --bed \"${bed}\" --nt=4 -o targetSeq"
#    run_on_grid "${cmdline}" targetSeq_${sample} targetSeq targetSeq/log 4
#
#    if [ ! -z "${targetcov_jobids}" ]; then
#        targetcov_jobids=${targetcov_jobids},targetSeq_${sample}
#    else
#        targetcov_jobids=targetSeq_${sample}
#    fi
#
#    ### NGSCat ###
#    mkdir NGSCat
#    cmdline="python /group/ngs/src/ngscat/ngscat.py --bams \"${sample}-ready.bam\" --bed \"${bed}\" --out NGSCat --reference /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa --saturation y"
#    run_on_grid "${cmdline}" NGSCat_${sample} NGSCat NGSCat/log 4

    ## QualiMap ##
    mkdir QualiMap
    run_on_grid "${qualimap_cmdline}" QualiMap_${sample} QualiMap QualiMap/log 8

    cd ..
done
#
### VarQC summary ##
#cmdline="python /gpfs/group/ngs/src/ngs_reporting/varqc_summary.py $bcbio_final_dir $samples varQC ${filtered_vcf_suffix}"
#run_on_grid "${cmdline}" VarQCSummary . ../work/log_varqc_summary 1 ${qc_jobids}
#
### Target coverage summary ##
#cmdline="python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py $bcbio_final_dir $samples targetSeq"
#run_on_grid "${cmdline}" targetSeqSummary . ../work/log_targetcov_summary 1 ${targetcov_jobids}
