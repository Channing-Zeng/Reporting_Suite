#!/bin/sh

bed=$1
finalDir=$2
samples=$3
vcf_suffix=$4

if [ -z "${samples}" ]; then
    samples="./samples.txt"
fi

if [ -z "${vcf_suffix}" ]; then
    vcf_suffix="-mutect"
fi

function run_on_grid {
    cmdline=$1
    name=$2
    output_dir=$3
    output_file=$4
    hold_jid=$5
    # if [ -z "$3" ]; then
    #         output="log"
    # fi
    # rm log
    if [ -f "${output_file}" ]; then
        rm ${output_file}
    fi

    runner_script=${output_dir}/runner.sh

    echo ${cmdline}
    echo "#!/bin/bash" > ${runner_script}
    echo "source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;" >> ${runner_script}
    echo ${cmdline} >> ${runner_script}
    chmod +x ${runner_script}

    if [ ! -z "${hold_jid}" ]; then
        hold_jid=" -hold_jid "${hold_jid}
    fi

    qsub_command="qsub ${hold_jid} -pe smp 8 -S /bin/bash -cwd -q batch.q -b y -o ${output_file} -e log -N ${name} bash ${runner_script}"
    echo "${qsub_command}"
    eval "${qsub_command}"
    echo ""
}

qc_jobids=""
targetcov_jobids=""

for sample in `cat ${samples}`
do
    echo "cd to ${finalDir}/${sample}"
    cd $finalDir/$sample

    rm -rf ${sample}${vcf_suffix}.filtered.vcf annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*

    if [ ! -f ${sample}${vcf_suffix}.vcf ]; then
        gunzip -c ${sample}${vcf_suffix}.vcf.gz > ${sample}${vcf_suffix}.vcf
    fi

    ### InDelFilter ###
    cmdline="python /group/ngs/bin/InDelFilter.py ${sample}${vcf_suffix}.vcf > ${sample}${vcf_suffix}.filtered.vcf"
    run_on_grid "${cmdline}" InDelFilter_${sample} . ${sample}${vcf_suffix}.filtered.vcf

    ### VarAnn ###
    mkdir annotation
    cmdline="python /group/ngs/src/varannotate.py --var ../${sample}${vcf_suffix}.filtered.vcf --bam ../"${sample}"-ready.bam --nt=4 -o annotation"
    run_on_grid "${cmdline}" VarAnn_${sample} annotation annotation/log InDelFilter_${sample}

    ### VarQC ###
    mkdir varQC
    cmdline="python /group/ngs/src/varqc.py --var ../${sample}${vcf_suffix}.filtered.vcf --nt=4 -o varQC"
    run_on_grid "${cmdline}" VarQC_${sample} varQC varQC/log InDelFilter_${sample}

    if [ ! -z "${qc_jobids}" ]; then
        qc_jobids=${qc_jobids},VarQC_${sample}
    else
        qc_jobids=VarQC_${sample}
    fi

    ### targetCov ###
    mkdir targetSeq
    cmdline="python /group/ngs/src/targetcov.py --bam ../${sample}-ready.bam --bed "${bed}" --nt=4 -o targetSeq"
    run_on_grid "${cmdline}" targetSeq_${sample} targetSeq targetSeq/log

    if [ ! -z "${targetcov_jobids}" ]; then
        targetcov_jobids=${targetcov_jobids},targetSeq_${sample}
    else
        targetcov_jobids=targetSeq_${sample}
    fi

    ### NGSCat ###
    mkdir NGSCat
    cmdline="python /group/ngs/src/ngscat/ngscat.py --bams ../${sample}-ready.bam --bed "${bed}" --out NGSCat --reference /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa --saturation y"
    run_on_grid "${cmdline}" NGSCat_${sample} NGSCat NGSCat/log

    ## QualiMap ##
    mkdir QualiMap
    cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam ../"${sample}"-ready.bam -outdir QualiMap -gff "${bed}" -c -gd HUMAN"
    run_on_grid "${cmdline}" QualiMap_${sample} QualiMap QualiMap/log
done

cmdline="python varqc_summary.py $finalDir $samples varQC"
run_on_grid "${cmdline}" VarQCSummary ${finalDir} log ${qc_jobids}

cmdline="python targetcov_summary.py $finalDir $samples targetSeq"
run_on_grid "${cmdline}" targetSeqSummary ${finalDir} log ${targetcov_jobids}
