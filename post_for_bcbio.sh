#!/bin/sh

bed=$1
bcbio_final_dir=$2
samples=$3
vcf_suffix=$4

if [ ! -z "${samples}" ]; then
    samples="./samples.txt"
fi

if [ ! -z "${vcf_suffix}" ]; then
    vcf_suffix="-mutect"
fi

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
        rm ${output_file}
    fi

    runner_script=${output_dir}/runner.sh

    echo ${cmdline}
    echo "#!/bin/bash" > ${runner_script}
    echo "source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;" >> ${runner_script}
    echo ${cmdline} >> ${runner_script}
    chmod +x ${runner_script}

    if [ ! -z "${hold_jid}" ]; then
        hold_jid="-hold_jid "${hold_jid}
    fi

    qsub_command="qsub ${hold_jid} -pe smp ${threads} -S /bin/bash -cwd -q batch.q -b y -o ${output_file} -e log -N ${name} bash ${runner_script}"
    echo "${qsub_command}"
    eval "${qsub_command}"
    echo ""

    rm ${runner_script}
}

qc_jobids=""
targetcov_jobids=""

for sample in `cat ${samples}`
do
    echo "cd to ${bcbio_final_dir}/${sample}"
    cd ${bcbio_final_dir}/${sample}

    rm -rf ${sample}${vcf_suffix}.filtered.vcf annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*

    if [ ! -f ${sample}${vcf_suffix}.vcf ]; then
        gunzip -c ${sample}${vcf_suffix}.vcf.gz > ${sample}${vcf_suffix}.vcf
    fi

    ### InDelFilter ###
    cmdline="python /group/ngs/bin/InDelFilter.py ${sample}${vcf_suffix}.vcf > ${sample}${vcf_suffix}.filtered.vcf"
    run_on_grid "${cmdline}" InDelFilter_${sample} . ${sample}${vcf_suffix}.filtered.vcf 1

    ### VarAnn ###
    mkdir annotation
    cmdline="python /group/ngs/src/varannotate.py --var ${sample}${vcf_suffix}.filtered.vcf --bam "${sample}"-ready.bam -o annotation"
    run_on_grid "${cmdline}" VarAnn_${sample} annotation annotation/log 1 InDelFilter_${sample}

    ### VarQC ###
    mkdir varQC
    cmdline="python /group/ngs/src/varqc.py --var ${sample}${vcf_suffix}.filtered.vcf -o varQC"
    run_on_grid "${cmdline}" VarQC_${sample} varQC varQC/log 1 InDelFilter_${sample}

    if [ ! -z "${qc_jobids}" ]; then
        qc_jobids=${qc_jobids},VarQC_${sample}
    else
        qc_jobids=VarQC_${sample}
    fi

    ### targetCov ###
    mkdir targetSeq
    cmdline="python /group/ngs/src/targetcov.py --bam ${sample}-ready.bam --bed "${bed}" --nt=4 -o targetSeq"
    run_on_grid "${cmdline}" targetSeq_${sample} targetSeq targetSeq/log 4

    if [ ! -z "${targetcov_jobids}" ]; then
        targetcov_jobids=${targetcov_jobids},targetSeq_${sample}
    else
        targetcov_jobids=targetSeq_${sample}
    fi

    ### NGSCat ###
    mkdir NGSCat
    cmdline="python /group/ngs/src/ngscat/ngscat.py --bams ${sample}-ready.bam --bed "${bed}" --out NGSCat --reference /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa --saturation y"
    run_on_grid "${cmdline}" NGSCat_${sample} NGSCat NGSCat/log 4

    ## QualiMap ##
    mkdir QualiMap
    cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam "${sample}"-ready.bam -outdir QualiMap -gff "${bed}" -c -gd HUMAN"
    run_on_grid "${cmdline}" QualiMap_${sample} QualiMap QualiMap/log 8
done

## VarQC summary ##
cmdline="python /gpfs/group/ngs/src/ngs_reporting/varqc_summary.py $bcbio_final_dir $samples varQC"
run_on_grid "${cmdline}" VarQCSummary ${bcbio_final_dir} ../work/log_varqc_summary 1 ${qc_jobids}

## Target coverage summary ##
cmdline="python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py $bcbio_final_dir $samples targetSeq"
run_on_grid "${cmdline}" targetSeqSummary ${bcbio_final_dir} ../work/log_targetcov_summary 1 ${targetcov_jobids}
