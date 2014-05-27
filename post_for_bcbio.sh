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
    output_file=$3
    hold_jid=$4
    # if [ -z "$3" ]; then
    #         output="log"
    # fi
    # rm log
    rm ${output_file}
    echo ${cmdline}
    echo "#!/bin/bash" > cmd.sh
    echo "source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load be${cmdline}e}module load samtools;" >> cmd.sh
    echo ${cmdline} >> cmd.sh
    chmod +x cmd.sh

    if [ -z "${hold_jid}" ]; then
        hold_jid=" -hold_jid "${hold_jid}
    fi

    qsub ${hold_jid} -pe smp 8 -S /bin/bash -cwd -q batch.q -b y -o ${output_file} -e log -N ${name} bash ./cmd.sh

    rm cmd.sh
}

for sample in `cat ${samples}`
do
    echo "cd to ${finalDir}/${sample}"
    cd $finalDir/$sample

    rm -rf cmd.sh log ${sample}${vcf_suffix}.filtered.vcf annotation varQC targetSeq NGSCat QualiMap *tmp* work *ready_stats*

    if [ ! -f ${sample}${vcf_suffix}.vcf ];
    then
        gunzip -c ${sample}${vcf_suffix}.vcf.gz > ${sample}${vcf_suffix}.vcf
    fi

    ### InDelFilter ###
    cmdline="python /group/ngs/bin/InDelFilter.py ${sample}${vcf_suffix}.vcf > ${sample}${vcf_suffix}.filtered.vcf"
    run_on_grid "${cmdline}" InDelFilter_${sample} ${sample}${vcf_suffix}.filtered.vcf

    ### VarAnn ###
    mkdir annotation
    cmdline="python /group/ngs/src/varannotate.py --var ../${sample}${vcf_suffix}.filtered.vcf --bam ../"${sample}"-ready.bam --nt=4 -o annotation"
    run_on_grid "${cmdline}" VarAnn_${sample} annotation/log InDelFilter_${sample}

    ### VarQC ###
    mkdir varQC
    cmdline="python /group/ngs/src/varqc.py --var ../${sample}${vcf_suffix}.filtered.vcf --nt=4 -o varQC"
    run_on_grid "${cmdline}" VarQC_${sample} varQC/log InDelFilter_${sample}

    ### targetCov ###
    mkdir targetSeq
    cmdline="python /group/ngs/src/targetcov.py --bam ../${sample}-ready.bam --bed "${bed}" --nt=4 -o targetSeq"
    run_on_grid "${cmdline}" targetSeq_${sample} targetSeq/log

    ### NGSCat ###
    mdkir NGSCat
    cmdline="python /group/ngs/src/ngscat/ngscat.py --bams ../${sample}-ready.bam --bed "${bed}" --out NGSCat --reference /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa --saturation y"
    run_on_grid "${cmdline}" NGSCat_${sample} NGSCat/log

    ## QualiMap ##
    mkdir QualiMap
    cmdline="/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam ../"${sample}"-ready.bam -outdir QualiMap -gff "${bed}" -c -gd HUMAN"
    run_on_grid "${cmdline}" QualiMap_${sample} QualiMap/log
done

qc_summary=qc_summary_report.txt