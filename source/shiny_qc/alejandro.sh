#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --exclusive
. /etc/bashrc
 

#tmp=/scratch/$LOGNAME.$RANDOM/
#mkdir $tmp
src=$1  #/ngs/projects/DS/fromNextSeq500/$1
res=/ngs/projects/DS/PreProcessed/$(basename $src)

#######################################################################
#################### RUN BCL2FASTQ CONVERSION #########################
#######################################################################
cpu=$(grep -c processor /proc/cpuinfo)
module add bcl2fastq
src=Data/Intensities/BaseCalls
sheet=$src/SampleSheet.csv
bcl2fastq -r $cpu -d $cpu -p $cpu -w $cpu -R $src  --no-lane-splitting -o $tmp --sample-sheet $sheet --fastq-compression-level 1 --interop-dir $tmp/InterOp
 
src="should not access src from this point onward"
 
#######################################################################
#################### RUN FASTQC REPORTS ###############################
#######################################################################
module add fastqc
fqs=$(find $tmp | grep fastq.gz)
mkdir $tmp/FastQC
fastqc --extract -o $tmp/FastQC -t $cpu $fqs
 


#######################################################################
################### METAMAPPING (CONTAMINATION) #######################
#######################################################################
   
    mkdir $tmp/MetaMapping
 
    fqs=$(find $tmp | grep fastq.gz | grep _R1_001.fastq.gz)
    i=1
    for x in $fqs
    do      
       let i+=1
       echo $x
       s=$(basename $x | sed s/_R1_001.fastq.gz//g)
       pigz -d -c $x \
        | head -n 1000000 \
        | awk -v s=$s '{print (NR%4 == 1) ? "@"s : $0}' > $tmp/MetaMapping/tmp.$i &
    done
    wait
  cat $tmp/MetaMapping/tmp.* > $tmp/MetaMapping/Sample.fq
 
 
    module add bwa
for ref in $(cat /ngs/projects/DS/Scripts/metaScreen.lst)
do
       refName=$(basename $ref | sed s/.fa//g)
       echo $refName $x
       bwa bwasw -t $cpu $ref $tmp/MetaMapping/Sample.fq | awk -v y=$x -v z=$refName '{
        if ($1 != "@SQ"){
         MQ="x";
         if ( $5 < 30 )         
            MQ="AMB";  
         if ( $5 > 0 ) 
            MQ="LOW";
         if ( $5 > 30 )
            MQ="HIGH";
         if ( $2 == 4 ) 
            MQ="UNK";
         MQ=$1"\t"$3"\t"MQ;
         freq[MQ]++
        }
       }
       END {
           for (word in freq)
             if (freq[word] > 5)
               printf "%s\t%s\t%s\t%d\n", y, z, word, freq[word]
       }' >> $tmp/MetaMapping/Mapping1M.tab
done
 
rm $tmp/MetaMapping/Sample.fq
rm $tmp/MetaMapping/tmp.*
 
mv $tmp $res