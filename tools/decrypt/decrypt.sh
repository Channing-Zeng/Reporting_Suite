#!/bin/bash
here=`dirname $0`
while read gpg_raw
do
	gpg=$(echo $gpg_raw | tr '#' '_');
	cmdl="qsub -cwd -pe smp 1 -N decrypt_$gpg -q all.q -j n -o log/decrypt_$gpg.log -e log/decrypt_$gpg.log $here/decrypt_one.sh $gpg_raw"
	# cmdl="qsub -cwd -pe smp 1 -N decrypt_2_$gpg -hold_jid decrypt_$gpg -q all.q -j n -o log/decrypt_$gpg_2.log -e log/decrypt_$gpg_2.err $here/decrypt.sh $gpg"
	echo $cmdl
	eval $cmdl
done <$1


# #!/bin/bash
# here=`dirname $0`
# for gpg in *.bam.gpg
# do
# 	bam=`echo $gpg | sed 's/\.gpg$//'`
# 	gpg_cmdl="gpg --batch --passphrase-file $here/passphrase.txt -o $bam --decrypt $gpg"
# 	gpg_job_name=$(basename $gpg);
# 	qsub_cmdl="qsub -cwd -pe smp 1 -q all.q -j n -N $gpg_job_name $here/runner.sh \"$gpg_cmdl\""
# 	echo $qsub_cmdl
# 	eval $qsub_cmdl

# 	r=`echo $bam | sed 's/\.bam$//'`
# 	r1=$r.1.fastq
# 	r2=$r.2.fastq

# 	bam2fastx_cmdl="bam2fastx --fastq --all -o $r.fastq --paired $bam; gzip $r1; gzip $r2"
# 	qsub_cmdl="qsub -cwd -pe smp 1 -q all.q -j n -hold_jid $gpg_job_name $here/runner.sh \"$bam2fastx_cmdl\""
# 	echo $qsub_cmdl
# 	eval $qsub_cmdl
# done
