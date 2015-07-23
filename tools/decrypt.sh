#!/bin/bash

date
echo "Processing $1"

gpg_raw=$1
gpg=$(echo $gpg_raw | tr '#' '_');
bam=`echo $gpg | sed 's/\.gpg$//'`

if [ ! -f $bam ]; then
	if [ ! -f $gpg ]; then
		echo "1. Transferring $gpg"
		cmd="aws s3 cp s3://az-ngs-collaborators/Sanger/$gpg_raw $gpg"
		echo $cmd
		eval $cmd
		echo ""
	fi

	if [ ! -f $gpg ]; then
		>&2 echo "$gpg does not exist"
		exit 1
	fi

	echo "2. Decrypting to $bam"
	cmd="gpg --batch --no-use-agent --passphrase-file /ngs/sanger_cell_line/passphrase.txt -o $bam --decrypt $gpg"
	echo $cmd
	eval $cmd
	echo ""
fi

if [ ! -f $bam ]; then
	>&2 echo "$bam does not exist"
	exit 1
fi

if [ -f $gpg ]; then
	echo "rm $gpg"
	rm $gpg
fi

r=`echo $bam | sed 's/\.bam$//'`
r1=$r.1.fastq
r2=$r.2.fastq
echo "3. Exracting fastq files $r1 and $r2"
cmd="/usr/local/share/bcbio/bin/bam2fastx --fastq --all -o $r.fastq --paired $bam"
echo $cmd
eval $cmd
echo ""

if [ ! -f $r1 ]; then
	>&2 echo "$r1 does not exist"
	exit 1
fi

if [ ! -f $r2 ]; then
	>&2 echo "$r2 does not exist"
	exit 1
fi

echo "rm $bam"
rm $bam

echo "gzip $r1"
gzip $r1
echo "gzip $r2"
gzip $r2

echo "4. Tranferring fastq files $r1.gz and $r3.gz back to s3"
cmd="aws s3 cp $r1.gz s3://az-ngs-collaborators/sanger-decrypted --sse --profile klpf990"
echo $cmd
eval $cmd
echo ""
cmd="aws s3 cp $r2.gz s3://az-ngs-collaborators/sanger-decrypted --sse --profile klpf990"
echo $cmd
eval $cmd
echo ""

echo "Done"
date


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

/gpfs/group/ngs/src/az.reporting-0.3/scripts/post_bcbio/varfilter.py --sys-cnf /gpfs/group/ngs/src/az.reporting-0.3/system_info_Waltham.yaml --run-cnf /ngs/oncology/analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/config/run_info.yaml --project-name Bio_0038_Bio_0038_150521_D00443_0159_AHK2KTADXX -t 10 --log-dir /ngs/oncology/analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/final/2015-05-29_bcbio/log/varFilter-ensemble /ngs/oncology/analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/final --caller ensemble --done-marker /ngs/oncology/analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/work/post_processing/done_markers/varFilter_-ZNIHQ=_Bio_0038_Bio_0038_150521_D00443_0159_AHK2KTADXX_ensemble.done --reuse

/gpfs/group/ngs/src/bcbio-nextgen/0.8.7/rhel6-x64/anaconda/bin/python2.7  /gpfs/group/ngs/src/az.reporting-0.3/sub_scripts/varannotate.py  --sys-cnf /gpfs/group/ngs/src/az.reporting-0.3/system_info_Waltham.yaml --run-cnf /gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/config/run_info.yaml -t 1 --log-dir /gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/final/2015-05-29_bcbio/log/varAnnotate --genome hg19 --project-name Bio_0038_150521_D00443_0159_AHK2KTADXX_bcbio  --vcf '/gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/final/AGS/var/AGS-ensemble.vcf.gz' --bam /gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/final/AGS/AGS-ready.bam  -o '/gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/final/AGS/varAnnotate' -s 'AGS' -c ensemble --work-dir '/gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/work/post_processing/varAnnotate_AGS_ensemble'  --done-marker /gpfs/ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio/work/post_processing/done_markers/varAnnotate_Ct6ycg=_Bio_0038_150521_D00443_0159_AHK2KTADXX_bcbio_AGS_ensemble.done

you want me to remove already downloaded files?
gpg --batch --passphrase-file /nfs/glusterfs/decrypt_gpg/passphrase.txt -o $(echo 6106_4#2.bam.gpg | sed 's/\.gpg$//') --decrypt 6106_4#2.bam.gpg
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


gpg: can't open `/opt/sge6/default/spool/exec_spool_local/bcbio-starcluster-Ronent-TestGluster-D2-030-node001/job_scripts/passphrase.txt': No such file or directory
6142_3#2.bam does not exists
gpg: CAST5 encrypted data
gpg: gpg-agent is not available in this session
gpg: can't query passphrase in batch mode
gpg: encrypted with 1 passphrase
gpg: decryption failed: bad key
6142_3_2.bam does not exists

