#!/bin/bash

date
echo "Processing $1"

gpg_raw=$1
gpg=$(echo $gpg_raw | tr '#' '_');
bam=`echo $gpg | sed 's/\.gpg$//'`
r=`echo $bam | sed 's/\.bam$//'`
r1=${r}.1.fastq
r2=${r}.2.fastq
zr1=${r1}.gz
zr2=${r2}.gz

if [ ! -f $zr1 ] || [ ! -f $zr2 ]; then
	if [ ! -f $bam ]; then
		if [ ! -f $gpg ]; then
			date
			echo "1. Transferring $gpg"
			cmd="aws s3 cp s3://az-ngs-collaborators/Sanger/$gpg_raw $gpg"
			echo $cmd
			eval $cmd
			echo ""
		fi

		if [ ! -f $gpg ]; then
			date
			>&2 echo "$gpg does not exist"
			exit 1
		fi

		date
		echo "2. Decrypting to $bam"
		cmd="gpg --batch --no-use-agent --passphrase-file /ngs/sanger_cell_line/passphrase.txt -o $bam --decrypt $gpg"
		echo $cmd
		eval $cmd
		echo ""
	fi

	if [ ! -f $bam ]; then
		date
		>&2 echo "$bam does not exist"
		exit 1
	fi

	if [ -f $gpg ]; then
		date
		echo "rm $gpg"
		rm $gpg
	fi

	date
	echo "3. Exracting fastq files $r1 and $r2"
	cmd="/usr/local/share/bcbio/bin/bam2fastx --fastq --all -o $r.fastq --paired $bam"
	echo $cmd
	eval $cmd
	echo ""

	if [ ! -f $r1 ]; then
		date
		>&2 echo "$r1 does not exist"
		exit 1
	fi

	if [ ! -f $r2 ]; then
		date
		>&2 echo "$r2 does not exist"
		exit 1
	fi

	if [ -f $bam ]; then
		date
		echo "rm $bam"
		rm $bam
	fi

	date
	echo "gzip $r1"
	gzip $r1

	date
	echo "gzip $r2"
	gzip $r2
fi

if [ ! -f $zr1 ]; then
	date
	>&2 echo "$zr1 does not exist"
	exit 1
fi

if [ ! -f $zr2 ]; then
	date
	>&2 echo "$zr2 does not exist"
	exit 1
fi

date
echo "4. Tranferring fastq files $zr1 and $zr2 back to s3"
cmd="aws s3 cp $zr1 s3://az-ngs-collaborators/sanger-decrypted/ --sse --profile klpf990"
echo $cmd
eval $cmd
ret_code=$?
echo ""
if [ ${ret_code} == 0 ] && [ -f $zr1 ]; then
	echo "rm $zr1"
	rm $zr1
fi

date
cmd="aws s3 cp $zr2 s3://az-ngs-collaborators/sanger-decrypted/ --sse --profile klpf990"
echo $cmd
eval $cmd
ret_code=$?
echo ""
if [ ${ret_code} == 0 ] && [ -f $zr2 ]; then
	echo "rm $zr2"
	rm $zr2
fi

echo "Done"
date

