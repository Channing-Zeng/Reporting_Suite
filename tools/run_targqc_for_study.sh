dir=$(readlink $1)
project="$(basename($(dirname $dir))_$(basename $dir)"

FLAGS=$2

TARGQC="python /gpfs/group/ngs/src/az.reporting/targqc.py"

$TARGQC $dir/*.bam --bed $dir/../../bed/UCSC_canonical_exons.key.bed -o $dir/tq_capture $FLAGS
if [ -f $dir/../capture.bed ];
then
	$TARGQC $dir/*.bam --bed $dir/../capture.bed -o $dir/tq_keygenes $FLAGS
fi

