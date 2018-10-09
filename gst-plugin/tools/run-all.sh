#!/bin/bash
#-----------------------------------------------------------------------------

if [ $# != 1 ]
then
	echo "Usage: $0 <video root path, with file LIST>"
	exit 0
fi

VIDEOROOT="$1"

if [ ! -x run-dsm.sh ]
then
	echo "File 'run-dsm.sh' not found!"
	exit 1
fi

if [ ! -e "$VIDEOROOT"/LIST ]
then
	echo "File '$VIDEOROOT/LIST' not found!"
	exit 1
fi
vlist=(`cat "$VIDEOROOT"/LIST`)
absvlist=()
for i in ${vlist[*]}
do
	VIDEOPATH="${VIDEOROOT}/${i}"
	absvlist=(${absvlist[*]} "$VIDEOPATH")
done

for i in ${vlist[*]}
do
	VIDEOPATH="${VIDEOROOT}/${i}"
	echo "Running for '$VIDEOPATH' ..."
	#./run-dsm.sh $VIDEOPATH | tee $VIDEOPATH/log.txt
	#./analyzeLog $VIDEOPATH >> $VIDEOPATH/log.txt
done
./analyzeLog ${absvlist[*]}

exit 0

