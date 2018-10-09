#!/bin/bash

if [ $# != 3 ]
then
	echo -e "Usage: $0"
	echo -e "\t\t<input file , e.g. radar.log>"
	echo -e "\t\t<output file, e.g. ground-truth.xml>"
	echo -e "\t\t<delay in frames>"
	exit 1
fi

TMP=/tmp/build-xml
FRAME_DELAY=$3
FIRST_FRAC=
VIDEO_MS_MULT=199
VIDEO_MS_DIV=6
OUTPUT_FILE=$2

mkdir -p $TMP
cat $1 | grep 'VEL_COR\|LOOP1_OUT' > $TMP/met.log

function check_mili {
	val=$1
	val=${val#0}
	sz=${#val}
	if [ "$sz" == "3" ]
	then
		echo $((10#$val))
	elif [ "$sz" == "2" ]
	then
		echo $((10#$val * 10))
	elif [ "$sz" == "1" ]
	then
		echo $((10#$val * 100))
	fi
}

rm -f $OUTPUT_FILE
function write_xml_line {
	echo "$1" >> $OUTPUT_FILE
}

TMS_I=$((VIDEO_MS_MULT / VIDEO_MS_DIV))
TMS_F=$(((VIDEO_MS_MULT * 1000) / VIDEO_MS_DIV - (TMS_I * 1000)))
STR_PERIOD_MS="$TMS_I.$TMS_F"

write_xml_line '<?xml version="1.0" encoding="UTF-8"?>'
write_xml_line '<tagset>'
write_xml_line "  <!-- Source file = $1 -->"
write_xml_line "  <!-- Frame delay = $FRAME_DELAY -->"
write_xml_line "  <!-- Frame interval [ms] = $STR_PERIOD_MS -->"
write_xml_line '  <GroundTruth>'

LAST_ID=("" "" "")
LAST_SPEED=("" "" "")
LAST_FSHOW=("" "" "")

while read line
do
	time=`echo "$line" | awk '{print $2}'`
	hou=`echo $time | awk -F '[/:]' '{print $1}'`
	hou=${hou#0}
	min=`echo $time | awk -F '[/:]' '{print $2}'`
	min=${min#0}
	sec=`echo $time | awk -F '[/:.]' '{print $3}'`
	sec=${sec#0}
	mil=`echo $time | awk -F '[/:.]' '{print $4}'`
	mil=`check_mili $mil`
	#echo "$time : $hou $min $sec $mil"
	frac=$(((((((hou * 60) + min) * 60) + sec) * 1000) + mil))
	if [ "$FIRST_FRAC" == "" ]
	then
		FIRST_FRAC=$frac
	fi
	frac=$((frac - FIRST_FRAC))
	frame=$((((frac * VIDEO_MS_DIV) / VIDEO_MS_MULT) + FRAME_DELAY))
	fshow=$((frame - (VIDEO_MS_MULT / (2 * VIDEO_MS_DIV))))
	fend=$((frame + (VIDEO_MS_MULT / VIDEO_MS_DIV)))

	info=`echo "$line" | awk '{print $3}'`
	lane=${info:3:1}
	id=${info:5}

	msg=`echo "$line" | awk '{print $4}'`
	aux=`echo $msg | tr "=" "\n"`

	case $lane in
		1 )
			ilane=0 ;;
		2 )
			ilane=1 ;;
		3 )
			ilane=2 ;;
		* )
			continue ;;
	esac

	if [ "${LAST_ID[$ilane]}" != "$id" ] && [ "${LAST_ID[$ilane]}" != "" ] && [ "${LAST_SPEED[$ilane]}" != "-1" ]
	then
		# than write the last value
		write_xml_line "    <vehicle frame=\"${LAST_FSHOW[$ilane]}\" lane=\"$lane\" speed=\"${LAST_SPEED[$ilane]}\" />"
	fi

	for ai in $aux
	do
		if [ "$ai" == "VEL_COR" ]
		then
			speed=`echo $msg | awk -F '[/=]' '{print $2}'`
			write_xml_line "    <vehicle frame=\"$fshow\" lane=\"$lane\" speed=\"$speed\" />"
			LAST_SPEED[$ilane]=-1
		elif [ "$ai" == "LOOP1_OUT" ]
		then
			LAST_SPEED[$ilane]=-666
			LAST_FSHOW[$ilane]=$fshow
		fi
	done
	LAST_ID[$ilane]=$id
done < $TMP/met.log

write_xml_line '  </GroundTruth>'
write_xml_line '</tagset>'
rm -rf $TMP

