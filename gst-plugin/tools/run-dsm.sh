#!/bin/bash

if [ $# != 1 ]
then
	echo "Usage: $0 <video path> # e.g. ../../../videos/video01"
	exit 0
fi

#-----------------------------------------------------------------------------
# OPTIONS:
#         fast:    run without video or file output
#         display: run with realtime video output
#         store:   run and generate one jpeg for each frame
#-----------------------------------------------------------------------------
opt=display

REBUILD="no"

# remove trailing slash from argument
UARG="$1"
il=$((${#UARG}-1))
if [ "${UARG:$il:1}" = "/" ]
then
	UARG=${UARG:0:$il}
fi
VIDEOPATH=$UARG

# change this to 'false' to display gstreamer prints
SILENT=false
IWIDTH=1920
IHEIGHT=1080
GTPATH_CMD="groundtruth-path=$VIDEOPATH/ground-truth.xml"
VEHICLES_XML="vehicles-xml=$VIDEOPATH/vehicles.xml"
OPATH_CMD="output-path=${VIDEOPATH}/output"
LIBS="../gstreamer/.libs/"

#-----------------------------------------------------------------------------
show_ok_here()
{
	echo "`pwd`: Ok here..."
}

check_training()
{
	pushd ../../training/icdar/
	if [ ! -e .run_training.stamp ]
	then
		touch .run_training.stamp
		./run.sh
	fi
	show_ok_here
	popd
}

prepare_svm_model()
{
	POS=../../training/icdar/positives.txt
	NEG=../../training/icdar/negatives.txt
	OUT=./saida

	rm -rf ${OUT}/positives.txt ${OUT}/positives.txt ./input/model.svm

	if [ ! -x train/train ] || [ ! -x train/svm-train ]
	then
		echo "Binary train/train or train/svm-train not found"
		exit 1
	fi

	./train/train ${POS} ${OUT}/positives.txt  1 ./input/1_7_9.txt
	./train/train ${NEG} ${OUT}/negatives.txt -1 ./input/1_7_9.txt

	cat ${OUT}/positives.txt ${OUT}/negatives.txt > ${OUT}/ALL.txt

	./train/svm-train -s 0 -c 1000 -t 2 -g 1 -r 1 -d 3 -b 1 \
			${OUT}/ALL.txt ./input/model.svm
}

check_svm_model()
{
	if [ ! -e input/model.svm ]
	then
		echo "File 'input/model.svm' not found!"
		#echo "Trying to run ./train/train and ./train/svm-train steps."
		#prepare_svm_model
		exit 1
	fi
}

clear_crete_dir()
{
	rm -rf "$1"
	mkdir -p "$1"
}

check_output_path()
{
	if [ -d "$VIDEOPATH" ]
	then
		clear_crete_dir "${VIDEOPATH}/output"
		mkdir "${VIDEOPATH}/output/snapshot"
		mkdir -p "${VIDEOPATH}/output/groundtruth/f1"
		mkdir -p "${VIDEOPATH}/output/groundtruth/f2"
		mkdir -p "${VIDEOPATH}/output/groundtruth/f3"
		mkdir -p "${VIDEOPATH}/output/detection/f1"
		mkdir -p "${VIDEOPATH}/output/detection/f2"
		mkdir -p "${VIDEOPATH}/output/detection/f3"
		clear_crete_dir "/tmp/jpeg"
		clear_crete_dir "/tmp/slopes"
		if [ ! -d saida ]
		then
			mkdir -p saida
		fi
	fi
}


check_before_run()
{
	if [ ! -e run-dsm.sh ]
	then
		echo "Please, run this script directly from gst-plugin/tools"
		exit 1
	fi
	# check_training
	check_svm_model
	check_output_path
}

func_main()
{
	check_before_run

	if [ $REBUILD = "yes" ]
	then
		echo "Rebuilding software..."
		if [ -x ./bake.sh ]
		then
			MAKE="./tools/bake.sh"
		else
			MAKE="make"
		fi
		pushd ..
			$($MAKE)
			bRET=$?
		popd
		if [ $bRET != 0 ]
		then
			exit 1
		fi
	fi

	case $opt in
	fast)
		gst-launch-1.0 --gst-plugin-path=${LIBS} \
			filesrc location="${VIDEOPATH}/video.h264" ! \
			decodebin ! videoconvert ! \
			video/x-raw,format=\(string\)NV12,\
			width=$IWIDTH,height=$IHEIGHT ! \
			vehicledsm silent=$SILENT $GTPATH_CMD $OPATH_CMD \
			$VEHICLES_XML ! fakesink
		;;
	display)
		gst-launch-1.0 --gst-plugin-path=${LIBS} \
			filesrc location="${VIDEOPATH}/video.h264" ! \
			decodebin ! videoconvert ! \
			video/x-raw,format=\(string\)NV12,\
			width=$IWIDTH,height=$IHEIGHT ! \
			vehicledsm silent=$SILENT $GTPATH_CMD $OPATH_CMD \
			$VEHICLES_XML ! videoconvert ! videoscale qos=false ! \
			ximagesink sync=false
		;;
	store)
		mkdir -p "${VIDEOPATH}/output/jpeg"
		gst-launch-1.0 --gst-plugin-path=${LIBS} \
			filesrc location="${VIDEOPATH}/video.h264" ! \
			decodebin ! videoconvert ! \
			video/x-raw,format=\(string\)NV12,\
			width=$IWIDTH,height=$IHEIGHT ! \
			vehicledsm silent=$SILENT $GTPATH_CMD $OPATH_CMD \
			$VEHICLES_XML ! videoconvert ! \
			jpegenc quality=70 idct-method=2 ! multifilesink \
			location="/tmp/jpeg/%05d.jpeg"
		;;
	*)
		echo "Bad option '${opt}'."
		echo "Valid options are: [fake] [display] [store]"
		exit 1
		;;
	esac
}

func_main
