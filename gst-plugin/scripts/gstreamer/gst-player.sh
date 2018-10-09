#!/bin/bash

if [ $# != 1 ]
then
	echo "Usage: $0 <video .avi | .h264 | .mp4>"
	exit 1
fi

UARG="$1"
IWIDTH=1920
IHEIGHT=1080

run_fast="fakesink"
run_store="videoconvert ! theoraenc ! oggmux ! fakesink"
run_display="videoconvert ! videoscale qos=false ! ximagesink sync=false"

gst-launch-1.0 filesrc location="$UARG" ! decodebin ! videoconvert ! \
	video/x-raw,format=\(string\)NV12,width=$IWIDTH,height=$IHEIGHT ! \
	${run_display}

