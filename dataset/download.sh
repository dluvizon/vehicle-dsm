#!/bin/bash

echo "Warning! You are going to dowload 7GB of videos. Ready to continue (y/n)?"
read opt
if [ "$opt" != "y" -a "$opt" != "Y" ]
then
  exit
fi

for s in 1 2 3 4 5
do
  setdir="set${s}"
  mkdir ${setdir}
  wget http://www.dainf.ct.utfpr.edu.br/~rminetto/projects/vehicle-speed/dataset/Set0${s}_video01.h264
  mv Set0${s}_video01.h264 ${setdir}/video.h264
  wget http://www.dainf.ct.utfpr.edu.br/~rminetto/projects/vehicle-speed/dataset/Set0${s}_video01.xml
  mv Set01_video0${s}.xml ${setdir}/vehicles.xml
done

