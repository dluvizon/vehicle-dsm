#!/bin/bash

echo "Warning! You are going to dowload >40GB of videos. Ready to continue (y/n)?"
read opt
if [ "$opt" != "y" -a "$opt" != "Y" ]
then
  exit
fi

declare -a arr=("01 02" "01 03" "01 04" \
  "02 02" "02 03" "02 04" "02 05" "02 06" "02 07" "02 08" "02 09" "02 10" "02 11" \
  "03 02" \
  "04 02")
unicamp_url="http://www.liv.ic.unicamp.br/~minetto/datasets/download-vehicle-speed"

for s in "${arr[@]}"
do
  sub=${s:0:2}
  vid=${s:3:5}
  subpath="subset${sub}"
  mkdir -p ${subpath}
  pushd ${subpath}

  videopath="${subpath}/video${vid}"
  mkdir -p "video${vid}"
  pushd "video${vid}"

  wget -c "${unicamp_url}/${videopath}/lanes.xml"
  wget -c "${unicamp_url}/${videopath}/matrix.txt"
  wget -c "${unicamp_url}/${videopath}/vehicles.xml"
  wget -c "${unicamp_url}/${videopath}/video.h264"

  popd
  popd
done

