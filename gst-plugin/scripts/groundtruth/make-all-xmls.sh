#!/bin/bash

INP=metrologico.log
OUT=ground-truth.xml

function buildxml() {
	VPATH=../../../../videos/video"$1"
	./build-xml.sh $VPATH/$INP $VPATH/$OUT $2
}

#     video	delay
buildxml 01	210
buildxml 02	164
buildxml 03	410
buildxml 04	-17
buildxml 05	105
buildxml 06	-80
buildxml 07	-38
buildxml 08	 20
buildxml 09	165
buildxml 10	-74
buildxml 11	-68
buildxml 12	-52
buildxml 13	255
buildxml 14	 84
buildxml 15	-65
buildxml 16	110
