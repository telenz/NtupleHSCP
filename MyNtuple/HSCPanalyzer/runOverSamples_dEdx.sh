#!/bin/bash

for m in {100,200,300,400,500,600,700,800}
do
    for b in {0,14}
    do

	filename="ntuple_signal_m${m}_width${b}.root" 
	echo "../"$filename
	echo "../"$filename > filelist.txt
	./analyzerdEdx
	mv analyzerdEdx_histograms.root analyzerdEdx_histograms_m${m}_width${b}.root
	


    done
    
done
