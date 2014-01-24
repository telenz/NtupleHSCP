#!/bin/bash

oldMass=100
newMass=100
oldWidth=0
for m in {1..7}
do
    for b in {0..6}
    do

	if [ $b -eq 6 ];then
	    if [ $m -eq 1 ];then
		break 
	    fi
	fi

	if [ ${oldWidth} -eq 6 ];then
	    newWidth=0
	else
	    newWidth=$((${oldWidth}+1))
	fi

        ./analyzer
	echo ""
	echo "ntuple_m${oldMass}_width${oldWidth}.root"
	echo ""
	mv analyzer_histograms.root analyzer_histograms_m${oldMass}_width${oldWidth}.root
        for i in filelist.txt ; sed -i "s/m${oldMass}_width${oldWidth}.root/m${newMass}_width${newWidth}.root/g" "$i"
	#echo "ersetze dies"
	#echo "m${oldMass}_width${oldWidth}.root"
	#echo "mit dem"
	#echo "m${newMass}_width${newWidth}.root"
	#echo ""
	oldWidth=${newWidth}

	if [ ${oldWidth} -eq 0 ];then
	    oldMass=${newMass}
	fi
	
	echo "ntuple_m${newMass}_width${newWidth}.root"
    done

    oldMass=${newMass}
    newMass=$((${oldMass}+100))
    
done
./analyzer
mv analyzer_histograms.root analyzer_histograms_m${oldMass}_width${newWidth}.root
echo "ntuple_m${oldMass}_width${newWidth}.root"
for i in filelist.txt ; sed -i "s/m${oldMass}_width${oldWidth}.root/m100_width0.root/g" "$i"