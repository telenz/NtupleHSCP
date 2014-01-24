#!/bin/bash


for m in {100,200,300,400,500,600,700,800}
do
    for b in {0..14}
    do

      printf "%s\n" /nfs/dust/cms/user/tlenz/mc/debug/test10/pMSSM12_MCMC1_30_549144_m${m}_width${b}_*.root > sampleList.txt

      k=0
      while read line1; do 
	  #echo ${line1}
	  while read line; do  
	      if [[ $line == *cms.untracked.vstring\(\"file:* ]]
		  then
		  echo  'process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:'${line1}'"))'
	      else
		  echo $line
	      fi
	  done < MyNtuple_cfg.py > newFile.py

	  cmsRun newFile.py
	  mv ntuple.root /nfs/dust/cms/user/tlenz/mc/HSCPntuples/SamplesFULL/ntuple_m${m}_width${b}_${k}.root
	  (( k+=1 ))

      done < sampleList.txt
    done
    
done

