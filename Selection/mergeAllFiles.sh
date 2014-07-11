#!/bin/bash
#------------------------------------------------------------
# Merge all ntuples
#
# Use the following command: ./mergeAllFiles.sh xsec.txt
#------------------------------------------------------------

  conf_file=$1

./mergeSelFiles.sh vbf /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV/ /afs/cern.ch/work/a/ariostas/private/vbf-bbbb/ 

while read line #loop over lines in ${conf_file}
do
  array=($line)
  #for conf in PhaseI/Configuration0 PhaseII/Configuration3 PhaseII/Configuration4v2 # to process all three configurations
  
   if [ "${array[0]}" != "#" ]; then 

	   ./mergeSelFiles.sh ${array[0]} /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/${array[0]}/ /afs/cern.ch/work/a/ariostas/private/vbf-bbbb/   

    fi

done < ${conf_file}
