#!/bin/bash
#------------------------------------------------------------
# Submit a batch of jobs to lxbtch
#
# example command:
#             root_script     conf_file output_location
# ./submit.sh selectDelphes.C xsec.txt  /afs/cern.ch/work/a/ariostas/private/HHToGGBB/
# 
# conf_file has format (no leading "#")
# sample_type cross_section
#
# Note: conf file needs empty last line... clunky I know.
# Jay Lawhorn 11/4/13
#------------------------------------------------------------

root_script=$1
  conf_file=$2

rm *_14TEV.txt
rm -r /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/*
rm -r LSFJOB_*

mkdir -v /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV
bsub -q 1nd -W 1440 -J Signal run.sh ${root_script} HHToBBBB_14TeV 0.000599

while read line #loop over lines in ${conf_file}
do
  array=($line)
  #for conf in PhaseI/Configuration0 PhaseII/Configuration3 PhaseII/Configuration4v2 # to process all three configurations
  for conf in /store/group/upgrade/delphes/ProdJun14 # or just one configuration
  do
    if [ "${array[0]}" != "#" ]; then 
	# get list of files in eos for that sample+configuration combination

	      mkdir -v /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/${array[0]}
	      /afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select ls ${conf}/${array[0]}/ | grep root > "${array[0]}.txt"
	      bsub -q 1nd -W 1440 -J ${array[0]} run.sh ${root_script} ${array[0]} ${array[1]} # submit to lxbtch

    fi
  done
done < ${conf_file}

