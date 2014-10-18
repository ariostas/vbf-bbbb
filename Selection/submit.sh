
rm -r /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/*
rm -r LSFJOB_*

mkdir -v /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV
cd /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV
/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select ls /store/group/upgrade/dihiggs_signal_4b/VBFHHTobbbb_TuneZ2_14TeV-madgraph/files/ | grep root > "HHToBBBB_14TeV.txt"
split -l 100 HHToBBBB_14TeV.txt
rm HHToBBBB_14TeV.txt
cd -
for file in /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV/*
do
    
    bsub -q 8nh -W 480 -J Signal run.sh selectDelphes.C $file HHToBBBB_14TeV 0.0006691873 1

done

#mkdir -v /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV
#cd /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV
#/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select ls /store/group/upgrade/delphes/dihiggs_signal_bbbb/gFHHTobbbb_TuneZ2_14TeV_madgraph/files/ | grep root > "HHToBBBB_14TeV.txt"
#split -l 150 HHToBBBB_14TeV.txt
#rm HHToBBBB_14TeV.txt
#cd -
#for file in /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV/*
#do
    
#    bsub -q 8nh -W 480 -J Signal run.sh selectDelphes.C $file HHToBBBB_14TeV 0.01331716 1

#done

while read line #loop over lines in ${conf_file}
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

	    mkdir -v /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/${array[0]}
        cd /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/${array[0]}
        /afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select ls /store/group/upgrade/delphes/PhaseII_140PU_ProdJul28/${array[0]}/ | grep root > "${array[0]}.txt"
        split -l 100 "${array[0]}.txt"
        rm "${array[0]}.txt"
        cd -
        for file in /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/${array[0]}/*
        do
	        bsub -q 8nh -W 480 -J ${array[0]} run.sh selectDelphes.C $file ${array[0]} ${array[1]} 0
        done

	fi
done < xsec.txt

