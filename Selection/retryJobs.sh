
for file in /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp_test/HHToBBBB_14TeV/x*
do
    
    bsub -q 8nh -W 480 -J Signal run.sh selectDelphes.C $file HHToBBBB_14TeV 0.0006691873 1

done


#for file in /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV/*
#do
    
#    bsub -q 8nh -W 480 -J Signal run.sh selectDelphes.C $file HHToBBBB_14TeV 0.01331716 1

#done

while read line
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

        for file in /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp_test/${array[0]}/x*
        do
	        bsub -q 8nh -W 480 -J ${array[0]} run.sh selectDelphes.C $file ${array[0]} ${array[1]} 0
        done

	fi
done < xsec.txt

