
./mergeSelFiles.sh vbf /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp_test/HHToBBBB_14TeV/ /afs/cern.ch/work/a/ariostas/public/vbf-bbbb_test/

mv /afs/cern.ch/work/a/ariostas/public/vbf-bbbb_test/vbf.root /afs/cern.ch/work/a/ariostas/public/vbf-bbbb_test/HHToBBBB_14TeV.root

while read line
do
  array=($line)
  
   if [ "${array[0]}" != "#" ]; then 

	   ./mergeSelFiles.sh ${array[0]} /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp_test/${array[0]}/ /afs/cern.ch/work/a/ariostas/public/vbf-bbbb_test/   

    fi

done < xsec.txt
