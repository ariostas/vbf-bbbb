
./mergeSelFiles.sh vbf /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/HHToBBBB_14TeV/ /afs/cern.ch/work/a/ariostas/public/vbf-bbbb/

mv /afs/cern.ch/work/a/ariostas/public/vbf-bbbb/vbf.root /afs/cern.ch/work/a/ariostas/public/vbf-bbbb/HHToBBBB_14TeV.root

while read line
do
  array=($line)
  
   if [ "${array[0]}" != "#" ]; then 

	./mergeSelFiles.sh ${array[0]} /afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/${array[0]}/ /afs/cern.ch/work/a/ariostas/public/vbf-bbbb/   

    fi

done < xsec.txt
