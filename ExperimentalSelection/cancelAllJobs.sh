
eval bjobs > cancel.txt

while read line
do
  array=($line)

    if [ "${array[0]}" != "No" ]; then 
		  
	bkill ${array[0]}

    fi
done < cancel.txt

rm cancel.txt
