for i in `cat anchors.list`; do 
	grep "A  "$i 4OWX-protr.pdb |egrep "CSE TRP|CSF [TYR|PHE]"
done


