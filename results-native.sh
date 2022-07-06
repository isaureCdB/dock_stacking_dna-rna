motif=TTT
n=`cat native-rest.list|wc -l`
for i in `seq $n`; do
	r=`awk -v i=$i 'NR==i{print $1}' native-rest.list`
	f=`awk -v i=$i 'NR==i{print $2}' native-rest.list`
	echo "----------------- $f"
	echo "rank    rmsd"
	sort -nk2 $motif-rest$i-$f.rmsd|head -n 3
done
