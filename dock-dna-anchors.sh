#!/usr/bin/env bash

d=`pwd`

protein=$1
aromatics=$2 #aromatics.list
np=$3   #number of CPU to use
nmodels=$4 #number of models per motif

################################################################################
# Check environment variables
################################################################################
if [ -z $ATTRACTDIR ];then
  echo "Define the ATTRACTDIR environment variable with:"
  echo "export ATTRACTDIR={your path to attract/bin}"
  echo "or add that line in your /home/.bashrc"
  exit
fi
if [ -z $ATTRACTTOOLS ];then
  echo "Define the ATTRACTTOOLS environment variable with:"
  echo "export ATTRACTTOOLS={your path to attract/tools}"
  echo "or add that line in your /home/.bashrc"
  exit
fi
if [ -z $TRIDNALIB ];then
  echo "Define the TRIDNALIB environment variable with:"
  echo "export TRIDNALIB={your path}/triDNAlib/ssDNA"
  echo "or add that line in your /home/.bashrc"
  exit
fi

echo '**************************************************************'
echo 'Convert receptor protein and reference dna'
echo '**************************************************************'

####Use "if false;then" and "fi" lines to disable parts of the protocol
###if false;then

python2 $ATTRACTDIR/../allatom/aareduce.py $protein protein-aa.pdb --dna --chain A --heavy --pdb2pqr > mapping
python2 $ATTRACTTOOLS/reduce.py protein-aa.pdb proteinr.pdb --dna --chain A > mapping
sed -i 's/25   0.000/20   0.000/g' proteinr.pdb

if [ -s dna.pdb ];then 
	python2 $ATTRACTDIR/../allatom/aareduce.py $dna dna-aa.pdb --dna --heavy > /dev/null
	python3 split-frag.py dna-aa.pdb frag > frag.list
fi
###fi
receptorr=proteinr.pdb

echo '**************************************************************'
echo 'Create restraints files'
echo '**************************************************************'

if [ ! -d restraints ];then
	mkdir restraints
fi
python3 create_stack_restraints.py $aromatics proteinr.pdb mapping restraints 5 6 > restraints.list

#restraintslist=restraints.list
# !!!!!!!!!!!!! This is to test the script on the "correct" ad hoc restraints !!!!!!!!
restraintslist=native-rest.list


echo '**************************************************************'
echo 'Generate starting positions'
echo '**************************************************************'
# Produce $nstart staring positions (translation+rotation) around receptor

nstart=100
python2  $ATTRACTTOOLS/randsearch.py 2 $nstart --fix-receptor --radius 50 > randsearch_$nstart.dat


######################################################################
dock(){
######################################################################
	motif=$1
	np=$2
	nmodels=$3

  ## Define dna ligand:
	ln -s $TRIDNALIB/${motif}
	list=$TRIDNALIB/${motif}r.list
	listaa=$TRIDNALIB/${motif}-aa.list
	ligandr=`head -n 1 $list`
	ligandaa=`head -n 1 $listaa`
	Nconf=`cat $list|wc -l`
	
	## Define general docking parameters:

	gridparams=" --grid 1 receptorgrid.gridheader"
	common="$receptorr $ligandr --ens 2 $list --np $np --chunks $np"
	
	# Disable RNA-protein force-field ("ghost") and translations ("only-rot"),
	# Minimise energy (max 100 iterations) only with restraints and rotations
	#    => Rotate the ligand toward the anchoring amino-acids.
	ghostparams="$ATTRACTDIR/../attract.par $common --vmax 100 --ghost --only-rot"
	
	# Minimise with both force-field and restraints, on precomputed grid.
	params="$ATTRACTDIR/../attract.par $common --fix-receptor --vmax 100 $gridparams"
	
	# Score the finale docking poses, without grid.
	scoreparams="$ATTRACTDIR/../attract.par $common --score --rcut 50"


	echo '**************************************************************'
	echo 'calculate receptorgrid grid'
	echo '**************************************************************'
	## Acceleration by pre-computed grid:

	$ATTRACTDIR/shm-clean #delete previous grid
		
	#print what coarse-grain beads (pseudo atom types) are in the ligand.
	awk '{print substr($0,58,2)}' $ligandr | sort -nu > receptorgrid.alphabet
	
	#Precompute forces felt by each ligand bead if placed on each point of a grid around the receptor.
	#During docking, forces will be computed from values on surrounding points,
	# instead on summing over on all receptor-ligand bead-bead pairs.
	$ATTRACTDIR/make-grid-omp $receptorr $ATTRACTDIR/../attract.par 5.0 7.0 receptorgrid.gridheader  --shm --alphabet receptorgrid.alphabet

	echo '**************************************************************'
	echo 'calculate starting poses'
	echo '**************************************************************'

	# All fragments of the library are placed on each starting position
	if [ ! -s randsearch_$nstart_$Nconf.dat ];then
		python2  $ATTRACTTOOLS/ensemblize.py randsearch_$nstart.dat $Nconf 2 all > randsearch_$nstart-$Nconf.dat
	fi
	start=randsearch_$nstart-$Nconf.dat

	# The next 2 lines are to quickly test the script
	# with a subset of the starting positions.
	# If the test runs fine, then remove the 2 next lines:
	#$ATTRACTTOOLS/top randsearch_$nstart-$Nconf.dat 100 > top.dat
	#start=top.dat

	echo '**************************************************************'
	echo "docking for each set of restraints"
	echo '**************************************************************'
	set -e
  i=0
	echo "dock"
  nrest=`cat $restraintslist|wc -l`
 
	for r in `seq $nrest`; do
		echo "rest=`awk -v r=$r 'NR==r{print $1}' $restraintslist` "
		
		rest=`awk -v r=$r 'NR==r{print $1}' $restraintslist` #get path to restraints file

		echo '**************************************************************'
		echo "sampling ${motif} rest$rest"
		echo '**************************************************************'
		python2  $ATTRACTDIR/../protocols/attract.py $start $ghostparams --rest $rest --output ${motif}-0.dat
		$ATTRACTDIR/fix_receptor ${motif}-0.dat 2 --ens 0 $Nconf > ${motif}-1.dat
		python2  $ATTRACTDIR/../protocols/attract.py ${motif}-1.dat $params --rest $rest --output ${motif}-2.dat

		echo '**************************************************************'
		echo "scoring ${motif} rest$rest"
		echo '**************************************************************'
		python2  $ATTRACTDIR/../protocols/attract.py ${motif}-2.dat $scoreparams --rest $rest --output ${motif}-2.score

		echo '**************************************************************'
		echo "sorting ${motif} rest$rest"
		echo '**************************************************************'
		python2  $ATTRACTTOOLS/fill-energies.py ${motif}-2.dat ${motif}-2.score > ${motif}-3.dat
		python2  select-dat.py ${motif}-3.dat --score 0 > ${motif}-4.dat
		python2  $ATTRACTTOOLS/sort.py ${motif}-4.dat > ${motif}-5.dat
		$ATTRACTDIR/fix_receptor ${motif}-5.dat 2 --ens 0 $Nconf | python2  $ATTRACTTOOLS/fill.py /dev/stdin ${motif}-5.dat > ${motif}-6.dat

		echo '**************************************************************'
		echo "removing redundant poses for ${motif} rest$rest"
		echo '**************************************************************'
		$ATTRACTDIR/deredundant ${motif}-6.dat 2 --ens 0 $Nconf --lim 0.2 | python2  $ATTRACTTOOLS/fill-deredundant.py /dev/stdin ${motif}-6.dat > $motif-rest$r.dat
		
		echo '**************************************************************'
		echo "convert ${motif} rest$rest to pdb"
		echo '**************************************************************'		
		$ATTRACTTOOLS/top $motif-rest$r.dat $nmodels >  $motif-top$nmodels.dat
		$ATTRACTDIR/collect $motif-top$nmodels.dat /dev/null $ligandaa --ens 2 $listaa > $motif-rest$r-top$nmodels.pdb
		
		echo '**************************************************************'
		echo 'compute RMSD'
		echo '**************************************************************'
		if [ -s frag.list ];then
			i=1
			for frag in `cat frag.list`; do
				python2 $ATTRACTDIR/lrmsd.py $motif-rest$r.dat $frag $frag --ens 2 $listaa --allatoms |awk '{print NR, $2}'> $motif-rest$r-frag$i.rmsd
				
				echo "rank  RMSD"
				sort -nk2 $motif-rest$r-frag$i.rmsd|head
				
			i=$(($i+1))
			done
		fi
	
	done

rm receptorgrid.gridheader
#rm ${motif}-[1-5].dat
}



#create a text list of the 3-nt motifs you want to dock,
#one motif per line.
#Add each motif only once in the list, even if it is
#present several times in the full sequence.
#ex: for AUUUUUG, the list is (AUU, UUU, UUG.)


for motif in `cat motifs.list`;do
	echo "docking $motif"
	dock ${motif} $np $nmodels
done
