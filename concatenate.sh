motif=$1
concat-dat.py $motif\_rest*.dat > $motif-allrest.dat
$ATTRACTDIR/deredundant $motif-allrest.dat --lim 0.2 > $motif-allrest-dr02.dat
$ATTRACTDIR/deredundant $motif-allrest.dat --lim 0.5 > $motif-allrest-dr05.dat
$ATTRACTDIR/deredundant $motif-allrest-dr02.dat --lim 0.5 > $motif-allrest-dr02-dr05.dat
