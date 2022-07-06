#!/usr/bin/env python3

import sys
from itertools import permutations, combinations

def combi(items, n):
	return [list(i) for i in combinations(items, n)]

def perm(items, n):
	return [list(i) for i in permutations(items, n)]

def printrest(toprint, alist, nlist, f, dist1, dist2):
	for l in toprint:
	  print(l.strip(), file=f),
	print("", file=f),
	for a,n in zip(alist, nlist):
		for k in [1, 2, 3]:
			print("%s1 %s%i 1 %.1f 2"%(a, n, k, dist1), file=f)
			print("%s2 %s%i 1 %.1f 2"%(a, n, k, dist2), file=f)

###############################################################"
#file with residues selection
anchres = [l.split() for l in open(sys.argv[1]).readlines()]

#protein pdb file in coarse-grain format
prot = [l for l in open(sys.argv[2]).readlines()]

#residues renumbering
mm = [l.split() for l in open(sys.argv[3]).readlines()]

#directory for restraints files
d = sys.argv[4] 

#max distances 
dist1 = float(sys.argv[5]) 
dist2 = float(sys.argv[6]) 
###############################################################"

#convert residue numbering
mapping = { m[0]:m[1] for m in mm}
aares = [y for x in anchres for y in x]
cgres = [mapping[r] for r in aares]

toprint=[]

#print protein anchoring beads
i = 1
anchors=[]
for r, res in enumerate(cgres):
	ll = [l for l in prot if l[22:26].strip() == res]
	if len(ll) == 0:
		raise AssertionError("residu %s does not exist"%aares[r]) 
	
	aa = ll[0][17:20]
	if aa not in ("TRP","PHE","TYR"):
		raise AssertionError("residu %s is not aromatic"%aares[r]) 
		
	bead1, bead2 = "CSE", "CSF"
	if ll[0][17:20] == "TRP":
		bead1, bead2 = "CG ", "CSE"
	at1 = [l[5:11].strip() for l in ll if l[13:16] == bead1][0]
	at2 = [l[5:11].strip() for l in ll if l[13:16] == bead2][0]
	toprint.append("a%i1 1 %i"%(i,int(at1)))
	toprint.append("a%i2 1 %i"%(i,int(at2)))
	anchors.append("a%i"%i)
	i+=1
#print dna anchored beads
nprot = len(prot)
j=4
for n in [1,2,3]:
	for b in [1,2,3]:
		toprint.append("n%i%i 1 %i"%(n,b,nprot+j))
		j+=1
	j+=3 # !!!! replace by 4 for purines !!!!!!

nucl = ["n1", "n2", "n3"]

for a in anchors:
	for n in nucl:
		filename = "%s/rest-%s%s"%(d,a,n)
		ff = open(filename, 'w')
		print(filename)
		printrest(toprint, [a], [n], ff, dist1, dist2)
		ff.close()

for [a1, a2] in combi(anchors, 2):
	for [n1, n2] in perm(nucl, 2):
		filename="%s/rest-%s%s-%s%s"%(d,a1,n1,a2,n2)
		print(filename)
		ff = open(filename, 'w')
		printrest(toprint, [a1, a2], [n1, n2], ff, dist1, dist2)
		ff.close()

for [a1, a2, a3] in combi(anchors, 3):
	for [n1, n2, n3] in perm(nucl, 3):
		filename="%s/rest-%s%s-%s%s-%s%s"%(d,a1,n1,a2,n2,a3,n3)
		print(filename)
		ff = open(filename, 'w')
		printrest(toprint, [a1, a2, a3], [n1, n2, n3], ff, dist1, dist2)
		ff.close()

