#!/usr/bin/env python3

import sys, os
pdbfile = sys.argv[1]
output_file_pattern = sys.argv[2]

from collections import deque
fragment_file_names = []
fragment_files = deque([],3)
fragment_counter = 0

def new_fragment_file():
    global fragment_counter
    if len(fragment_files) == 3:
        f = fragment_files.popleft()
        f.close()
    fragment_counter += 1
    fname = "{}-{}.pdb".format(output_file_pattern, fragment_counter)
    f = open(fname, "w")
    fragment_files.append(f)
    fragment_file_names.append(fname)

last_resid = None
for l in open(pdbfile):
    if not l.startswith("ATOM"):
        continue
    resid = l[21:26]
    if resid != last_resid:
        new_fragment_file()
        last_resid = resid
    for f in fragment_files:
        f.write(l)

for f in fragment_files:
    f.close()
for fname in fragment_file_names[-2:]:
    os.remove(fname)
for fname in fragment_file_names[:-2]:
    print(fname)
