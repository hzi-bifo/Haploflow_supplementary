#!/usr/bin/env python

import sys
import re

inp = sys.argv[1]
out = sys.argv[2]

lengths = {}
with open(inp, 'r') as i:
    ct = 0
    for line in i:
        if (ct % 3 == 0): # node name
            to_set = line.strip()
        elif (ct % 3 == 1):
            snps = line.strip().split()
            count = 0
            k = 0
            for snp in snps:
                if k+1 < len(snps) and snp[-1] == "-" or snp[0] == "-": # treat longer indels as single mutation
                    if int(snps[k+1][1:-1]) - int(snp[1:-1]) == 1:
                        count += 1
                k += 1
            lengths[to_set] = len(line.strip().split()) - count
        else:
            pass
        ct += 1

print(lengths)

to_write = ""
with open(out, 'r') as o:
    for line in o:
        y = re.finditer('[\,,)]', line)
        z = re.finditer('[\,,()]', line)
        l1 = list(y)
        l2 = [it.span() for it in z]
        last = 0
        i = 0
        for it in l1:
            prev = l2.index(it.span()) - 1
            name = line[l2[prev][1]:it.span()[0]] # get name since last match
            if name in lengths:
                to_write += line[last: it.span()[0]]
                to_write += ":" + str(lengths[name])
            elif (name.isdigit()):
                to_write += line[last: it.span()[0] - len(name)]
                to_write += ":0.0"
            else:
                to_write += line[last: it.span()[0]] # write leaves
                to_write += ":0.1" # give them a little distance for style
            to_write += line[it.span()[0]:it.span()[1]]
            last = it.span()[1]
            i += 1
        to_write += ":0.0"
        to_write += ";"

print(to_write)
#with open(out, 'w') as o:
#    o.write(to_write)
