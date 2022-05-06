#!/usr/bin/env python

import sys
import os
import subprocess
import math
from sklearn.cluster import KMeans as kmeans
import numpy as np
import matplotlib.pyplot as plt

contig_file = sys.argv[1] # HaploFlow fasta file
reference_file = sys.argv[2]
snp_file = sys.argv[3]
coords_file = sys.argv[4]
duplication_ratio_file = sys.argv[5]
bam = sys.argv[6]
out = sys.argv[7]

with open(duplication_ratio_file, 'r') as dr:
    for line in dr:
        if (not line.startswith("Duplication ratio")):
            continue
        name, duplication_ratio = line.strip().split('\t')
        break

duplication_ratio = float(duplication_ratio)
if duplication_ratio - math.floor(duplication_ratio) > 0.15:
    nr_clusters = math.ceil(duplication_ratio)
else:
    nr_clusters = math.floor(duplication_ratio)
min_len = 500
seed = 1234

flow_map = {} # map contigs to their flow
with open(contig_file, 'r') as c:
    for line in c:
        if (not line.startswith('>')):
            length = len(line)
            if (length < min_len):
                del flow_map[nr]
            continue
        contig, nr, f, flow, cc, ccnr = line.strip().split("_")
        flow_map[nr] = float(flow)
      
# sklearn cluster
data = np.fromiter(flow_map.values(), dtype=float)
reshaped = data.reshape(-1,1)
if (len(flow_map) > 1):
    clusters = kmeans(n_clusters = int(nr_clusters), random_state = seed).fit(reshaped)
    labels = clusters.labels_
    centers = clusters.cluster_centers_
else:
    labels = [0]

i = 0
label_map = {}
for nr in flow_map:
    for i in range(len(data)):
        if (flow_map[nr] == data[i]):
            label_map[nr] = labels[i]

# make sure that the highest flow is always the "0" label
maxf = 0.
maxnr = -1
for nr in label_map:
    flow = flow_map[nr]
    if flow > maxf:
        maxf = flow
        maxnr = nr

if label_map[maxnr] == 1:
    for nr in label_map:
        if label_map[nr] == 0:
            label_map[nr] = 1
        else:
            label_map[nr] = 0

def get_next_min(covered, val):
    temp = []
    for v in covered:
        if (v[0] >= val):
            temp.append(v)
    return min(temp)

def get_covered(coords_file):
    covered = []
    # get the positions of the contigs
    with open(coords_file, 'r') as coords:
        for line in coords:
            if line.startswith("="):
                continue
            sep = line.strip().split()
            if (sep[0].startswith("[")):
                continue
            start = int(sep[0])
            end = int(sep[1])
            contig = [s for s in sep if "flow" in s][0]
            c, nr, f, flow, cc, ccnr = contig.split("_") 
            if label_map[nr] == 1:
                covered.append((start, end))

    cov_sorted_start = sorted(covered) # sorts by start positions
    cov_sorted_end = sorted(covered, key = lambda x: x[1])
    next_start = 0
    next_end = 0
    min_v = 0
    covered_final = []
    while next_end != max(cov_sorted_end, key = lambda x: x[1])[1]:
        next_start, next_end = get_next_min(cov_sorted_start, next_end)
        for v in cov_sorted_start:
            if (v[0] <= next_end and v[1] >= next_end):
                next_end = v[1]
        covered_final += [next_start, next_end]
    covered = covered_final
    return covered

if (nr_clusters > 1):
    covered = get_covered(coords_file)
else:
    covered = []

def is_covered(position, covered):
    i = 0
    for c in covered:
        if (i % 2 == 0):
            if (position >= c and position <= covered[i + 1]):
                return (True, c, covered[i + 1])
            elif (position < c and i > 0):
                return (False, covered[i - 1], c)
            elif (position < c and i == 0):
                return (False, 0, c)
        i += 1
    return (False, covered[-1], -1)

positions = []
snp_map = {}
with open(snp_file, 'r') as snps:
    for line in snps:
        snp_list = line.strip().split('\t')
        pos = int(snp_list[2]) #changed in latest QUAST!
        orig = snp_list[3]
        new = snp_list[4]
        contig = snp_list[1]
        c, nr, f, flow, cc, ccnr = contig.split("_") 
        if pos in snp_map:
            snp_map[pos].append((nr, orig, new))
        else:
            positions.append(pos)
            snp_map[pos] = [(nr, orig, new)]

ref = ""
with open(reference_file, 'r') as reference:
    for line in reference:
        if (line.startswith('>')):
            continue
        ref += line.strip()

positions = sorted(list(set(positions)))
js = [0 for i in range(nr_clusters)]
out_file = os.path.join(out, "strains_cds.fa")
#with open(out_file, 'a+') as o:
#    o.write(">Wuhan_reference\n")
#    o.write(ref)

def no_homopolymer(ref, pos, length):
    for i in range(length):
        seq_pre = ref[pos - length + i + 1:pos + i + 1]
        if (seq_pre != []):
            if (seq_pre == len(seq_pre) * seq_pre[0]): # consists of only one char
                return False
    return True

def write(to_write, ref, js, pos, new, label, prev_v, covered, indel):
    if (label == 0):
        if (indel == 0):
            ret = to_write[:pos - 1 + js[0]] + to_write[pos + js[0]:]
            js[0] -= 1
        elif (indel == 1):
            ret = to_write[:pos - 1 + js[0]] + new + to_write[pos - 1 + js[0]:]
            js[0] += 1
        elif (indel == 2):
            ret = to_write[:pos - 1 + js[0]] + new + to_write[pos + js[0]:]
        else:
            print("indel must be 0,1 or 2")
            raise ValueError
    else:
        is_covd, start, end = is_covered(pos, covered)
        prev = prev_v[0]
        prev_t = prev_v[1]
        prev_covd, prev_start, prev_end = is_covered(prev, covered) #both times the bool must be true
        #print("Previous/Current: %s/%s" % (prev_t,indel))
        #print("%s < %s" % (prev_end, pos + js[1]))
        if (is_covd and prev_covd or prev == 0):
            if (prev_t == 0):
                ret = to_write[:prev + js[1]]
                nextp = prev
            elif (prev_t == 1):
                ret = to_write[:prev + 2 + js[1]]
                nextp = prev + 2
            else:
                ret = to_write[:prev + 1 + js[1]]
                nextp = prev + 1
            if (indel == 0):
                if (prev_end < pos + js[1]):
                    ret += ref[nextp:prev_end] + to_write[prev_end + js[1]:start] + ref[start - js[1]:pos - 1] + ref[pos:]
                else:
                    ret += ref[nextp:pos - 1] + ref[pos:]
                js[1] -= 1
            elif (indel == 1):
                if (prev_end < pos + js[1]):
                    ret += ref[nextp:prev_end] + to_write[prev_end + js[1]:start] + ref[start - js[1]:pos] + new + ref[pos - 1:]
                else:
                    ret += ref[nextp:pos - 1] + new + ref[pos - 1:]
                js[1] += 1
            elif (indel == 2):
                if (prev_end < pos + js[1]):
                    ret += ref[nextp:prev_end] + to_write[prev_end + js[1]:start] + ref[start - js[1]:pos - 1] + new + ref[pos:]
                else:
                    ret += ref[nextp:pos - 1] + new + ref[pos:]
            else:
                print("indel must be 0,1 or 2")
                raise ValueError
        else:
            print(prev)
            print(pos)
            print(covered)
            print(is_covered(pos, covered))
            print(is_covered(prev, covered))
            print("Last variant was not covered")
            raise ValueError
    return (ret, js)

# only works for SARS-CoV-2
def strand_bias(ref, alt, pos, bam):
    proc = subprocess.Popen("samtools mpileup -f %s -r \"NC_045512.2:%s-%s\" -d 80000 %s" % (ref, pos, pos, bam), shell=True, stdout=subprocess.PIPE)
    out = proc.stdout.read()
    code = out.split('\t')[4] # . = match forward , = match reverse CHAR = mismatch forward char = mismatch reverse
    ref_f = code.count(".")
    ref_r = code.count(",")
    if (alt == "."):
        alt = "*"
    alt_f = code.count(alt.upper())
    alt_r = code.count(alt.lower())
    return (ref_f, ref_r, alt_f, alt_r)
    

prev_ref = ref
for i in range(nr_clusters):
    #out_file = os.path.join(out, "%s_%s.fa" % (os.path.split(contig_file)[-2],i))
    print(i)
    to_write = prev_ref
    prev = (0,2)
    for pos in positions:
        snps = snp_map[pos]
        for k in range(len(snps)):
            nr, orig, new = snps[k]
            label = label_map[nr]
            if label == i:
                if (orig != "." and new != "." and to_write[pos - 1 + js[i]] != orig):
                    #if (to_write[pos - 1 + js[i]] != new):
                    #    print("mismatched SNP for contig %s (label %s)" % (nr,label))
                    #    print("Expected %s at pos %s, got %s" % (orig, pos, to_write[pos - 1 + js[i]]))
                    #    print("%s %s %s" % (ref[pos - 5: pos - 1], ref[pos - 1], ref[pos: pos + 4]))
                    #    print("%s %s %s" % (to_write[pos - 5 + js[i]: pos - 1 + js[i]], to_write[pos - 1 + js[i]], to_write[pos + js[i]: pos + 4 + js[i]]))
                    #else:
                    #    print("Already replaced for contig %s (%s -> %s) at pos %s" % (nr, orig, new, pos))
                    #    #print("Sequence: %s" % to_write[pos - 5 + js[i]: pos + 4 + js[i]])
                    #    #print("Original: %s" % (ref[pos - 5: pos + 4])) 
                    continue
                else:
                    if (new == "."):
                        if (no_homopolymer(ref, pos, 4)):
                            #print("deleted %s at pos %s for contig %s" % (orig, pos, nr))
                            to_write, js = write(to_write, ref, js, pos, ".", label, prev, covered, 0)
                            #print(pos)
                            #print(strand_bias(reference_file, new, pos, bam))
                            prev = (pos, 0)
                        else:
                            pass
                            #print("did not perform deletion of %s in homopolymer at pos %s for contig %s" % (orig, pos, nr))
                            #print("sequence is %s" % (ref[pos - 4 + js[i]: pos + 5 + js[i]]))
                    elif (orig == "."):
                        if (no_homopolymer(ref, pos, 4)):
                            #print("added %s at pos %s for contig %s" % (new, pos, nr))
                            to_write, js = write(to_write, ref, js, pos, new, label, prev, covered, 1)
                            #print(pos)
                            #print(strand_bias(reference_file, new, pos, bam))
                            prev = (pos, 1)
                        else:
                            pass
                            #print("did not perform addition of %s in homopolymer at pos %s for contig %s" % (new, pos, nr))
                            #print("sequence is %s" % (ref[pos - 4 + js[i]: pos + 5 + js[i]]))
                    else:
                        if (no_homopolymer(ref, pos, 4)):
                            #print("replaced %s with %s at %s for contig %s" % (orig, new, pos, nr))
                            to_write, js = write(to_write, ref, js, pos, new, label, prev, covered, 2)
                            #print(pos)
                            #print(strand_bias(reference_file, new, pos, bam))
                            prev = (pos, 2)
                        else:
                            pass
                            #print("did not perform change of %s to %s in homopolymer at pos %s for contig %s" % (orig, new, pos, nr))
                            #print("sequence is %s" % (ref[pos - 4 + js[i]: pos + 5 + js[i]]))
    with open(out_file, 'a+') as o:
        o.write(">%s_%s\n" % (os.path.split(contig_file)[-1][:-3],i)
        o.write(to_write)
        o.write('\n')
    prev_ref = to_write
