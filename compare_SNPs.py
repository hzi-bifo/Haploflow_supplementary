#!/usr/bin/env python

import sys

haploflow_snps = sys.argv[1]
wastewater_snps = sys.argv[2]

no_snps = 0
per_sample = {}
with open(wastewater_snps, 'r') as wsnps:
    for line in wsnps:
        if line.startswith("sample"):
            continue
        sample, genome_pos, pos_cov, allele_count, ref_base, consensus_base, var_base, ref_freq, con_freq, var_freq, a_count, c_count, t_count, g_count = line.strip().split()
        if ((float(var_freq) > 0.1
            and ((var_base == 'A' and int(a_count) > 2)
            or (var_base == 'C' and int(c_count) > 2)
            or (var_base == 'T' and int(t_count) > 2)
            or (var_base == 'G' and int(g_count) > 2)))
            or ((consensus_base != ref_base and float(con_freq) > 0.1)
            and ((consensus_base == 'A' and int(a_count) > 2)
            or (consensus_base == 'C' and int(c_count) > 2)
            or (consensus_base == 'G' and int(g_count) > 2)
            or (consensus_base == 'T' and int(t_count) > 2)))):
            no_snps += 1
            if sample in per_sample:
                per_sample[sample].append((genome_pos, ref_base, var_base, pos_cov))
            else:
                per_sample[sample] = [(genome_pos, ref_base, var_base, pos_cov)]

unique_hap = set()
unique_var = set()
haplo_sample = {}
with open(haploflow_snps, 'r') as hsnps:
    for line in hsnps:
        try:
            ref, contig, pos, ref_base, var_base, contig_pos, sample = line.strip().split()
        except Exception as e:
            print(e)
            print(line.strip().split())
        if sample in haplo_sample:
            haplo_sample[sample].append((pos, ref_base, var_base))
        else:
            haplo_sample[sample] = [(pos, ref_base, var_base)]
        if sample in per_sample:
            if True not in [pos in x for x in per_sample[sample]]:
                unique_hap.add((sample, pos))

for sample in per_sample:
    for pos, ref, var, cov in per_sample[sample]:
        if True not in [pos in x for x in haplo_sample[sample]]:
            unique_var.add((sample, int(pos), int(cov)))

#print([(x, len(per_sample[x])) for x in per_sample])
#print([(x, len(haplo_sample[x])) for x in haplo_sample])
for elem in unique_hap:
    print(elem)
#print("X")
#for elem in sorted(unique_var):
#    print(elem)
print(sum([len(haplo_sample[x]) for x in haplo_sample]))
print(sum([len(per_sample[x]) for x in per_sample]))
