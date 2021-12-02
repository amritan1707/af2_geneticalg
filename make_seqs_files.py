def write_seq_csv_pd1(fname, helices, pd1_seqs):
    with open(fname, "w") as f:
        for helix in helices:
            for pd1_seq in pd1_seqs:
                f.write(","+pd1_seq+","+helix+"\n")


pd1_seqs = []
helix_seqs = []

with open("helix_seqs.txt", "r") as hf:
    with open("pd1_seqs.txt","r") as pf:

        for linepd1 in pf:
            pd1_seqs.append(linepd1.strip())
        for linehelix in hf:
            helix_seqs.append(linehelix.strip())


print(pd1_seqs)
print(helix_seqs)

desired_pd1_seqs = [pd1_seqs[2]]
desired_helix_seqs = []
for hs in helix_seqs:
    if len(hs)==18:
        desired_helix_seqs.append(hs)

with open("inputs/sequences.csv","w") as seqsf:
    for dps in desired_pd1_seqs:
        for dhs in desired_helix_seqs:
            seqsf.write(","+dps+","+dhs+"\n")
