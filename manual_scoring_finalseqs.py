from score_pd1 import score_pd1
from run_geneticalg import run_af2_pd1
import os
import subprocess

path = "/home/amrita/pd1/run_geneticalg/starting_odessa/iter_20_all_var/"
final_seqs = []
pd1seq = "WNPPTFSPALLVVTEGDNATFTCSFSNTSESFHVVWHRESPSGQTDTLAAFPEDRSQPGQDSRFRVTQLPNGRDFHMSVVRARRNDSGTYVCGVISLAPKIQIKESLRAELRVTERRAE"
pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
with open(path+"final_seqs.txt","r") as inpf:
    for l in inpf:
        final_seqs.append(l.strip())

#run_af2_pd1(final_seqs, pd1seq, path+"final_sequences/", "/home/amrita/pd1/", "flags_froome.txt", "/home/nzrandol/alphafold/run/")
oppath = path+"final_sequences/outputs/"
files = os.listdir(oppath)
for f in files:
    if f.endswith("pbz2"):
        subprocess.run(["bunzip2", oppath+f])

files = os.listdir(oppath)
final_scores = []
for f in files:
    if f.endswith("out"):
        seqnum = int(f.split("_")[1])
        pdbf = f.partition('results')[0]+"unrelaxed.pdb"
        print(f, pdbf)
        contacts, contactscore, confscore = score_pd1(oppath+pdbf, oppath+f, pocketresidues)
        final_scores.append((-contactscore*100 + confscore, -contactscore*100, confscore))
print(final_seqs, final_scores)

with open(oppath+"final_seqs_scores.txt", "w") as opfile:
    for seq, score in zip(final_seqs, final_scores):
        opfile.write(str(seq)+"\t"+str(score)+"\n")
