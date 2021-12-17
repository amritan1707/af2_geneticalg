import random
from Sequence import Sequence
from score_pd1 import score_pd1, run_af2_pd1
import glob
import os
import subprocess

path = "../run1/"
oppath = path+"outputs/"
pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
files = os.listdir(oppath)
for f in files:
    if f.endswith("pbz2"):
        subprocess.run(["bunzip2", f])
files = os.listdir(oppath)
scored_seqs = {}
scoring_pool = [Sequence("PSREFLILALQIALTLKA"), Sequence("QFWNLLIYLMRVYLQKHA"), Sequence("EAKNILISLLIYWAQMLD"), Sequence("FMWNILVTIARVMAQQLD"), Sequence("TAWELLIKIARYMAQQLD"), Sequence("PKREFLILALLIALKLES"), Sequence("KLTEIMLSIGLVFMWRKS"), Sequence("PEETFHRLLWEYMERLLA"), Sequence("EEEELWIQFLRLALKIAL"), Sequence("AYEMFQILFMWYLEMKDA"), Sequence("SYERMIELMLKWLEKHLA"), Sequence("LEYLLWILAMQYLEKHLA"), Sequence("TEREVTELLKIWRELFMA"), Sequence("ECRLLHILHIRYAKAWTA"), Sequence("PSKNIFLSLAWWIAQVLT")]
pool = [Sequence("PSREFLILALQIALTLKA"), Sequence("QFWNLLIYLMRVYLQKHA"), Sequence("EAKNILISLLIYWAQMLD"), Sequence("FMWNILVTIARVMAQQLD"), Sequence("TAWELLIKIARYMAQQLD"), Sequence("PKREFLILALLIALKLES"), Sequence("KLTEIMLSIGLVFMWRKS"), Sequence("PEETFHRLLWEYMERLLA"), Sequence("EEEELWIQFLRLALKIAL"), Sequence("AYEMFQILFMWYLEMKDA"), Sequence("SYERMIELMLKWLEKHLA"), Sequence("LEYLLWILAMQYLEKHLA"), Sequence("TEREVTELLKIWRELFMA"), Sequence("ECRLLHILHIRYAKAWTA"), Sequence("PSKNIFLSLAWWIAQVLT")]
files = os.listdir(oppath)
for f in files:
    if f.endswith("out"):
        seqnum = int(f.split("_")[1])
        pdbf = f.partition('results')[0]+"unrelaxed.pdb"
        print(f, pdbf)
        contacts, contactscore, confscore = score_pd1(oppath+pdbf, oppath+f, pocketresidues)
        scored_seqs[scoring_pool[seqnum]] = (-contactscore*100 + confscore, -contactscore*100, confscore)

scored_pool = {}
for p in pool:
    scored_pool[p] = scored_seqs[p]

print(scored_pool)

sorted_scored_pool = sorted(scored_pool.items(), key=lambda x: x[1])
print(sorted_scored_pool)
sorted_scores = []
for s in sorted_scored_pool:
    print(str(s[0]))
    print(scored_pool[s[0]])

#seqs_per_iteration.append(sorted_scored_pool)
#scores_per_iteration.append(sorted_scores)

newpool = []
for sp in sorted_scored_pool[:round(len(sorted_scored_pool)/2)]:
    newpool.append(sp)


print(newpool)
for p in newpool:
    print(str(p))
