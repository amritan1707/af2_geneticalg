import random
from Sequence import Sequence
from score_pd1 import score_pd1, run_af2_pd1
import glob
import os
import subprocess

def create_new_seqs(startseqs, num_seqs, crossoverpercent = 0.2):
    """takes starting sequences and creates a pool of size num_seqs by mutation and crossover"""
    pool = startseqs
    while len(pool)<num_seqs*(1-crossoverpercent):
        obj = random.choice(startseqs)
        newseq = obj.mutate(resnums = [0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16], restypes = ["alpha" for x in range(13)])
        newseq = Sequence("".join(newseq))
        if newseq not in pool:
            pool.append(newseq)
    #print(pool)
    while len(pool)<num_seqs:
        oldseqs = random.sample(startseqs, 2)
        newseq = oldseqs[0].crossover(oldseqs[1])
        newseq = Sequence("".join(newseq))
        if newseq not in pool:
            pool.append(newseq)

    return pool

def run_genetic_alg_pd1(pa, pd1seq, startingseqs, pocketresidues, poolsize = 50, num_iter = 20):
    scored_seqs = {}
    curr_iter = 1
    while curr_iter <= num_iter:
        path = pa + "run"+str(curr_iter)+"/"
        if curr_iter == 1:
            pool = create_new_seqs(startingseqs, poolsize)

        else:
            pool = create_new_seqs(pool, poolsize)


        scoring_pool = [p for p in pool if p not in scored_seqs.keys()]
        #uncomment when running
        #def run_af2_pd1(pool, pd1seq, directory, fpath, flagsfile, af2path):
        run_af2_pd1(scoring_pool, pd1seq, path, "/home/amrita/pd1/", "flags_froome.txt", "/home/nzrandol/alphafold/run/") 
        print("done running af2")

        oppath = path+"outputs/"

        files = os.listdir(oppath)
        for f in files:
            if f.endswith("pbz2"):
                subprocess.run(["bunzip2", oppath+f])
        files = natsort.natsorted(os.listdir(oppath))
        print("done unzipping results files")

        for f in files:
            if f.endswith("out"):
                seqnum = int(f.split("_")[1])
                pdbf = f.partition('results')[0]+"unrelaxed.pdb"
                print(f, pdbf)
                contacts, contactscore, confscore = score_pd1(oppath+pdbf, oppath+f, pocketresidues)
                scored_seqs[scoring_pool[seqnum]] = contactscore*100 - confscore

        print("done scoring sequences")

        scored_pool = {}
        for p in pool:
            scored_pool[p] = scored_seqs[p]
        
        print(scored_pool)

        sorted_scored_pool = sorted(scored_pool, key=scored_pool.get, reverse=True)
        newpool = []
        for sp in sorted_scored_pool[:round(len(sorted_scored_pool)/2)]:
            newpool.append(sp)

        pool = newpool
        curr_iter+=1

    for p in pool:
        print(str(p))

    return pool


def generate_random_seqs(aalist, num_aas, num_seqs):
    oplist = []
    for i in range(num_seqs):
        oplist.append(Sequence("".join(random.choices(aalist, k=num_aas))))

    return oplist
    

if __name__=="__main__":
    startingseqs = [Sequence("PSREFLILALQIALTLKA"), Sequence("QFWNLLIYLMRVYLQKHA"), Sequence("EAKNILISLLIYWAQMLD"), Sequence("FMWNILVTIARVMAQQLD"), Sequence("TAWELLIKIARYMAQQLD"), Sequence("TMKEYLILALILYELQLS"), Sequence("PKREFLILALLIALKLES"), Sequence("KLTEIMLSIGLVFMWRKS"), Sequence("PEETFHRLLWEYMERLLA"), Sequence("EEEELWIQFLRLALKIAL"), Sequence("AYEMFQILFMWYLEMKDA"), Sequence("SYERMIELMLKWLEKHLA"), Sequence("LEYLLWILAMQYLEKHLA"), Sequence("TEREVTELLKIWRELFMA"), Sequence("ECRLLHILHIRYAKAWTA"), Sequence("PSKNIFLSLAWWIAQVLT")]
    pd1seq = "WNPPTFSPALLVVTEGDNATFTCSFSNTSESFHVVWHRESPSGQTDTLAAFPEDRSQPGQDSRFRVTQLPNGRDFHMSVVRARRNDSGTYVCGVISLAPKIQIKESLRAELRVTERRAE"
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    run_genetic_alg_pd1("/home/amrita/pd1/run_geneticalg/", pd1seq, startingseqs, pocketresidues, poolsize = 20, num_iter = 20)

