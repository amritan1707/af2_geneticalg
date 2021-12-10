from score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists
from pdb_parser import get_coordinates_pdb
import os
import subprocess
import shutil

def run_af2_pd1(pool, pd1seq, directory, fpath, flagsfile, af2path):
    """runs af2 on inputted list of Sequences and writes outputs to files"""
    os.makedirs(directory+"inputs")
    path = directory+"inputs/"
    with open(path+"sequences.csv", "w") as opf:
        for p in pool:
            opf.write(","+str(pd1seq)+","+str(p)+"\n")

    shutil.copy(fpath+flagsfile, directory)
    os.chdir(directory)
    return subprocess.run(["python", af2path+"run_af2.py", "@"+flagsfile])
    

def score_pd1(pdbfile, resultsfile, pocketresidues):
    residuelist = []
    chains, residues, resindices = get_coordinates_pdb(pdbfile)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdbfile, reslist1, reslist2)
    respairs1 = [x[0] for x in contacts]
    respairs2 = [x[1] for x in contacts]
    confidencescore = score_confidence_lists(resultsfile, reslist1, reslist2, resindices)
    return contacts, contactscore, confidencescore

if __name__ == "__main__":
    path = "/home/amrita/pd1/run1/outputs/"
    pdbfiles = []
    resultsfiles = []
    scores = []
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    for i in range(15):
        pdbfiles.append("sequences_"+str(i)+"_model_1_multimer_seed_0_unrelaxed.pdb")
        resultsfiles.append("sequences_"+str(i)+"_model_1_multimer_seed_0_results.pbz2.out")
    for pdb, res in zip(pdbfiles, resultsfiles):
        scores.append(score_pd1(path+pdb, path+res, pocketresidues)[1:])

    sorted_scores = sorted(scores, key=lambda x: (x[1], -x[0]))
    print(sorted_scores)
