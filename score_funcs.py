from pdb_parser import get_coordinates_pdb
import math
import pickle
import numpy as np

def distance(p1, p2):
    """returns the distance between two 3D points represented as tuples"""

    dist = math.sqrt((float(p2[0])-float(p1[0]))**2+(float(p2[1])-(float(p1[1])))**2+(float(p2[2])-float(p1[2]))**2)
    return dist


def score_contacts(pdbfile, reslist1, reslist2):
    chains, residues = get_coordinates_pdb(pdbfile)
    score = 0
    pairs = []
    for res1 in reslist1:
        for res2 in reslist2:
            contact = 0
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=8:
                        contact=1
                        pair = (res1, res2)
                        if pair not in pairs:
                            pairs.append(pair)
            score = score+contact

    return pairs, score

def score_confidence(resultsfile, reslist1, reslist2, resindices):
    """calculates confidence score of pairwise residue interactions"""

    with open(resultsfile,'rb') as f:
        p = pickle.load(f)
        pae = p['pae_output']
        score = 0
        for res1 in reslist1:
            res1_id = resindices[res1]
            for res2 in reslist2:
                res2_id = resindices[res2]
                score = score + pae[res1_id][res2_id]
    return score


        print(len(pae))
        print(pae[1])


if __name__ == "__main__":
    score_confidence("/home/amrita/pd1/outputs/sequences_0_model_1_multimer_seed_0_results.pbz2.out",["A_VAL_35", "B_LEU_8"])
