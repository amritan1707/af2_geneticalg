from score_funcs import score_contacts, score_confidence
from pdb_parser import get_coordinates_pdb

def score_pd1(pdbfile, resultsfile):
    residuelist = []
    chains, residues, resindices = get_coordinates_pdb(pdbfile)
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdbfile, reslist1, reslist2)
    respairs1 = [x[0] for x in contacts]
    respairs2 = [x[1] for x in contacts]
    confidencescore = score_confidence(resultsfile, respairs1, respairs2, resindices)
    return contacts, contactscore, confidencescore


if __name__ == "__main__":
    print(score_pd1("/home/amrita/pd1/outputs/sequences_0_model_1_multimer_seed_0_unrelaxed.pdb", "/home/amrita/pd1/outputs/sequences_0_model_1_multimer_seed_0_results.pbz2.out"))
