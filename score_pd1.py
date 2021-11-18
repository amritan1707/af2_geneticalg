from score_funcs import score_contacts, score_confidence

def score_pd1(pdbfile, resultsfile, reslist1, reslist2, resindices):
    residuelist = []
    contacts, contactscore = score_contacts(pdbfile,["A_VAL_35"], ["B_LEU_8"])
    print(contacts, contactscore)
    confidencescore = score_confidence(resultsfile, reslist1, reslist2, resindices)
    return contacts, contactscore+confidencescore

if __name__ == "__main__":
    score_pd1("/home/amrita/pd1/outputs/sequences_0_model_1_multimer_seed_0_unrelaxed.pdb", "/home/amrita/pd1/outputs/sequences_0_model_1_multimer_seed_0_results.pbz2.out")