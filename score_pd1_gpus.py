from score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists
from pdb_parser import get_coordinates_pdb
import os
import subprocess
from alphafold.common import protein
import shutil

def score_pd1_func(results):
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    print(resindices)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    #respairs1 = [x[0] for x in contacts]
    #respairs2 = [x[1] for x in contacts]
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    return contacts, contactscore, confidencescore

if __name__ == "__main__":
    score_pd1_func()
