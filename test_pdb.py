import pickle
from alphafold.common import protein
from pdb_parser import get_coordinates_pdb

filename = "../../test-af2/run2/outputs/sequences_0_model_1_multimer_0_results.pbz2.out"
f = pickle.load(open(filename, "rb"))
pdb = protein.to_pdb(f['unrelaxed_protein'])
c, r, ri = get_coordinates_pdb(pdb, fil = False)
print(c, r, ri)
