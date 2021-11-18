def get_coordinates_pdb(filename):
	lines = []
	chains = []
	residues = {}
	with open(filename,"r") as f:
		for lin in f:
			l = lin.strip().split()
			if 'ATOM' in l[0] or 'HETATM' in l[0]:
				resid = l[4]+'_'+l[3]+'_'+l[5]
				atominfo = (l[1], l[2], (l[6], l[7], l[8]))
				if l[4] not in chains:
					chains.append(l[4])
				if resid not in residues:
					residues[resid] = [atominfo]
				else:
					residues[resid].append(atominfo)
				# print(l)
	return chains, residues

if __name__ == "__main__":
	c, r = get_coordinates_pdb("/Users/amritanallathambi/Desktop/unc/kuhlmanlab/pd1/sequences_0_model_1_multimer_seed_0_unrelaxed.pdb")
	print(r)
