import random

class Sequence:

    def __init__(self, seq):
        self.seq = seq
        #self.score = score(self, type="PD1")

    def __eq__(self, other):
        return self.seq == other.seq

    def __str__(self):
        return self.seq

    def __hash__(self):
        return hash(self.seq)

    def mutate(self, resnums = None, mutpercent = 0.125, restypes = None, use_blosum = True):
        if resnums is None:
            resnums = [x for x in range(len(str(self)))]
        if restypes is None:
            restypes = ["all" for x in range(len(str(self)))]

        if not len(resnums)==len(restypes):
            print("Length of resnums and restypes lists need to be the same")
            return

        if max(resnums)>len(str(self)) or min(resnums)<0:
            print("invalid resnum")
            return

        #calculating number of mutants from mutpercent and maybe adding or subtracting one mutant for stochasticity
        num_mut = round(len(str(self))*mutpercent) + random.choice([-1, 0, 1])

        similarities = {}
        if use_blosum:
            matrix = []
            with open("blosum62.txt","r") as blosumf:
                for lin in blosumf:
                    l = lin.strip().split(" ")
                    l = [x for x in l if x]
                    matrix.append(l)
            aas = matrix[0]
            for m in matrix[1:]:
                aa = m[0]
                temp = {}
                vals = [float(x) for x in m[1:]]
                for weight, corraa in zip(vals, aas):
                    norm_weight = (float(weight)-min(vals))/(max(vals)-min(vals))
                    temp[corraa] = norm_weight

                similarities[aa]=temp

        newseq = list(str(self))
        oldseq = list(str(self))
        hydrophobic = ["G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]
        hydphob_weights = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        hydrophilic = ["Q", "D", "K", "R", "E", "N", "S", "T"]
        hydphil_weights = [1, 1, 1, 1, 1, 1, 1, 1]
        alpha = ["A", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "Q", "R", "S", "T", "V", "W", "Y"]
        alpha_weights = [1, 0.31, 0.6, 0.46, 0, 0.39, 0.59, 0.74, 0.79, 0.74, 0.35, 0.61, 0.79, 0.5, 0.34, 0.39, 0.51, 0.47]
        allres = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        all_weights = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        resnums_types = {}
        for resnum, restype in zip(resnums,restypes):
            resnums_types[resnum] = restype

        mut_ids = random.sample(resnums, num_mut)
        #print(mut_ids)
        for mut_id in mut_ids:
            w = []
            res = []
            if resnums_types[mut_id] is "alpha":
                w = alpha_weights
                res = alpha
            if resnums_types[mut_id] is "hydphob":
                w = hydphob_weights
                res = hydrophobic
            if resnums_types[mut_id] is "hydphil":
                w = hydphil_weights
                res = hydrophilic
            if resnums_types[mut_id] is "all":
                w = all_weights
                res = allres

            weights = []
            if use_blosum:
                for aa, aw in zip(res, w):
                    sw = similarities[newseq[mut_id]][aa]
                    weights.append(aw*sw)
            else:
                weights = w
            newseq[mut_id] = random.choices(alpha, weights)[0]
            while newseq[mut_id]==oldseq[mut_id]:
                newseq[mut_id] = random.choices(alpha, weights)[0]
        return newseq 
    
    def crossover(self, otherSeq, ncross = 1, p=None):
        seq1 = list(str(self))
        seq2 = list(str(otherSeq))
        newseq = list(seq1)

        if p:
            if p>=0 and p<=1:
                s=0
                for i, aa1, aa2 in zip(range(len(seq1)), seq1, seq2):
                    options = [aa1, aa2]
                    if p > random.random():
                        if s==0:
                            s=1
                        else:
                            s=0
                    newseq[i] = options[s]
            else:
                print("invalid p value")
                return
        else:
            if ncross<len(seq1):
                points = random.sample(range(len(seq1)), ncross)
                s=0
                for i, aa1, aa2 in zip(range(len(seq1)), seq1, seq2):
                    options = [aa1, aa2]
                    if i in points:
                        if s==0:
                            s=1
                        else:
                            s=0
                    newseq[i] = options[s]
            else:
                print("invalid ncross value")
                return
        #print(newseq)
        return newseq

        
if __name__ == "__main__":
    test1 = Sequence("PSREFLILALQIALTLKA")
    test2 = Sequence("AAAAAAAAAAAAAAAAAA")
    test3 = Sequence("BBBBBBBBBBBBBBBBBB")
    resnums_test = [0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16]
    restypes_test = ["alpha" for x in resnums_test]
    print(str(test2))
    print("".join(test2.mutate(resnums = resnums_test, restypes = restypes_test)))
    print(test2.crossover(test3, ncross = 1))

