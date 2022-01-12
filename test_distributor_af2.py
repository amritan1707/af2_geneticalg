import random
import time
from score_pd1_gpus import score_pd1_func
import sys

def generate_random_seqs(aalist, num_aas, num_seqs):
    oplist = []
    for i in range(num_seqs):
        oplist.append("".join(random.choices(aalist, k=num_aas)))

    return oplist

def shout_my_name(proc_id):
    random.seed(proc_id)
    print("initialization of process", proc_id)
    return run_af2_and_score

def run_af2_and_score(seq):
    results = run_af2_gpu(seq)
    return score_pd1_func(results)

def naive_fitness(result):
    return result['pae_output']

if __name__ == "__main__":

    n_workers = 2
    len_seq = 16
    num_seqs = 16
    pd1seq = "WNPPTFSPALLVVTEGDNATFTCSFSNTSESFHVVWHRESPSGQTDTLAAFPEDRSQPGQDSRFRVTQLPNGRDFHMSVVRARRNDSGTYVCGVISLAPKIQIKESLRAELRVTERRAE"

    aalist = ["A", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "Q", "R", "S", "T", "V", "W", "Y"]
    sequences = generate_random_seqs(aalist, len_seq, num_seqs)
    #sys.path.append("/proj/kuhl_lab/alphafold/run/")
    sys.path.append("/nas/longleaf/home/nzrandol/for_amrita/")
    #from run_af2 import af2_init
    from af_design import Distributor
    from af_design import af2_init
    dist = Distributor(n_workers, af2_init, "flags.txt", [[len(pd1seq), len_seq]], naive_fitness)
    work_list = [[[pd1seq, seq]] for seq in sequences]
    print(work_list)
    results = dist.churn(work_list)
    dist.spin_down()
   
    for i, j in zip(work_list, results):
        print("result:", i, j)
