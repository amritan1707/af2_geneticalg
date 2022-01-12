import random
import time
import sys
import numpy as np

sys.path.append('/proj/kuhl_lab/alphafold/run')

from distributor import Distributor
from run_af2 import af2_init

def generate_random_seqs(length, num_seq, aalist=None):
    if aalist == None:
        aalist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
                  'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
    return [[''.join(random.choices(aalist, k=length))] for _ in range(num_seq)]


def naive_scoring(result_dict):
    return np.mean(result_dict['plddt'])


if __name__ == "__main__":

    n_workers = 2

    dist = Distributor(n_workers, af2_init, init_lens=[[20]], arg_file='flags.txt', fitness_fxn=naive_scoring)

    for i in range(2):
        work_list = generate_random_seqs(length=20, num_seq=4)
        results = dist.churn(work_list)
        print(f'Work list {i+1} has been completed.')
    
    dist.spin_down()
