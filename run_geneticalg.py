import random

def run_af2_pd1():
    """runs af2 on inputted list of Sequences and writes outputs to files"""

def create_new_seqs(startseqs, num_seqs, crossoverpercent = 0.1):
    """takes starting sequences and creates a pool of size num_seqs by mutation and crossover"""
    pool = startseqs
    while len(pool)<num_seqs:
        seq = pool[random.randrange(len(pool))]
        newseq = seq.mutate()
        pool.append(newseq)


if __name__=="__main__":
    startingseqs = []
    pool = create_new_seqs(startingseqs, 100)
    

