import random
from Sequence import Sequence
from score_pd1 import score_pd1, run_af2_pd1

def create_new_seqs(startseqs, num_seqs, crossoverpercent = 0.2):
    """takes starting sequences and creates a pool of size num_seqs by mutation and crossover"""
    pool = startseqs
    while len(pool)<num_seqs*(1-crossoverpercent):
        obj = random.choice(startseqs)
        newseq = obj.mutate(resnums = [0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16], restypes = ["alpha" for x in range(13)])
        newseq = Sequence("".join(newseq))
        if newseq not in pool:
            pool.append(newseq)
    #print(pool)
    while len(pool)<num_seqs:
        oldseqs = random.sample(startseqs, 2)
        newseq = oldseqs[0].crossover(oldseqs[1])
        newseq = Sequence("".join(newseq))
        if newseq not in pool:
            pool.append(newseq)

    return pool

def run_genetic_alg_pd1(path, pd1seq, startingseqs, poolsize = 50, num_iter = 20):
    pool = create_new_seqs(startingseqs, poolsize)
    scored_pool = {}
    curr_iter = 1
    while curr_iter <= num_iter:

        #uncomment when running
        #run_af2_pd1(pool, pd1seq, path) 

        oppath = path+"outputs/"

        #iterate through outputs and read
        #add score for each sequence (not already in scored_pool) into scored_pool
        #contacts, contactscores, confscores = score_pd1(oppath)
        #sort by score
        #get rid of lower scoring half

        curr_iter+=1

    for p in pool:
        print(str(p))
        #scored_pool[p] = score_pd1(p)
    
    #print(scored_pool)

if __name__=="__main__":
    startingseqs = [Sequence("PSREFLILALQIALTLKA"), Sequence("QFWNLLIYLMRVYLQKHA"), Sequence("EAKNILISLLIYWAQMLD"), Sequence("FMWNILVTIARVMAQQLD"), Sequence("TAWELLIKIARYMAQQLD"), Sequence("TMKEYLILALILYELQLS"), Sequence("PKREFLILALLIALKLES"), Sequence("KLTEIMLSIGLVFMWRKS"), Sequence("PEETFHRLLWEYMERLLA"), Sequence("EEEELWIQFLRLALKIAL"), Sequence("AYEMFQILFMWYLEMKDA"), Sequence("SYERMIELMLKWLEKHLA"), Sequence("LEYLLWILAMQYLEKHLA"), Sequence("TEREVTELLKIWRELFMA"), Sequence("ECRLLHILHIRYAKAWTA"), Sequence("PSKNIFLSLAWWIAQVLT")]
    pd1seq = "WNPPTFSPALLVVTEGDNATFTCSFSNTSESFHVVWHRESPSGQTDTLAAFPEDRSQPGQDSRFRVTQLPNGRDFHMSVVRARRNDSGTYVCGVISLAPKIQIKESLRAELRVTERRAE"
    run_genetic_alg_pd1("/home/amrita/pd1/run_geneticalg/test", pd1seq, startingseqs, poolsize = 20, num_iter = 20)

