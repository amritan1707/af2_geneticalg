import multiprocessing as mp
import collections
import random
import time
import sys, os
from typing import Sequence, Union

sys.path.append('/proj/kuhl_lab/alphafold/run')

from run_af2 import af2
from functools import partial
import numpy as np

class Distributor:
    """This class will distribute work to sub-processes where
    the same function is run repeatedly with different inputs.
    The Distributor is initialized with the number of worker
    processes and a function, f_init, whose job it is to create
    a function "f" that will be run repeatedly by the worker
    processes.
    """

   
    def __init__(self, n_workers, f_init, arg_file, lengths, fitness_fxn):
        """
        Construct a Distributor that manages n_workers sub-processes.
        The distributor will give work to its sub-processes in the
        form of function inputs.

        f_init should be a function and return a function. It will
        be called by each sub-process once as that process gets started.
        It should return the worker function that will do the heavy
        lifting.
        """
        self.n_workers = n_workers
        self.lock = mp.Lock()
        self.qs_out = [mp.Queue() for _ in range(n_workers)]
        self.q_in = mp.Queue()
        self.processes = [
            mp.Process(
                target=Distributor._worker_loop,
                args=(
                    f_init,
                    i,
                    self.lock,
                    self.qs_out[i],
                    self.q_in,
		    arg_file,
		    lengths,
                    fitness_fxn)
            )
            for i in range(n_workers)
        ]
   
        for p in self.processes:
            p.start()


    def spin_down(self):
        """When all work is done, send out a spin-down signal to all the
        subprocesses and join them
        """
        for i in range(self.n_workers):        
            # spin down the worker
            self.qs_out[i].put((False, None))
        for i in range(self.n_workers):
            self.processes[i].join()
           

    def churn(self, work_list):
        """Process the work in the work list, farming out work
        to the subprocesses
        """
        print("Distributor::churn", len(work_list))
        work_queue = collections.deque(work_list)
        n_jobs = len(work_queue)
        job_ind_for_worker = [-1] * self.n_workers
        job_output = [None] * n_jobs


        count_job = 0
        count_completed = 0
       
        # step 1: put all the work we already have into the work queues
        if self.n_workers < n_jobs:
            for i in range(self.n_workers):
                w = work_queue.popleft()
                self.qs_out[i].put((True, w))
                job_ind_for_worker[i] = count_job
                count_job += 1
        else:
            count = 0
            while len(work_queue) > 0:
                w = work_queue.popleft()
                self.qs_out[count].put((True, w))
                job_ind_for_worker[count] = count_job
                count += 1
                count_job += 1

        while len(work_queue) > 0:
            proc_id, val = self.q_in.get()
            count_completed += 1
            #print("work_list loop count completed:", count_completed)
            job_ind = job_ind_for_worker[proc_id]
            assert job_ind != -1
            job_output[job_ind] = val
   
            w = work_queue.popleft()
            self.qs_out[proc_id].put((True, w))
            job_ind_for_worker[proc_id] = count_job
            count_job += 1

        while count_completed < n_jobs:
            proc_id, val = self.q_in.get()
            #print("wait-for-jobs to finish loop. count completed:", count_completed)
            count_completed += 1
            job_ind = job_ind_for_worker[proc_id]
            assert job_ind != -1
            job_ind_for_worker[proc_id] = -1
            job_output[job_ind] = val

        return job_output

           
    @staticmethod
    def _worker_loop(f_init, proc_id, lock, q_in, q_out, arg_file, lengths, fitness_fxn):

        f = f_init(proc_id, arg_file, lengths, fitness_fxn)
   
        is_job, val = q_in.get()
        while is_job:
            result = f(val)
       
            lock.acquire()
            try:
                q_out.put((proc_id, result))
            finally:
                lock.release()
            is_job, val = q_in.get()
        print("spinning down worker", proc_id)


def shout_my_name(proc_id):
    random.seed(proc_id)
    print("initialization of process", proc_id)
    return sleepy_add_ten


def af2_init(proc_id: int, arg_file: str, lengths: Sequence[Union[str, Sequence[str]]], fitness_fxn):
    print('initialization of process', proc_id)

    os.environ['TF_FORCE_UNITED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '2.0'
    os.environ['TF_XLA_FLAGS'] = '--tf_xla_cpu_global_jit'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(proc_id)

    import jax
    from features import (getRawInputs, getChainFeatures, getInputFeatures)
    from setup import (getAF2Parser, QueryManager, getOutputDir)
    from model import (getModelNames, getModelRunner, predictStructure, getRandomSeeds)

    parser = getAF2Parser()
    args = parser.parse_args([f'@{arg_file}'])    

    output_dir = getOutputDir(out_dir=args.output_dir)

    print(f'lengths: {lengths}')
    sequences = []
    for protein in lengths:
        if type(protein) == int:
            sequences.append( ['A'*protein] )
        elif type(protein) == list:
            sequence = []
            for chain in protein:
                sequence.append( 'A'*chain )
            sequences.append(sequence)
    print(f'determined sequences: {sequences}')

    qm = QueryManager(
	input_dir=args.input_dir,
	sequences=sequences,
	min_length=args.min_length,
	max_length=args.max_length,
	max_multimer_length=args.max_multimer_length)
    qm.parse_files()
    qm.parse_sequences()
    queries = qm.queries
    del qm

    raw_inputs = getRawInputs(
        queries=queries,
        msa_mode=args.msa_mode,
        use_templates=args.use_templates,
        output_dir=output_dir)

    print('finished generating raw inputs')
    print(raw_inputs)

    model_names = getModelNames(
        first_n_seqs=len(queries[0][1]),
        last_n_seqs=len(queries[-1][1]),
        use_ptm=args.use_ptm, num_models=args.num_models)

    query_features = []
    for query in queries:
        sequences = query[1]
    
        features_for_chain = getChainFeatures(
            sequences=sequences,
            raw_inputs=raw_inputs,
            use_templates=args.use_templates)

        input_features = getInputFeatures(
            sequences=sequences,
            chain_features=features_for_chain)

        query_features.append( (sequences, input_features) )

    results_list = []
    for model_name in model_names:
        model_runner = getModelRunner(
            model_name=model_name,
            num_ensemble=args.num_ensemble,
            is_training=args.is_training,
            num_recycle=args.max_recycle,
            recycle_tol=args.recycle_tol,
            params_dir=args.params_dir)

        run_multimer = False
        if 'multimer' in model_name:
            run_multimer = True
        
        for query in query_features:
            sequences = query[0]

            if len(sequences) > 1 and not run_multimer:
                continue
            elif len(sequences) == 1 and run_multimer:
                continue

            input_features = query[1]      
            
            t = time.time()
            result = predictStructure(
                model_runner=model_runner,
                feature_dict=input_features,
                run_multimer=run_multimer)
            print(f'Model {model_name} took {time.time()-t} sec on GPU {proc_id}.')

    af2_partial = partial(af2, arg_file=arg_file, proc_id=proc_id, fitness_fxn=fitness_fxn, compiled_runner=model_runner)

    return af2_partial


def naive_fitness(result):
    return np.mean(result['plddt'])


def generate_random_monomers(length: int, num_seq: int, aalist=None) -> Sequence[Sequence[str]]:
    if aalist == None:
        aalist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    return [[''.join(random.choices(aalist, k=length))] for _ in range(num_seq)]


def generate_random_multimers(lengths: Sequence[int], num_seq: int, aalist=None) -> Sequence[Sequence[Sequence[str]]]:
    if aalist == None:
        aalist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    seqs_list = []
    for _ in range(num_seq):
        seqs = []
        for length in lengths:
            seq = ''.join(random.choices(aalist, k=length))
            seqs.append(seq)

        seqs_list.append([seqs]) #seqs_list.append(seqs)

    return seqs_list #return [seqs_list]


if __name__ == "__main__":

    n_workers = 2
    init_len = 50

    lengths = [[init_len, init_len]]
    dist = Distributor(n_workers, af2_init, 'flags.txt', lengths, naive_fitness)
    
    all_work = []
    all_results = []
    for _ in range(30):
        work_list = generate_random_multimers(lengths[0], 4)
        print(work_list)
        results = dist.churn(work_list)
        for seq in work_list:
            all_work.append(seq)
        for score in results:
            all_results.append(score)

    dist.spin_down()
   
    for i, j in zip(all_work, all_results):
        print("result:", i[0], j[0])
