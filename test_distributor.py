from distributor2 import Distributor
import random
import time


def shout_my_name(proc_id):
    random.seed(proc_id)
    print("initialization of process", proc_id)
    return sleepy_add_ten

def sleepy_add_ten(val):
    time.sleep(random.random())
    return val+10


if __name__ == "__main__":

    n_workers = 5

    dist = Distributor(n_workers, shout_my_name)
    work_list = [i for i in range(20)]
    results = dist.churn(work_list)
    dist.spin_down()
   
    for i, j in zip(work_list, results):
        print("result:", i, j)
