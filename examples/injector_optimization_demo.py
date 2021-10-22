# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 21:30:54 2020

@author: lixiangk
"""

# Import the goal function
from  injector_demo import *

import timeit

### Use the default version of NSGAII or NSAGIII
#from platypus import *
###

### or use the modified version that can save intermediate results
from interface.NSGAPlus import *
###

from platypus.mpipool import MPIPool
import logging

logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    # define the problem definition
    problem = Problem(4, 3)
    problem.types[0] = Real(2.4, 4.8)
    problem.types[1] = Real(-15, 15)
    problem.types[2] = Real(-30, -5)
    problem.types[3] = Real(340, 380)
    problem.function = obj_THz4nC_MOGA
    
    time1 = timeit.default_timer()
    
    pool = MPIPool()
    print(pool)
    # only run the algorithm on the master process
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
        
    # instantiate the optimization algorithm to run in parallel
    with PoolEvaluator(pool) as evaluator:
        # Note population_size = number_of_cores-1
        ### Either from zero
        algorithm = NSGAII(problem, 
                           evaluator=evaluator, 
                           population_size = 15)
        
        # With NSGAIII
        #algorithm = NSGAIII(problem,
        #                    divisions_outer = 12,
        #                    evaluator=evaluator,
        #                    population_size = 15)
        ###
        
        ### Or from an intermediate result
        #algorithm = initialize_algorithm('algorithm@0031.001.pkl')
        #algorithm.evaluator = evaluator
        #algorithm.nth_run = 2
        ###
        
        
        algorithm.run(30)
    
    # display the results
    for solution in algorithm.result:
        print(solution.objectives)

    pool.close()
    
    time2 = timeit.default_timer()
    print('time elapsed: ', time2-time1)