# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:38:47 2020

@author: lixiangk
"""

from platypus import *
#from platypus.mpipool import MPIPool

import pickle
import numpy as np

from platypus.config import PlatypusConfig
default_evaluator = PlatypusConfig.default_evaluator

def initialize_algorithm(fullname):
    '''
    Initialize the algorithm from a saved population.

    Parameters
    ----------
    fullname : string
        The full name of a saved intermediate data, e.g., 'algorithm@0095.001.pkl'

    Returns
    -------
    algo : NSGAII or NSGAIII algorithm object
        DESCRIPTION.

    '''
    #fullname = 
    with open(fullname, 'rb') as f_handle:
        algo = pickle.load(f_handle)

    #print(algo.__dict__)
    return algo

class NSGAII(AbstractGeneticAlgorithm):
    
    def __init__(self, problem,
                 population_size = 100,
                 generator = RandomGenerator(),
                 selector = TournamentSelector(2),
                 variator = None,
                 archive = None,
                 nth_run = 1,
                 save_all = True,
                 **kwargs):
        '''
        New NSGAII class modified from the default one, enabling the saving of 
        intermediate results

        Parameters
        ----------
        nth_run : int, optional
            An integer number to indicating the nth run of the algorithm.
            The default is 1.
        save_all : bool, optional
            If saving all intermediate results. The default is True.
            
        Returns
        -------
        None.

        '''
        super(NSGAII, self).__init__(problem, population_size, generator, **kwargs)
        self.selector = selector
        self.variator = variator
        self.archive = archive
        self.nth_run = nth_run
        self.save_all = save_all
        
    def step(self):
        if self.nfe == 0:
            self.initialize()
        else:
            self.iterate()

        if self.archive is not None:
            self.result = self.archive
        else:
            self.result = self.population

        if self.save_all:
            self.dump()

    def initialize(self):
        super(NSGAII, self).initialize()
        
        if self.archive is not None:
            self.archive += self.population
        
        if self.variator is None:
            self.variator = default_variator(self.problem)

    def initialize_to_resume(self, generator = None):
        if generator is None:
            generator = self.generator
        
        self.population = [generator.generate(self.problem) for _ in range(self.population_size)]

        self.evaluate_all(self.population)
        
        if self.archive is not None:
            self.archive += self.population
        
        if self.variator is None:
            self.variator = default_variator(self.problem)

        if self.save_all:
            self.dump()

    def iterate(self):
        offspring = []
        
        while len(offspring) < self.population_size:
            parents = self.selector.select(self.variator.arity, self.population)
            offspring.extend(self.variator.evolve(parents))

        #retrieve the simulation manually
        #offspring = [self.generator.generate(self.problem) for _ in range(self.population_size)]
        
        self.evaluate_all(offspring)

        ### print/save the current offspring, added on September 23, 2019
        if self.save_all:
            current = []
            for k in np.arange(len(offspring)):
                r_vars = offspring[k].variables[:]
                r_objs = offspring[k].objectives[:]
                current.append(list(r_vars)+list(r_objs))
            with open('offspring@%04d.%03d' % (self.nfe, self.nth_run), 'w') as f_handle:
                print(self.nfe, np.atleast_2d(current)[-1])
                np.savetxt(f_handle, np.atleast_2d(current), fmt='%14.6f')
        ###
        
        offspring.extend(self.population)
        nondominated_sort(offspring)
        self.population = nondominated_truncate(offspring, self.population_size)

        if self.archive is not None:
            self.archive.extend(self.population)

    def dump(self):

        ### save the current population
        current = []
        for k in np.arange(len(self.population)):
            r_vars = self.population[k].variables[:]
            r_objs = self.population[k].objectives[:]
            current.append(list(r_vars)+list(r_objs))
        with open('population@%04d.%03d' % (self.nfe, self.nth_run), 'w') as f_handle:
            print(self.nfe, np.atleast_2d(current)[-1])
            np.savetxt(f_handle, np.atleast_2d(current), fmt='%14.6f')
        ###

        ### print/save the current result, added on September 13, 2019
        current = []
        for k in np.arange(len(self.result)):
            r_vars = self.result[k].variables[:]
            r_objs = self.result[k].objectives[:]
            current.append(list(r_vars)+list(r_objs))
        with open('result@%04d.%03d' % (self.nfe, self.nth_run), 'w') as f_handle:
            print(self.nfe, np.atleast_2d(current)[-1])
            np.savetxt(f_handle, np.atleast_2d(current), fmt='%14.6f')
        ###
        
        ### save the attributes except evaluator
        fullname = 'attribute@%04d.%03d.pkl' % (self.nfe, self.nth_run)
        att = self.__dict__.copy()
        att['evaluator'] = []
        #print('att is: ', att)
        with open(fullname, 'wb') as f_handler:
            pickle.dump(att, f_handler)
        ###
        
        ### save the algorithm w/o evaluator
        evaluator = self.evaluator
        self.evaluator = None
        
        fullname = 'algorithm@%04d.%03d.pkl' % (self.nfe, self.nth_run)
        with open(fullname, 'wb') as f_handler:
            pickle.dump(self, f_handler)
        
        self.evaluator = evaluator


class NSGAIII(AbstractGeneticAlgorithm):
    
    def __init__(self, problem,
                 divisions_outer,
                 divisions_inner = 0,
                 generator = RandomGenerator(),
                 selector = TournamentSelector(2),
                 variator = None,
                 nth_run = 1,
                 **kwargs):
        '''
        New NSGAIII class modified from the default one, enabling the saving of 
        intermediate results

        Parameters
        ----------
        nth_run : int, optional
            An integer number to indicating the nth run of the algorithm.
            The default is 1.
        save_all : bool, optional
            If saving all intermediate results. The default is True.
            
        Returns
        -------
        None.

        '''
        super(NSGAIII, self).__init__(problem, generator = generator, **kwargs)
        self.selector = selector
        self.variator = variator
        self.nth_run = nth_run
        
        self.population_size = choose(problem.nobjs + divisions_outer - 1, divisions_outer) + \
                (0 if divisions_inner == 0 else choose(problem.nobjs + divisions_inner - 1, divisions_inner))
        self.population_size = int(math.ceil(self.population_size / 4.0)) * 4

        print('population size: ', self.population_size)
        
        self.ideal_point = [POSITIVE_INFINITY]*problem.nobjs
        self.reference_points = normal_boundary_weights(problem.nobjs, divisions_outer, divisions_inner)
        
        # NSGAIII currently only works on minimization problems
        if any([d != Problem.MINIMIZE for d in problem.directions]):
            raise PlatypusError("NSGAIII currently only works with minimization problems")
        

    def step(self):
        if self.nfe == 0:
            self.initialize()
            self.result = self.population
        else:
            self.iterate()
            self.result = self.population
            
        self.dump() 
        
    def _find_extreme_points(self, solutions, objective):
        nobjs = self.problem.nobjs
        
        weights = [0.000001]*nobjs
        weights[objective] = 1.0
        
        min_index = -1
        min_value = POSITIVE_INFINITY
        
        for i in range(len(solutions)):
            objectives = solutions[i].normalized_objectives
            value = max([objectives[j]/weights[j] for j in range(nobjs)])
            
            if value < min_value:
                min_index = i
                min_value = value
                
        return solutions[min_index]

    def _associate_to_reference_point(self, solutions, reference_points):
        result = [[] for _ in range(len(reference_points))]
        
        for solution in solutions:
            min_index = -1
            min_distance = POSITIVE_INFINITY
            
            for i in range(len(reference_points)):
                distance = point_line_dist(solution.normalized_objectives, reference_points[i])
        
                if distance < min_distance:
                    min_index = i
                    min_distance = distance
                    
            result[min_index].append(solution)
            
        return result
    
    def _find_minimum_distance(self, solutions, reference_point):
        min_index = -1
        min_distance = POSITIVE_INFINITY
            
        for i in range(len(solutions)):
            solution = solutions[i]
            distance = point_line_dist(solution.normalized_objectives, reference_point)
        
            if distance < min_distance:
                min_index = i
                min_distance = distance
                    
        return solutions[min_index]
        
    def _reference_point_truncate(self, solutions, size):
        nobjs = self.problem.nobjs
        
        if len(solutions) > size:
            result, remaining = nondominated_split(solutions, size)

            # update the ideal point
            for solution in solutions:
                for i in range(nobjs):
                    self.ideal_point[i] = min(self.ideal_point[i], solution.objectives[i])
                    
            # translate points by ideal point
            for solution in solutions:
                solution.normalized_objectives = [solution.objectives[i] - self.ideal_point[i] for i in range(nobjs)]
            
            # find the extreme points
            extreme_points = [self._find_extreme_points(solutions, i) for i in range(nobjs)]
            
            # calculate the intercepts
            degenerate = False
            
            try:
                b = [1.0]*nobjs
                A = [s.normalized_objectives[:] for s in extreme_points]
                x = lsolve(A, b)
                intercepts = [1.0 / i for i in x]
            except:
                degenerate = True
                
            if not degenerate:
                for i in range(nobjs):
                    if intercepts[i] < 0.001:
                        degenerate = True
                        break
                    
            if degenerate:
                intercepts = [-POSITIVE_INFINITY]*nobjs
                
                for i in range(nobjs):
                    intercepts[i] = max([s.normalized_objectives[i] for s in solutions] + [EPSILON])
    
            # normalize objectives using intercepts
            for solution in solutions:
                solution.normalized_objectives = [solution.normalized_objectives[i] / intercepts[i] for i in range(nobjs)]
    
            # associate each solution to a reference point
            members = self._associate_to_reference_point(result, self.reference_points)
            potential_members = self._associate_to_reference_point(remaining, self.reference_points)
            excluded = set()
            
            while len(result) < size:
                # identify reference point with the fewest associated members
                min_indices = []
                min_count = sys.maxsize
                
                for i in range(len(members)):
                    if i not in excluded and len(members[i]) <= min_count:
                        if len(members[i]) < min_count:
                            min_indices = []
                            min_count = len(members[i])
                        min_indices.append(i)
                
                # pick one randomly if there are multiple options
                min_index = random.choice(min_indices)
                
                # add associated solution
                if min_count == 0:
                    if len(potential_members[min_index]) == 0:
                        excluded.add(min_index)
                    else:
                        min_solution = self._find_minimum_distance(potential_members[min_index], self.reference_points[min_index])
                        result.append(min_solution)
                        members[min_index].append(min_solution)
                        potential_members[min_index].remove(min_solution)
                else:
                    if len(potential_members[min_index]) == 0:
                        excluded.add(min_index)
                    else:
                        rand_solution = random.choice(potential_members[min_index])
                        result.append(rand_solution)
                        members[min_index].append(rand_solution)
                        potential_members[min_index].remove(rand_solution)
                        
            return result
        else:
            return solutions
        
    def initialize(self):
        super(NSGAIII, self).initialize()
        
        if self.variator is None:
            self.variator = default_variator(self.problem)
    
    def iterate(self):
        offspring = []
        
        while len(offspring) < self.population_size:
            parents = self.selector.select(self.variator.arity, self.population)
            offspring.extend(self.variator.evolve(parents))
            
        self.evaluate_all(offspring)
        
        ### print/save the current offspring, added on September 23, 2019
        current = []
        for k in np.arange(len(offspring)):
            r_vars = offspring[k].variables[:]
            r_objs = offspring[k].objectives[:]
            current.append(list(r_vars)+list(r_objs))
        with open('offspring@%04d.%03d' % (self.nfe, self.nth_run), 'w') as f_handle:
            print(self.nfe, np.atleast_2d(current)[-1])
            np.savetxt(f_handle, np.atleast_2d(current), fmt='%14.6f')
        ###
        
        offspring.extend(self.population)
        nondominated_sort(offspring)
        self.population = self._reference_point_truncate(offspring, self.population_size)

    def dump(self):

        ### save the current population
        current = []
        for k in np.arange(len(self.population)):
            r_vars = self.population[k].variables[:]
            r_objs = self.population[k].objectives[:]
            current.append(list(r_vars)+list(r_objs))
        with open('population@%04d.%03d' % (self.nfe, self.nth_run), 'w') as f_handle:
            print(self.nfe, np.atleast_2d(current)[-1])
            np.savetxt(f_handle, np.atleast_2d(current), fmt='%14.6f')
        ###

        ### print/save the current result, added on September 13, 2019
        current = []
        for k in np.arange(len(self.result)):
            r_vars = self.result[k].variables[:]
            r_objs = self.result[k].objectives[:]
            current.append(list(r_vars)+list(r_objs))
        with open('result@%04d.%03d' % (self.nfe, self.nth_run), 'w') as f_handle:
            print(self.nfe, np.atleast_2d(current)[-1])
            np.savetxt(f_handle, np.atleast_2d(current), fmt='%14.6f')
        ###
        
        ### save the attributes except evaluator
        fullname = 'attribute@%04d.%03d.pkl' % (self.nfe, self.nth_run)
        att = self.__dict__.copy()
        att['evaluator'] = []
        #print('att is: ', att)
        with open(fullname, 'wb') as f_handler:
            pickle.dump(att, f_handler)
        ###
        
        ### save the algorithm w/o evaluator
        evaluator = self.evaluator
        self.evaluator = None
        
        fullname = 'algorithm@%04d.%03d.pkl' % (self.nfe, self.nth_run)
        with open(fullname, 'wb') as f_handler:
            pickle.dump(self, f_handler)
        
        ### restore the evaluator
        self.evaluator = evaluator
