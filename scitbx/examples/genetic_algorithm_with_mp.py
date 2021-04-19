# -*- coding: utf-8 -*-

from __future__ import division, print_function
from libtbx import adopt_init_args
from libtbx import group_args
import sys
from numpy.random import poisson, uniform, random_integers
from libtbx.easy_mp import run_jobs_with_large_fixed_objects

##############################################################################
#############  INSTRUCTIONS FOR GENETIC ALGORITHM WITH MULTIPROCESSING #######
##############################################################################
'''
  Simple generic implementation of genetic algorithm with multiprocessing

  To use, you need to create the following methods (You can copy examples
    from below to get started):

    new_gene_method(params)  #  returns a new gene
    mutation_method(params, gene) # returns a mutated gene
    recombination_method(params, gene, other_gene)  # returns recombined gene
    scoring_method(params, gene)  # returns a score

  Here a gene looks like (see example methods below):

  gene = group_args(
     group_args_type = 'default gene class as group_args',
     gene_id = gene_id,
     length = length,
     values = values,
     score = score)

  Then you can edit the example params below and run with:

  ga = genetic_algorithm(
    params = params,
    new_gene_method = new_gene_method,
    mutation_method = mutation_method,
    recombination_method = recombination_method,
    scoring_method = scoring_method,
   )

  and your result will be available with:

  ga.get_best_result()
'''
##############################################################################
##########  EXAMPLE PARAMS FOR GENETIC ALGORITHM WITH MULTIPROCESSING #######
##############################################################################

''' Default parameters for genetic algorithm'''

example_params = group_args(
    group_args_type = 'parameters for genetic algorithm',
    nproc = 1,
    mutation_rate = 0.5,
    recombination_rate = 0.5,
    number_of_variants = None,
    number_of_cycles  = None,
    top_fraction_of_variants_to_keep = 0.2,
    number_of_variants_per_gene_unit = 10,
    number_of_cycles_per_gene_unit = 100,
    typical_gene_length = 10,
  )
##############################################################################
##########  EXAMPLE METHODS FOR GENETIC ALGORITHM WITH MULTIPROCESSING #######
##############################################################################

def example_new_gene_method(params):
  '''
  Default gene for genetic algorithm
  You can use anything but it needs to have values for:
    length
    values
    score
    gene_id
   '''
  from scitbx.array_family import flex
  gene_id = 0
  length = params.typical_gene_length
  values = flex.random_double(length)
  score = None

  return group_args(
     group_args_type = 'default gene class as group_args',
     gene_id = gene_id,
     length = length,
     values = values,
     score = score)

def example_mutation_method(params, gene):
  ''' Must take params and a gene and return a new gene '''
  from copy import deepcopy
  new_gene = deepcopy(gene)
  new_gene.values[random_integers(0,gene.length-1,1)[0]] += uniform(-4,4,1)[0]
  return new_gene

def example_recombination_method(params,gene_1, gene_2):
  ''' Must take params and two genes and return a new gene '''
  from copy import deepcopy
  new_gene = deepcopy(gene_1)
  i = random_integers(0,gene_1.length-1,1)[0]
  k = random_integers(0,1,1)[0]
  if k == 0:
    new_gene.values = gene_1.values[:i]
    new_gene.values.extend(gene_2.values[i:])
  else:
    new_gene.values = gene_2.values[:i]
    new_gene.values.extend(gene_1.values[i:])
  return new_gene

def example_scoring_method(params, gene):
  ''' Must take params and a gene and return a number'''
  from scitbx.array_family import flex
  perfect = flex.double(list(range(gene.length)))  #
  return params.typical_gene_length - (gene.values - perfect).rms()


def example():
  '''  Exercise genetic algorithm '''

  from scitbx.genetic_algorithm_with_mp import  genetic_algorithm,\
    example_params, example_new_gene_method, example_mutation_method, \
    example_recombination_method, example_scoring_method

  params = example_params
  new_gene_method = example_new_gene_method
  mutation_method = example_mutation_method
  recombination_method = example_recombination_method
  scoring_method = example_scoring_method

  result = genetic_algorithm(
    params = params,
    new_gene_method = new_gene_method,
    mutation_method = mutation_method,
    recombination_method = recombination_method,
    scoring_method = scoring_method,
   )
  best_gene = result.get_best_gene()
  print (best_gene.score)
  print (list(best_gene.values))
  assert abs(best_gene.score  - 10 ) < 0.1

if __name__ == "__main__":
  example()

