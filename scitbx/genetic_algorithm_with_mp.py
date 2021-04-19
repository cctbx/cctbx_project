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

##############################################################################
##################  GENETIC ALGORITHM WITH MULTIPROCESSING ###################
##############################################################################

class genetic_algorithm:
  def __init__(self,
      params,
      mutation_method,
      recombination_method,
      scoring_method,
      new_gene_method = None,
      log = sys.stdout):

    # Save inputs
    adopt_init_args(self,locals())

    self.cycle = 0
    self.run()

  def run(self):

    # Create genes using supplied new_gene_method
    self.generate_new_genes()

    n_cycles = self.get_number_of_cycles()
    print("\nRunning total of %s cycles" %(n_cycles), file = self.log)
    for cycle in range(n_cycles):
      print("Cycle %s (current top score: %s" %(
        cycle, self.get_best_score_as_text()), file = self.log)
      self.one_cycle()

    # All done, get best result
    best_gene = self.get_best_gene()
    print("Best gene with score of %s: \n%s" %(
      self.get_best_score_as_text(),self.get_best_gene()), file = self.log)

  def one_cycle(self):
    '''  Run one cycle of optimization '''

    # Mutate genes to create new genes
    self.create_new_genes_by_mutation_or_recombination(mutate = True)

    # Recombine genes to create new genes
    self.create_new_genes_by_mutation_or_recombination(recombine = True)

    # Select the best genes
    self.select_top_variants()

  def create_new_genes_by_mutation_or_recombination(self,
     mutate = None,
     recombine = None):

    all_new_genes = []

    nproc = self.params.nproc
    end_number = -1
    n_tot = len(self.genes)
    n = n_tot//nproc
    if n * nproc < n_tot:
      n = n + 1
    assert n * nproc >= n_tot

    runs_to_carry_out = []
    for run_id in range(nproc):
      start_number = end_number + 1
      end_number = min(n_tot-1, start_number + n -1)
      if end_number < start_number: continue
      runs_to_carry_out.append(group_args(
        run_id = run_id,
        start_number = start_number,
        end_number = end_number,
        ))

    kw_dict = {
      'params':self.params,
      'genes':self.genes,
      'mutate':mutate,
      'recombine':recombine,
      'mutation_method':self.mutation_method,
      'recombination_method':self.recombination_method,
      'scoring_method':self.scoring_method,
     }

    runs_carried_out = run_jobs_with_large_fixed_objects(
      nproc = nproc,
      verbose = False,
      kw_dict = kw_dict,
      run_info_list = runs_to_carry_out,
      job_to_run = group_of_mutate_or_recombine_one_gene,
      log = self.log)

    for run_info in runs_carried_out:
      new_genes = run_info.result.new_genes
      if new_genes:
        all_new_genes += new_genes


    next_gene_id = len(self.genes)
    for gene in all_new_genes:
      gene.gene_id = next_gene_id
      next_gene_id += 1

    self.genes += all_new_genes


  def get_best_score_as_text(self):
    score = self.get_best_score()
    if score is None:
      return "None"
    else:
      return "%.2f" %(score)

  def get_best_score(self):
    best_gene = self.get_best_gene()
    if not best_gene:
      return None
    else:
      return best_gene.score

  def get_best_gene(self):
    if not self.genes:
      return None
    self.sort_genes()
    return self.genes[0]

  def sort_genes(self):
    self.genes = sorted(self.genes,
      key = lambda gene: gene.score, reverse = True)

  def select_top_variants(self):
    if not self.genes: return
    self.sort_genes()

    n = self.get_number_of_variants_to_keep()
    print("Selecting top %s of %s genes (top score = %s )" %(
      n, len(self.genes), self.get_best_score_as_text()), file = self.log)
    self.genes = self.genes[:n]

  def score_genes(self):
    ''' Score all genes. If they already have a score...keep it'''
    for gene in self.genes:
      if gene.score is None:
         gene.score = self.scoring_method(self.params, gene)
    self.sort_genes()

  def generate_new_genes(self):
    self.genes = []
    for i in range(self.get_number_of_variants_to_make()):
      gene = self.new_gene_method(self.params)
      gene.gene_id = i
      self.genes.append(gene)
    print("Total of %s new genes created" %(len(self.genes)), file = self.log)

  def get_number_of_variants_to_keep(self):
    if self.params.number_of_variants is not None:
      return max(1,int(0.5+self.params.number_of_variants * \
       self.params.top_fraction_of_variants_to_keep))
    else:
      return max(1,int(0.5+self.params.typical_gene_length * \
         self.params.number_of_variants_per_gene_unit *\
       self.params.top_fraction_of_variants_to_keep))

  def get_number_of_variants_to_make(self):
    if self.params.number_of_variants is not None:
      return self.params.number_of_variants
    else:
      return self.params.typical_gene_length * \
         self.params.number_of_variants_per_gene_unit

  def get_number_of_cycles(self):
    if self.params.number_of_cycles is not None:
      return self.params.number_of_cycles
    else:
      return self.params.typical_gene_length * \
         self.params.number_of_cycles_per_gene_unit

def group_of_mutate_or_recombine_one_gene(
        run_info,
        params,
        genes = None,
	recombine = None,
        mutate = None,
        mutation_method = None,
        recombination_method = None,
        scoring_method = None,
        log = sys.stdout):
  new_genes = []
  for index in range(run_info.start_number, run_info.end_number + 1):
    info = mutate_or_recombine_one_gene(
      params,
      genes,
      index,
      recombine,
      mutate,
      mutation_method,
      recombination_method,
      scoring_method,
      log = log)
    if info and info.new_genes:
      new_genes += info.new_genes

  return group_args(
      group_args_type = ' one set of recombined/mutated genes',
      new_genes = new_genes,
    )

def mutate_or_recombine_one_gene(
        params,
        genes = None,
        index = None,
	recombine = None,
        mutate = None,
        mutation_method = None,
        recombination_method = None,
        scoring_method = None,
        log = sys.stdout):

    assert (mutate, recombine).count(True) == 1
    new_genes = []
    if mutate:
      # Mutate gene to create new gene
      gene = genes[index]
      new_gene = create_new_gene_by_mutation(
        params,
        gene,
        mutation_method,
        log = log)
      if new_gene:
        new_genes.append(new_gene)

    elif recombine:
      # Recombine gene to create new gene
      gene = genes[index]
      other_gene_index = random_integers(0,len(genes) - 2,1)[0]
      if other_gene_index >= index:
        other_gene_index += 1
      other_gene = genes[other_gene_index]
      new_gene = create_new_gene_by_recombination(
        params,
        gene,
        other_gene,
        recombination_method,
        log = log)
      if new_gene:
        new_genes.append(new_gene)

    for gene in new_genes:
      gene.score = scoring_method(params, gene)

    return group_args(
      group_args_type = ' one set of recombined/mutated genes',
      new_genes = new_genes,
    )


def create_new_gene_by_mutation(params,
    gene,
    mutation_method,
    log = sys.stdout):
  n = poisson(params.mutation_rate,1)[0]
  if n == 0:
    return # nothing to do
  new_gene = gene
  for i in range(n):
    new_gene = mutation_method(params, new_gene)
  return new_gene

def create_new_gene_by_recombination(params,
   gene,
   other_gene,
   recombination_method,
   log = sys.stdout):

  if uniform(0,1,1)[0] > params.recombination_rate:
    return # nothing to do
  return recombination_method(params, gene, other_gene)

##############################################################################
##################  END GENETIC ALGORITHM WITH MULTIPROCESSING ###############
##############################################################################



def tst_01():
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
  tst_01()

