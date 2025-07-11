"""Simple generic implementation of genetic algorithm with multiprocessing"""

# -*- coding: utf-8 -*-

from __future__ import division, print_function
from libtbx import adopt_init_args
from libtbx import group_args
from libtbx.utils import null_out
import sys
import numpy.random as np_random
from libtbx.easy_mp import run_jobs_with_large_fixed_objects

##############################################################################
#############  INSTRUCTIONS FOR GENETIC ALGORITHM WITH MULTIPROCESSING #######
##############################################################################
'''
  Simple generic implementation of genetic algorithm with multiprocessing

  To use, you need to create the following methods (You can copy examples
    from below to get started):

    new_gene_method(params, n)  #  returns n new genes
    mutation_method(params, gene) # returns a mutated gene
    recombination_method(params, gene, other_gene)  # returns recombined gene
    scoring_method(params, gene)  # returns a score
    genes_are_identical_method(gene, other_gene)  # returns True if same

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
    genes_are_identical_method = genes_are_identical_method,
   )

  and your result will be available with:

  ga.get_best_gene()

  NOTE: you may need to import some functions...see the top of this file
  Example:
    import numpy.random as np_random

'''
##############################################################################
##########  EXAMPLE PARAMS FOR GENETIC ALGORITHM WITH MULTIPROCESSING #######
##############################################################################

''' Default parameters for genetic algorithm'''

example_params = group_args(
    group_args_type = 'parameters for genetic algorithm',
    nproc = 1,
    random_seed = 784321,
    mutation_rate = 0.5,
    recombination_rate = 0.5,
    number_of_variants = None,
    total_number_of_cycles = 1000,
    number_of_cycles = None,
    number_of_macro_cycles  = None,
    top_fraction_of_variants_to_keep = 0.2,
    number_of_variants_per_gene_unit = 10,
    total_number_of_cycles_per_gene_unit = 10,
    number_of_tries_for_mutations_and_crossovers = 2,
    typical_gene_length = 10,
    end_cycles_if_no_improvement_for_n_cycles = 2,
    min_fraction_of_cycles_to_run = 0.1,
  )
##############################################################################
##########  EXAMPLE METHODS FOR GENETIC ALGORITHM WITH MULTIPROCESSING #######
##############################################################################

def example_new_gene_method(params, n):
  '''
  Default gene for genetic algorithm. Create n new genes
  You can use anything but it needs to have values for:
    length
    values
    score
    gene_id
   '''
  from scitbx.array_family import flex
  genes = []
  for i in range(n):
    gene_id = i
    length = params.typical_gene_length
    values = flex.random_double(length)
    score = None
    gene = group_args(
        group_args_type = 'default gene class as group_args',
        gene_id = gene_id,
        length = length,
        values = values,
        score = score)
    genes.append(gene)
  return genes

def example_mutation_method(params, gene):
  ''' Must take params and a gene and return a new gene '''
  from copy import deepcopy
  new_gene = deepcopy(gene)
  new_gene.values[np_random.randint(0,gene.length,1)[0]] += np_random.uniform(-4,4,1)[0]
  return new_gene

def example_recombination_method(params,gene_1, gene_2):
  ''' Must take params and two genes and return a new gene '''
  from copy import deepcopy
  new_gene = deepcopy(gene_1)
  i = np_random.randint(0,gene_1.length,1)[0]
  k = np_random.randint(0,2,1)[0]
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

def example_genes_are_identical_method(gene, other_gene):
  if gene.length != other_gene.length:
    return False
  for v1,v2 in zip(gene.values, other_gene.values):
    if v1 != v2:
      return False
  return True

##############################################################################
##################  GENETIC ALGORITHM WITH MULTIPROCESSING ###################
##############################################################################

class genetic_algorithm:
  def __init__(self,
      params,
      new_gene_method,
      mutation_method,
      recombination_method,
      scoring_method,
      genes_are_identical_method,
      genes = None,
      log = sys.stdout):

    # Save inputs
    adopt_init_args(self,locals())

    self.check_params()

    self.cycle = 0
    if not genes:
      self.genes = []
    self.run()
    self.update_random_seed()

  def check_params(self):
    '''  Make sure all expected values are there '''
    expected_keys = example_params().keys()
    found_keys = self.params().keys()
    for key in expected_keys:
      if not key in found_keys:
        print("Missing parameter: %s" %(key), file = self.log)
        assert key in found_keys # Missing parameter

  def update_random_seed(self):
    np_random.seed(self.params.random_seed)
    self.params.random_seed = np_random.randint(0,100000)

  def run(self):
    if not self.genes:
      self.get_genes()
    if len(self.genes)< 1:
      print("No genes available...quitting", file = self.log)
      return
    else:
      print("\nTotal working genes: %s" %(len(self.genes)),
         file = self.log)

    # How many cycles to try
    n_macro_cycles = self.get_number_of_macro_cycles()
    n_cycles = self.get_number_of_cycles()

    if n_macro_cycles:
      self.run_macro_cycles(n_macro_cycles, n_cycles)
    else:
      self.run_cycles(n_cycles)

    # Done, get best result
    best_gene = self.get_best_gene()
    print("Best gene with score of %s: \n%s" %(
      self.get_best_score_as_text(),self.get_best_gene()), file = self.log)


  def run_cycles(self, n_cycles):
      # Run n_cycles, but if number of genes is smaller than target
      #   then run more to generate more genes, up to 3x
      n_extra_requested = 0
      n_genes_target = self.get_number_of_variants_to_keep()
      last_improvement_info = group_args(
        group_args_type = 'last improvement info',
        last_improvement_cycle = None,
        last_improvement_score = None,
       )
      total_cycles_to_carry_out = n_cycles
      max_cycles = n_cycles * 3
      if self.params.min_fraction_of_cycles_to_run is None or \
          self.params.end_cycles_if_no_improvement_for_n_cycles is None:
        max_cycles_without_improvement = None
      else:
        max_cycles_without_improvement = max(
          n_cycles * self.params.min_fraction_of_cycles_to_run,
          2 *self.params.end_cycles_if_no_improvement_for_n_cycles)  # more here
      for cycle in range(max_cycles): # maximum
        self.one_cycle()
        if len(self.genes) < n_genes_target:
          total_cycles_to_carry_out += 1
        if cycle > min(max_cycles, total_cycles_to_carry_out): break

        new_best_score = self.get_best_score()
        if last_improvement_info.last_improvement_score is None or \
            new_best_score > last_improvement_info.last_improvement_score or \
            len(self.genes) < n_genes_target:  # doesn't count yet
          last_improvement_info.last_improvement_score = new_best_score
          last_improvement_info.last_improvement_cycle = cycle
        elif max_cycles_without_improvement is not None \
          and cycle - last_improvement_info.last_improvement_cycle > \
           max_cycles_without_improvement:
          print("\nEnding cycles as no improvement in past %s cycles " %(
          cycle - last_improvement_info.last_improvement_cycle),
              file = self.log)
          break

  def run_macro_cycles(self, n_macro_cycles, n_cycles):
      total_number_of_cycles = self.get_total_number_of_cycles()
      print("\nRunning about ",
       "%s macro_cycles of %s cycles each on %s processors" %(
        n_macro_cycles,n_cycles, self.params.nproc), file = self.log)
      max_macro_cycles = n_macro_cycles * 2
      working_macro_cycles = n_macro_cycles
      n_genes_target = self.get_number_of_variants_to_keep()

      print("Approximate total of %s cycles" %(
       total_number_of_cycles), file = self.log)
      last_improvement_info = group_args(
        group_args_type = 'last improvement info',
        last_improvement_cycle = None,
        last_improvement_score = None,
       )
      if self.params.min_fraction_of_cycles_to_run is None or \
          self.params.end_cycles_if_no_improvement_for_n_cycles is None:
        max_cycles_without_improvement = None
      else:
        max_cycles_without_improvement = max(
          n_macro_cycles * self.params.min_fraction_of_cycles_to_run,
          self.params.end_cycles_if_no_improvement_for_n_cycles)
      for macro_cycle in range(max_macro_cycles):
        if len(self.genes) < n_genes_target:
          working_macro_cycles += 1 # need another cycle

        if macro_cycle > min(working_macro_cycles, max_macro_cycles):
          break # done
        print( "\nMacro-cycle ",
         "%s (Genes: %s current top score: %s length: %s nproc: %s)" %(
          macro_cycle, len(self.genes), self.get_best_score_as_text(),
          self.get_best_gene().length if self.get_best_gene() else None,
           self.params.nproc),
          file = self.log)
        local_params = group_args(**self.params().copy())
        local_params.number_of_macro_cycles = 0
        local_params.total_number_of_cycles = n_cycles
        self.run_something(params = local_params,
           macro_cycle = True)
        self.score_genes()
        # Select the best genes
        self.select_top_variants()

        new_best_score = self.get_best_score()
        if last_improvement_info.last_improvement_score is None or \
            new_best_score > last_improvement_info.last_improvement_score:
          last_improvement_info.last_improvement_score = new_best_score
          last_improvement_info.last_improvement_cycle = macro_cycle
        elif max_cycles_without_improvement is not None \
           and macro_cycle - last_improvement_info.last_improvement_cycle > \
           max_cycles_without_improvement:
          print("\nEnding macro-cycles as no improvement in past %s cycles " %(
          macro_cycle - last_improvement_info.last_improvement_cycle),
              file = self.log)
          break # done with cycles


  def get_genes(self):
      # Create genes using supplied new_gene_method
      n = self.get_number_of_variants_to_make()
      self.update_random_seed()
      new_genes = self.new_gene_method(self.params, n)
      if new_genes:
        self.genes = new_genes
      self.score_genes()
      self.select_top_variants()


  def one_cycle(self):
    '''
     Run one cycle of optimization
     '''
    # Mutate genes to create new genes
    self.run_something(mutate = True)
    self.genes = make_unique(self.genes, self.genes_are_identical_method)

    if len(self.genes) > 1:
      # Recombine genes to create new genes
      self.run_something(recombine = True)
      self.genes = make_unique(self.genes, self.genes_are_identical_method)

    # Select the best genes
    self.select_top_variants()

  def run_something(self,
     params = None,
     genes = None,
     macro_cycle = None,
     create = None,
     mutate = None,
     recombine = None,
     score_only = None,
      ):

    if not params:
      params = self.params
    if not genes:
      genes = self.genes

    all_new_genes = []

    nproc = params.nproc
    end_number = -1
    if create:
      n_tot = self.get_number_of_variants_to_make()
    else:
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
        random_seed = np_random.randint(0,100000),
        start_number = start_number,
        end_number = end_number,
        ))

    local_params = group_args(**params().copy())
    local_params.nproc = 1 # Required

    kw_dict = {
      'params':local_params,
      'genes':genes,
      'create':create,
      'mutate':mutate,
      'macro_cycle':macro_cycle,
      'recombine':recombine,
      'score_only':score_only,
      'new_gene_method':self.new_gene_method,
      'mutation_method':self.mutation_method,
      'recombination_method':self.recombination_method,
      'genes_are_identical_method':self.genes_are_identical_method,
      'scoring_method':self.scoring_method,
     }

    runs_carried_out = run_jobs_with_large_fixed_objects(
      nproc = nproc,
      verbose = False,
      kw_dict = kw_dict,
      run_info_list = runs_to_carry_out,
      job_to_run = group_of_run_something,
      log = self.log)
    for run_info in runs_carried_out:
      new_genes = run_info.result.new_genes
      if new_genes:
        all_new_genes += new_genes
    all_new_genes = make_unique(all_new_genes, self.genes_are_identical_method)

    if score_only or create or macro_cycle:  # keep id and replace genes
      self.genes = all_new_genes

    else:  # usual
      self.set_gene_id_values(all_new_genes)

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


  def set_gene_id_values(self, new_genes):
    if not new_genes:
      return # nothing to do
    gene_id = self.get_highest_gene_id() + 1
    for gene in new_genes:
      gene.gene_id = gene_id
      gene_id += 1

  def get_highest_gene_id(self):
    highest = None
    for gene in self.genes:
      if highest is None or gene.gene_id > highest:
        highest = gene.gene_id
    if highest is None:
      highest = 0
    return highest
  def get_best_gene(self):
    if not self.genes:
      return None
    self.sort_genes()
    return self.genes[0]

  def _sort_genes_key(self, gene):
    if not gene.score:
      return -sys.maxsize
    else:
      return gene.score

  def sort_genes(self):
    self.genes = sorted(self.genes,
      key = self._sort_genes_key, reverse = True)

  def select_top_variants(self):
    if not self.genes: return
    self.sort_genes()

    n = self.get_number_of_variants_to_keep()
    self.genes = self.genes[:n]

  def score_genes(self):
    ''' Score all genes. If they already have a score...keep it'''
    for gene in self.genes:
      if gene.score is None:
         gene.score = self.scoring_method(self.params, gene)
    self.sort_genes()

  def get_number_of_variants_to_keep(self, minimum_variants = 2):
    if self.params.number_of_variants is not None:
      return max(minimum_variants,int(0.5+self.params.number_of_variants * \
       self.params.top_fraction_of_variants_to_keep))
    else:
      return max(minimum_variants,int(0.5+self.params.typical_gene_length * \
         self.params.number_of_variants_per_gene_unit *\
       self.params.top_fraction_of_variants_to_keep))

  def get_number_of_variants_to_make(self):
    if self.params.number_of_variants is not None:
      return self.params.number_of_variants
    else:
      return max(1,int(0.5+self.params.typical_gene_length * \
         self.params.number_of_variants_per_gene_unit))

  def get_total_number_of_cycles(self):
    #  Run a total of self.params.total_number_of_cycles, but if we have nproc>1
    #    run them in nproc groups of total_number_of_cycles/nproc
    #  total cycles = nproc * n_macro_cycles * n_cycles

    if self.params.total_number_of_cycles is not None:
      return self.params.total_number_of_cycles
    else:
      return max(1,int(0.5+self.params.typical_gene_length * \
         self.params.total_number_of_cycles_per_gene_unit))

  def get_number_of_cycles(self):
    if self.params.number_of_cycles is not None:
      return self.params.number_of_cycles
    else:
      #  Run a total of self.params.total_number_of_cycles,
      #    in nproc groups of ncycles
      nproc = self.params.nproc if self.params.nproc is not None else 1
      n_total_cycles = self.get_total_number_of_cycles()
      n_macro_cycles = self.get_number_of_macro_cycles()
      n_cycles = max(1,
        int(0.999 + n_total_cycles//max(1,n_macro_cycles*nproc)))
      return n_cycles

  def get_number_of_macro_cycles(self):
    if self.params.number_of_macro_cycles is not None:
      return self.params.number_of_macro_cycles
    else:
      nproc = self.params.nproc if self.params.nproc is not None else 1
      if nproc == 1:
         return 1
      else:
        n_total = self.get_total_number_of_cycles()
        n_macro_cycles = int(0.9999+(n_total/nproc)**0.5)
        return n_macro_cycles

def group_of_run_something(
        run_info,
        params,
        genes = None,
        macro_cycle= None,
        recombine = None,
        create = None,
        mutate = None,
        score_only = None,
        new_gene_method = None,
        mutation_method = None,
        recombination_method = None,
        scoring_method = None,
        genes_are_identical_method = None,
        log = sys.stdout):

  np_random.seed(run_info.random_seed)  # different for each run
  params = group_args(**params().copy())
  params.random_seed = np_random.randint(0,100000)

  if macro_cycle:  # Run a macro-cycle
    ga = genetic_algorithm(
      genes = genes,
      params = params,
      new_gene_method = new_gene_method,
      mutation_method = mutation_method,
      recombination_method = recombination_method,
      scoring_method = scoring_method,
      genes_are_identical_method = genes_are_identical_method,
      log = null_out(),
     )
    return group_args(
      group_args_type = 'set of genes after running one macro_cycle',
      new_genes = ga.genes,
     )

  # Usual
  new_genes = []
  for index in range(run_info.start_number, run_info.end_number + 1):
    params.random_seed = np_random.randint(0,100000)
    info = run_one_something(
      params,
      genes,
      index,
      recombine,
      create,
      mutate,
      score_only,
      new_gene_method,
      mutation_method,
      recombination_method,
      scoring_method,
      genes_are_identical_method,
      log = log)
    if info and info.new_genes:
      new_genes += info.new_genes

  # Make sure all new ones are unique
  if (not score_only) and genes:  # we have existing ones and not just scoring
    new_genes = make_unique(new_genes, genes_are_identical_method)

  return group_args(
      group_args_type = ' one set of recombined/mutated genes',
      new_genes = new_genes,
    )

def make_unique(new_genes, genes_are_identical_method):
  unique_new_genes = []
  for new_gene in new_genes:
    is_dup = False
    for gene in unique_new_genes:
      if genes_are_identical_method(gene, new_gene):
         is_dup = True
         break
    if not is_dup:
      unique_new_genes.append(new_gene)
  return unique_new_genes

def remove_duplicates_of_existing(existing_genes = None,
    new_genes = None,
    genes_are_identical_method = None):
  unique_new_genes = []
  for new_gene in new_genes:
    is_dup = False
    for gene in existing_genes:
      if genes_are_identical_method(gene, new_gene):
        is_dup = True
        break
    if not is_dup:
      unique_new_genes.append(new_gene)
  return unique_new_genes

def run_one_something(
        params,
        genes = None,
        index = None,
        recombine = None,
        create = None,
        mutate = None,
        score_only = None,
        new_gene_method = None,
        mutation_method = None,
        recombination_method = None,
        scoring_method = None,
        genes_are_identical_method = None,
        log = sys.stdout):

    assert (create, mutate,  recombine, score_only).count(True) == 1
    new_genes = []
    if not genes and not create:
      return new_genes

    if create: # create a gene
      new_gene = new_gene_method(params)
      if new_gene:
        new_genes.append(new_gene)

    elif score_only:
      new_genes.append(genes[index])

    elif mutate:
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
      other_gene_index = np_random.randint(0,len(genes) - 1,1)[0]
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
      group_args_type = ' one set of recombined/mutated/scored genes',
      new_genes = new_genes,
    )


def create_new_gene_by_mutation(params,
    gene,
    mutation_method,
    log = sys.stdout):
  n = np_random.poisson(params.mutation_rate,1)[0]
  if n == 0:
    return # nothing to do
  n_mutations_made = 0
  original_gene = gene
  for i in range(params.number_of_tries_for_mutations_and_crossovers * n):
    new_gene = mutation_method(params, gene)
    if new_gene:
      gene = new_gene
      n_mutations_made += 1
    if n_mutations_made >= n:
      break
  if n_mutations_made >= n:
    return gene
  else:
    return None

def create_new_gene_by_recombination(params,
   gene,
   other_gene,
   recombination_method,
   log = sys.stdout):

  if np_random.uniform(0,1,1)[0] > params.recombination_rate:
    return # nothing to do

  return recombination_method(params, gene, other_gene)

##############################################################################
##################  END GENETIC ALGORITHM WITH MULTIPROCESSING ###############
##############################################################################
