# -*- coding: utf-8 -*-

from __future__ import division, print_function

def tst_01():
  '''  Exercise genetic algorithm '''

  from scitbx.genetic_algorithm_with_mp import  genetic_algorithm,\
    example_params, example_new_gene_method, example_mutation_method, \
    example_recombination_method, example_scoring_method,\
    example_genes_are_identical_method

  params = example_params
  new_gene_method = example_new_gene_method
  mutation_method = example_mutation_method
  recombination_method = example_recombination_method
  scoring_method = example_scoring_method
  genes_are_identical_method = example_genes_are_identical_method

  result = genetic_algorithm(
    params = params,
    new_gene_method = new_gene_method,
    mutation_method = mutation_method,
    recombination_method = recombination_method,
    scoring_method = scoring_method,
    genes_are_identical_method = genes_are_identical_method,
   )
  best_gene = result.get_best_gene()
  print (best_gene.score)
  print (list(best_gene.values))
  assert abs(best_gene.score  - 10 ) < 0.1

if __name__ == "__main__":
  tst_01()
