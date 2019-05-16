from __future__ import absolute_import, division, print_function
from six.moves import range
'''
Author      : Uervirojnangkoorn, M.
Created     : 12/1/2014
Description : Genetic algorithm.
'''
import random

class ga_handler(object):
  '''
  Genetic algorithm class.
  '''
  def __init__(self):
    '''
    Constructor
    '''


  def initialse(self, pop_size, idv_length, cdf_from_pdf_set, phi_for_hl):
    '''
    generate population from cdf of HL
    '''
    pop=[];
    for i_pop in range(pop_size):
      idv=[0]*idv_length;

      for i_idv in range(idv_length):
        tmp_rand=random.random();
        for i_phi in range(len(phi_for_hl)):
          if cdf_from_pdf_set[i_idv][i_phi]>= tmp_rand:
            idv[i_idv]=phi_for_hl[i_phi]
            break

      pop.append(idv)

    return pop;


  def crossover(self,parent1,parent2,ratio_cross):

    num_point_cross=int(round(ratio_cross*len(parent1)))
    child1=parent1[:]
    child2=parent2[:]

    ''' unicross '''
    cross_template=random.sample(range(len(parent1)),num_point_cross)

    for i_cross in range(num_point_cross):
      child1[cross_template[i_cross]]=parent2[cross_template[i_cross]]
      child2[cross_template[i_cross]]=parent1[cross_template[i_cross]]


    ''' normal cross
    max_i=int(round(len(parent1)/(num_point_cross*2)))

    for i in range(max_i):
      i_st=num_point_cross*i*2
      i_en=i_st+num_point_cross
      child1[i_st:i_en]=parent2[i_st:i_en]
      child2[i_st:i_en]=parent1[i_st:i_en]
    '''


    return child1,child2,cross_template

  def mutation(self,parent,prob_of_mut,num_point_mut,cdf_from_pdf_set,phi_for_hl):

    child=parent[:];
    if random.random() < prob_of_mut:
      mut_template=random.sample(range(len(parent)),num_point_mut);

      for i_mut in range(num_point_mut):
        tmp_rand=random.random();
        dif_abs=[0]*len(phi_for_hl);
        for i_phi in range(len(phi_for_hl)):
          dif_abs[i_phi]=abs(cdf_from_pdf_set[mut_template[i_mut]][i_phi]-tmp_rand);


        child[mut_template[i_mut]]=phi_for_hl[dif_abs.index(min(dif_abs))];


    return child;
