from __future__ import division, print_function

class pdb_handler(object):
  '''
  Author      : Uervirojnangkoorn, M.
  Created     : 7/13/2014
  Determines number of atoms in a given pdb or sequence file.
  '''

  def __init__(self, file_name_pdb):
    '''
    Constructor
    '''
    file_pdb = open(file_name_pdb,'r')
    self.data_pdb = file_pdb.read().split("\n")


  def print_file_content(self):

    for line in self.data_pdb:
      print(line)

  def atom_stats(self):

    cn_C = 0
    cn_H = 0
    cn_N = 0
    cn_O = 0
    cn_S = 0


    for line in self.data_pdb:
      if line.find('ATOM') == 0 or line.find('HETATM') == 0:
        col = line.split(' ')
        atom = col[len(col)-3].strip()
        if atom == 'N':
          cn_N += 1
        elif atom == 'O':
          cn_O += 1
        elif atom == 'C':
          cn_C += 1
        elif atom == 'S':
          cn_S += 1
        elif atom == 'H':
          cn_H += 1


    #print 'C ', cn_C, 'H', cn_H, 'N', cn_N, 'O', cn_O, 'S', cn_S

    return (cn_C, cn_H, cn_N, cn_O, cn_S)
