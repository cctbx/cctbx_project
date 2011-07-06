## Peter Zwart Mai 10, 2005
from cctbx.array_family import flex
from scitbx.math import chebyshev_polynome
from iotbx import data_plots
import scitbx.math
import math
import sys



class p_vm_calculator(object):
  """ solvent content/matthews probability calculator"""
  def __init__(self, miller_array, n_residues, n_bases=None, out=None, verbose=0):
    self.unit_cell_volume = miller_array.unit_cell().volume()
    self.z = miller_array.space_group().order_z()

    if  n_residues is not None:
      self.n_residues = n_residues
    else:
      self.n_residues = 0.0
    if n_bases is not None:
      self.n_bases = n_bases
    else:
      self.n_bases=0.0
    assert (self.n_bases+self.n_residues>0)

    self.mw_residue = 112.5
    self.mw_base = 311.0

    self.rho_spec = 0.74 # only for protein, not DNA/needs to be modified

    self.vm_prop = []

    if out is None:
      out=sys.stdout
    self.out = out
    self.verbose = verbose

    coefs = flex.double([-14.105436736742137,    -0.47015366358636385,
                          -2.9151681976244639,   -0.49308859741473005,
                           0.90132625209729045,   0.033529051311488103,
                           0.088901407582105796,  0.10749856607909694,
                           0.055000918494099861, -0.052424473641668454,
                          -0.045698882840119227,  0.076048484096718036,
                          -0.097645159906868589,  0.03904454313991608,
                          -0.072186667173865071])

    self.log_p_solc = chebyshev_polynome(15,0,1,coefs)
    self. vm_prop_table()
    self.best_guess = self.guesstimate()

  def vm(self,copies):
    result = None
    den = ((self.n_residues*self.mw_residue+self.n_bases*self.mw_base)*\
      self.z*copies)
    if(den != 0.0):
      result = self.unit_cell_volume/den
    return(result)

  def solc(self,vm):
    return (1.0-self.rho_spec/(0.602*vm))

  def p_solc_calc(self, sc):
    if (sc>1):
      sc=1
    if (sc<0):
      sc=0
    tmp = self.log_p_solc.f(sc)
    return(math.exp(tmp))

  def vm_prop_table(self):
    solc = 1.0
    n_copies = 0.0
    tot_p = 0.0
    while solc > 0:
      n_copies+=1.0
      vm = self.vm(n_copies)
      solc = self.solc(vm)
      p_vm = self.p_solc_calc(solc)
      entry = [ n_copies, solc, vm, p_vm ]
      tot_p += p_vm
      self.vm_prop.append(entry)

    tmp = self.vm_prop.pop() ## The last item has a negative solvent content
    tot_p -= tmp[3]
    if (int(n_copies)==1):
      print >> self.out, "Too many residues to fit in the ASU"
      print >> self.out, "  resetting numer of residues in monomer to %5.0f" \
            %(self.n_residues/10.0)
      self.n_residues/=10.0
      self.n_bases/=10.0
      self.vm_prop_table()

    for ii in range(int(n_copies)-1):
      self.vm_prop[ii][3] = self.vm_prop[ii][3]/tot_p

  def guesstimate(self):
    max = 0.0
    guess = 1;
    for ii in range(len(self.vm_prop)):
      if (self.vm_prop[ii][3]>max):
        max = self.vm_prop[ii][3]
        guess = ii+1
    return guess

def matthews_rupp(miller_array,
                  n_residues=None,
                  n_bases=None,
                  out=None,verbose=0):
  """Probabilistic estimation of number of copies in the asu"""
  if out is None:
    out = sys.stdout

  if (n_residues==None):
    if (n_bases==None):
      print >> out
      print >> out, "Number of residues unknown, assuming 50% solvent content"
      n_residues=1
      n_bases=0
      verbose=0
      print >> out

  vm_estimator = p_vm_calculator(miller_array, n_residues, n_bases)

  if verbose>0:
    print >> out,"----------------------------------------------------------------"
    print >> out,"| Copies | Solvent content | Matthews Coef. | P(solvent cont.) |"
    print >> out,"|--------|-----------------|----------------|------------------|"
    for ii in range( len(vm_estimator.vm_prop) ):
      print >> out,"|%7.0f" %(vm_estimator.vm_prop[ii][0]),"|" \
            "%11.3f" %(vm_estimator.vm_prop[ii][1]),"     |" \
            "%11.3f" %(vm_estimator.vm_prop[ii][2]),"    |" \
            "%12.3f" %(vm_estimator.vm_prop[ii][3]),"     |"
  print >> out,"----------------------------------------------------------------"
  if verbose>0:
    print >> out,"|              Best guess : %4.0f" %(vm_estimator.best_guess),\
          " copies in the asu            |"
  if verbose<=0:
    print >> out,"|              Best guess : %4.0f" %(vm_estimator.best_guess),\
          " residues in the asu          |"
  print >> out, "----------------------------------------------------------------"
  if verbose==0:
    vm_estimator.n_residues = vm_estimator.best_guess
    vm_estimator.n_bases=0.0
    vm_estimator.best_guess = 1

  matthews_table = data_plots.table_data(
    title="Solvent content analysis",
    column_labels=["Copies", "Solvent content", "Matthews coeff.",
                   "P(solvent content)"],
    graph_names=["Solvent content", "Matthews coeff."],
    graph_columns=[[0,1], [0,2]])
  for ii in range( len(vm_estimator.vm_prop) ):
    matthews_table.add_row(vm_estimator.vm_prop[ii])

  return( [vm_estimator.n_residues,
           vm_estimator.n_bases,
           vm_estimator.best_guess,
           vm_estimator.vm_prop[vm_estimator.best_guess-1][1],
           matthews_table] )
