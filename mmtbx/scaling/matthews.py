## Peter Zwart Mai 10, 2005
## refactored by Gabor & Nat 20111104
import math
import sys

log_p_solc = None
def get_log_p_solc () :
  global log_p_solc
  from cctbx.array_family import flex
  from scitbx.math import chebyshev_polynome
  if (log_p_solc is None) :
    coeffs = flex.double([
      -14.105436736742137,    -0.47015366358636385,
      -2.9151681976244639,   -0.49308859741473005,
      0.90132625209729045,   0.033529051311488103,
      0.088901407582105796,  0.10749856607909694,
      0.055000918494099861, -0.052424473641668454,
      -0.045698882840119227,  0.076048484096718036,
      -0.097645159906868589,  0.03904454313991608,
      -0.072186667173865071])
    log_p_solc = chebyshev_polynome( 15, 0, 1, coeffs)
  return log_p_solc

def p_solc_calc(sc):
  """
  Calculate solvent fraction probability
  """
  if sc < 0 or 1.0 < sc:
    raise ValueError, "solvent content out of range"
  log_p_solc = get_log_p_solc()
  return math.exp(log_p_solc.f(sc))

class density_calculator(object):
  """
  Calculate Matthews coefficient and solvent fraction
  """
  def __init__(self, crystal):
    self.asu_volume = (
      crystal.unit_cell().volume() / crystal.space_group().order_z())

  def vm(self, weight):
    if weight <= 0:
      raise ValueError, "Incorrect weight"
    return self.asu_volume / weight

  def macromolecule_fraction(self, weight, rho_spec):
    return rho_spec / ( 0.602 * self.vm( weight = weight ) )

  def solvent_fraction(self, weight, rho_spec):
    return 1.0 - self.macromolecule_fraction(
      weight = weight,
      rho_spec = rho_spec)

class component(object):
  """
  Macromolecule component
  """
  def __init__(self, mw, rho_spec):
    self.mw = mw
    self.rho_spec = rho_spec

  def protein(cls, nres):
    return cls( mw = nres * 112.5, rho_spec = 0.74 )

  def nucleic(cls, nres):
    # Kantardjieff & Rupp, Prot. Sci. 12, 1865-1871
    return cls( mw = nres * 311.0, rho_spec = 0.50 )

  protein = classmethod( protein )
  nucleic = classmethod( nucleic )

def number_table(components, density_calculator):
  if not components:
    raise ValueError, "Empty component list"
  unit_mfrac = sum([
    density_calculator.macromolecule_fraction(
      weight = c.mw,
      rho_spec = c.rho_spec) for c in components ])
  results = []
  count = 1
  while True:
    sc = 1.0 - count * unit_mfrac
    assert sc <= 1.0
    if sc < 0:
      break
    results.append( ( count, p_solc_calc( sc = sc ) ) )
    count += 1
  return results

class p_vm_calculator(object):
  """solvent content/matthews probability calculator"""
  def __init__(self, crystal_symmetry, n_residues, n_bases=None, out=None,
      verbose=0):
    assert ((crystal_symmetry.unit_cell() is not None) and
            (crystal_symmetry.space_group() is not None))
    self.unit_cell_volume = crystal_symmetry.unit_cell().volume()
    self.z = crystal_symmetry.space_group().order_z()
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
    return p_solc_calc( sc = max( min( sc, 1.0 ), 0 ) )

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

def matthews_rupp(crystal_symmetry,
                  n_residues=None,
                  n_bases=None,
                  out=None,
                  verbose=0):
  """Probabilistic estimation of number of copies in the asu"""
  from iotbx import data_plots
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
  vm_estimator = p_vm_calculator(crystal_symmetry, n_residues, n_bases)
  if verbose>0:
    print >> out,"-" * 64
    print >> out,"| Copies | Solvent content | Matthews Coef. | P(solvent cont.) |"
    print >> out,"|--------|-----------------|----------------|------------------|"
    for ii in range( len(vm_estimator.vm_prop) ):
      print >> out,"|%7.0f" %(vm_estimator.vm_prop[ii][0]),"|" \
            "%11.3f" %(vm_estimator.vm_prop[ii][1]),"     |" \
            "%11.3f" %(vm_estimator.vm_prop[ii][2]),"    |" \
            "%12.3f" %(vm_estimator.vm_prop[ii][3]),"     |"
  print >> out,"-" * 64
  if verbose>0:
    print >> out,"|              Best guess : %4.0f" %(vm_estimator.best_guess),\
          " copies in the asu            |"
  if verbose<=0:
    print >> out,"|              Best guess : %4.0f" %(vm_estimator.best_guess),\
          " residues in the asu          |"
  print >> out, "-" * 64
  if verbose==0:
    vm_estimator.n_residues = vm_estimator.best_guess
    vm_estimator.n_bases=0.0
    vm_estimator.best_guess = 1 # XXX why?
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
           matthews_table])

########################################################################
# REGRESSION
def exercise () :
  from cctbx import crystal
  from libtbx.test_utils import approx_equal, contains_substring
  from libtbx.utils import null_out
  from cStringIO import StringIO
  symm = crystal.symmetry(
    unit_cell=(77.3,107.6,84.4,90,94.2,90),
    space_group_symbol="P21")
  out = StringIO()
  results = matthews_rupp(
    crystal_symmetry=symm,
    n_residues=153,
    out=out,
    verbose=1)
  assert contains_substring(out.getvalue(),
    """|      1 |      0.940      |     20.337     |       0.000      |""")
  assert (results[2] == 8)
  assert approx_equal(results[3], 0.5165, eps=0.0001)
  results = matthews_rupp(
    crystal_symmetry=symm,
    out=null_out(),
    verbose=1)
  assert (results[0] == 1281)

if (__name__ == "__main__") :
  exercise()
  print "OK"
