import math,exceptions,sys
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from rstbx.dps_core import dps_core
from rstbx.dps_core.sampling import hemisphere_shortcut
from rstbx.dps_core.basis_choice import SelectBasisMetaprocedure

# reciprocal space vectors expressed internally in units of inverse Angstroms:
# eV per KeV *
# Joules per eV /
# Joule-seconds (Planck's const) /
# meters/second (speed of light) /
# Angstroms per meter
keV_to_inv_Angstrom = 1E3             * \
                      1.602176487E-19 / \
                      6.62606896E-34  / \
                      299792458.E0    / \
                      1E10

def parse_input(filename):
  G = open(filename,"r")
  reciprocal_vectors = flex.vec3_double()
  lines = G.readlines()
  for line in lines[0:len(lines)-0]:
    tokens = line.strip().split('\t')
    assert len(tokens)==7
    q_vector = matrix.col([float(i) for i in tokens[4:7]])
    checklength = q_vector.dot(q_vector)
    assert approx_equal(checklength, 1.0, eps=0.001)
    energy_keV = float(tokens[2])
    theta = (math.pi/180.) * float(tokens[3])
    inv_lambda_Angstrom = keV_to_inv_Angstrom * energy_keV
    inv_d_Angstrom = 2. * math.sin(theta) * inv_lambda_Angstrom
    extended_q_Angstrom = inv_d_Angstrom * q_vector
    reciprocal_vectors.append(extended_q_Angstrom.elems)
  return reciprocal_vectors

def parse_synthetic(filename):
  G = open(filename,"r")
  reciprocal_vectors = flex.vec3_double()
  #Example of Silicon F d 3-bar m unit cell oriented along lab axes
  A_mat = matrix.sqr((5.43,0.,0.,0.,5.43,0.,0.,0.,5.43))
  A_star = A_mat.inverse()
  for line in G.readlines()[1:]:
    tokens = line.strip().split('\t')
    assert len(tokens)==5
    miller = matrix.col([float(i) for i in tokens[0:3]])
    reciprocal_vector = A_star*miller
    reciprocal_vectors.append(reciprocal_vector.elems)
  return reciprocal_vectors

class algorithm_parameters:
  max_cell_edge_for_dps_fft = 100.0 # in Angstroms.
                                    # 100 Angstrom should be suitable for small molecular work
                                    # use Bragg spot spacing to estimate the maximum
                                    # cell edge for macromolecular cases
  directional_sampling_granularity = 0.029 # in radians; 0.029 radian granularity
                                           # is both fast enough and granular enough
                                           # to sample the hemisphere for small molecular work.
                                           # Tradeoff between performance and coverage occurs
                                           # for macromolecular work; DPS paper uses 0.029
  max_cell_edge_basis_choice = 8.0  # in Angstroms. This input parameter is important.
                                    # choice of basis vector is extremely sensitive to
                                    # this parameter--for silicon, must choose a value less
                                    # than twice the cell edge
def do_index(reciprocal_space_vectors):
  D = dps_core()
  D.setMaxcell(algorithm_parameters.max_cell_edge_for_dps_fft)
  D.setXyzData(reciprocal_space_vectors)
  hemisphere_shortcut(ai = D,
    characteristic_sampling = algorithm_parameters.directional_sampling_granularity,
    max_cell = algorithm_parameters.max_cell_edge_basis_choice)
  M = SelectBasisMetaprocedure(D)
  print D.getOrientation().unit_cell(), D.getOrientation().unit_cell().volume()
    # de novo subgroup calc with iotbx code
  from rstbx.dps_core.lepage import iotbx_converter
  L = iotbx_converter(D.getOrientation().unit_cell().minimum_cell(),5.0)
  for subgroup in L:
      print subgroup.digest()
  supergroup = L[0]

  print "before", D.getOrientation().unit_cell(),D.getOrientation().unit_cell().volume()
  cb_op = supergroup['cb_op_inp_best'].c().as_double_array()[0:9]
  orient = D.getOrientation()
  print cb_op
  orient.change_basis(matrix.sqr(cb_op).transpose())
  constrain_orient = orient.constrain(supergroup['system'])
  D.setOrientation(constrain_orient)
  print "after", D.getOrientation().unit_cell(),D.getOrientation().unit_cell().volume()
  A_star = matrix.sqr(D.getOrientation().reciprocal_matrix())
  a_star = matrix.col((A_star[0],A_star[3],A_star[6]))
  print "+++++++++++++++++++++++++++++++"
  print "RMSDEV:",D.rmsdev()
  print "-------------------------------"

  for hkl,obs in zip(D.hklobserved(),D.observed()):
    displace = matrix.col(hkl) - matrix.col(obs)
    diff = math.sqrt(displace.dot(displace))
    print hkl,diff

def test_case_obs_data():
  R = parse_input("./si_brief.txt")
  do_index(R)

def test_case_synthetic_data():
  R = parse_synthetic("./si_synthetic.txt")
  do_index(R)

if __name__=='__main__':
  print "Test autoindexing on synthetic reciprocal space positions"
  test_case_synthetic_data()
  print "Test autoindexing on monoscan data"
  test_case_obs_data()
