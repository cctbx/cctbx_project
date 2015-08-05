from __future__ import division
import cctbx.array_family.flex # import dependency
import boost.python
ext = boost.python.import_ext("mmtbx_ncs_ext")
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import sgtbx

pdb_str = """
CRYST1   21.954   18.566   12.975  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       7.622  14.852   6.274  1.00 20.00           N
ATOM      2  CA  ALA A   1       7.323  14.430   7.635  1.00 20.00           C
ATOM      3  C   ALA A   1       6.442  13.183   7.656  1.00 20.00           C
ATOM      4  O   ALA A   1       6.848  12.171   8.227  1.00 20.00           O
ATOM      5  CB  ALA A   1       6.675  15.566   8.414  1.00 20.00           C
ATOM      6  N   ALA A   2       5.256  13.245   7.012  1.00 20.00           N
ATOM      7  CA  ALA A   2       4.288  12.143   6.916  1.00 20.00           C
ATOM      8  C   ALA A   2       4.847  10.948   6.125  1.00 20.00           C
ATOM      9  O   ALA A   2       4.553   9.797   6.462  1.00 20.00           O
ATOM     10  CB  ALA A   2       3.000  12.642   6.276  1.00 20.00           C
ATOM     11  N   ALA A   3       5.664  11.231   5.084  1.00 20.00           N
ATOM     12  CA  ALA A   3       6.313  10.229   4.231  1.00 20.00           C
ATOM     13  C   ALA A   3       7.391   9.457   4.994  1.00 20.00           C
ATOM     14  O   ALA A   3       7.616   8.281   4.687  1.00 20.00           O
ATOM     15  CB  ALA A   3       6.915  10.892   3.000  1.00 20.00           C
ATOM     16  N   ALA A   4       8.063  10.119   5.976  1.00 20.00           N
ATOM     17  CA  ALA A   4       9.089   9.491   6.822  1.00 20.00           C
ATOM     18  C   ALA A   4       8.441   8.418   7.703  1.00 20.00           C
ATOM     19  O   ALA A   4       8.888   7.272   7.671  1.00 20.00           O
ATOM     20  CB  ALA A   4       9.785  10.535   7.683  1.00 20.00           C
ATOM     21  N   ALA A   5       7.338   8.773   8.423  1.00 20.00           N
ATOM     22  CA  ALA A   5       6.573   7.863   9.298  1.00 20.00           C
ATOM     23  C   ALA A   5       6.042   6.663   8.513  1.00 20.00           C
ATOM     24  O   ALA A   5       5.981   5.554   9.048  1.00 20.00           O
ATOM     25  CB  ALA A   5       5.433   8.608   9.975  1.00 20.00           C
ATOM     26  N   ALA A   6       5.715   6.893   7.221  1.00 20.00           N
ATOM     27  CA  ALA A   6       5.258   5.889   6.260  1.00 20.00           C
ATOM     28  C   ALA A   6       6.440   4.995   5.865  1.00 20.00           C
ATOM     29  O   ALA A   6       6.361   3.780   6.060  1.00 20.00           O
ATOM     30  CB  ALA A   6       4.664   6.568   5.032  1.00 20.00           C
ATOM     31  N   ALA A   7       7.556   5.609   5.374  1.00 20.00           N
ATOM     32  CA  ALA A   7       8.788   4.918   4.962  1.00 20.00           C
ATOM     33  C   ALA A   7       9.411   4.069   6.075  1.00 20.00           C
ATOM     34  O   ALA A   7       9.954   3.000   5.782  1.00 20.00           O
ATOM     35  CB  ALA A   7       9.806   5.915   4.433  1.00 20.00           C
TER
ATOM     36  N   ALA B   1      16.622  14.862   6.586  1.00 20.00           N
ATOM     37  CA  ALA B   1      16.323  14.369   7.923  1.00 20.00           C
ATOM     38  C   ALA B   1      15.442  13.123   7.879  1.00 20.00           C
ATOM     39  O   ALA B   1      15.848  12.082   8.396  1.00 20.00           O
ATOM     40  CB  ALA B   1      15.675  15.463   8.761  1.00 20.00           C
ATOM     41  N   ALA B   2      14.256  13.218   7.239  1.00 20.00           N
ATOM     42  CA  ALA B   2      13.288  12.123   7.086  1.00 20.00           C
ATOM     43  C   ALA B   2      13.847  10.971   6.233  1.00 20.00           C
ATOM     44  O   ALA B   2      13.553   9.804   6.509  1.00 20.00           O
ATOM     45  CB  ALA B   2      12.000  12.655   6.473  1.00 20.00           C
ATOM     46  N   ALA B   3      14.664  11.308   5.208  1.00 20.00           N
ATOM     47  CA  ALA B   3      15.313  10.352   4.304  1.00 20.00           C
ATOM     48  C   ALA B   3      16.391   9.541   5.026  1.00 20.00           C
ATOM     49  O   ALA B   3      16.616   8.383   4.658  1.00 20.00           O
ATOM     50  CB  ALA B   3      15.915  11.078   3.110  1.00 20.00           C
ATOM     51  N   ALA B   4      17.063  10.151   6.041  1.00 20.00           N
ATOM     52  CA  ALA B   4      18.089   9.479   6.853  1.00 20.00           C
ATOM     53  C   ALA B   4      17.441   8.362   7.677  1.00 20.00           C
ATOM     54  O   ALA B   4      17.888   7.219   7.585  1.00 20.00           O
ATOM     55  CB  ALA B   4      18.785  10.477   7.767  1.00 20.00           C
ATOM     56  N   ALA B   5      16.338   8.679   8.414  1.00 20.00           N
ATOM     57  CA  ALA B   5      15.573   7.724   9.240  1.00 20.00           C
ATOM     58  C   ALA B   5      15.042   6.567   8.394  1.00 20.00           C
ATOM     59  O   ALA B   5      14.981   5.431   8.870  1.00 20.00           O
ATOM     60  CB  ALA B   5      14.433   8.433   9.955  1.00 20.00           C
ATOM     61  N   ALA B   6      14.715   6.864   7.115  1.00 20.00           N
ATOM     62  CA  ALA B   6      14.258   5.912   6.103  1.00 20.00           C
ATOM     63  C   ALA B   6      15.440   5.040   5.662  1.00 20.00           C
ATOM     64  O   ALA B   6      15.361   3.816   5.793  1.00 20.00           O
ATOM     65  CB  ALA B   6      13.664   6.654   4.912  1.00 20.00           C
ATOM     66  N   ALA B   7      16.556   5.678   5.204  1.00 20.00           N
ATOM     67  CA  ALA B   7      17.788   5.010   4.756  1.00 20.00           C
ATOM     68  C   ALA B   7      18.411   4.104   5.823  1.00 20.00           C
ATOM     69  O   ALA B   7      18.954   3.052   5.475  1.00 20.00           O
ATOM     70  CB  ALA B   7      18.806   6.033   4.280  1.00 20.00           C
TER
END
"""

def run():
  #
  # Read PDB file from string above, create xray_structure object
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  pdb_inp.write_pdb_file(file_name="model.pdb")
  xray_structure = pdb_inp.xray_structure_simple()
  #
  # Calculate "Fobs" from this model
  #
  f_obs = abs(xray_structure.structure_factors(d_min=2).f_calc())
  f_obs.set_sigmas(sigmas = flex.double(f_obs.data().size(), 1.0))
  #
  # Create resolution bins
  #
  f_obs.setup_binner(reflections_per_bin = 250)
  binner = f_obs.binner()
  n_bins = binner.n_bins_used()
  for i_bin in binner.range_used():
    bin_sel = f_obs.binner().selection(i_bin)
    f_obs_bin = f_obs.select(bin_sel)
    print "bin: %d n_refl.: %d" % (i_bin, f_obs_bin.data().size()), \
      "%6.3f-%-6.3f"%f_obs_bin.d_max_min()
  #
  # Construct rbin array used in tncs_eps_factor_refinery.
  # Very inefficient, ok for now..
  #
  rbin = flex.int(f_obs.data().size(), -1)
  for i_bin in binner.range_used():
    for i_seq in binner.array_indices(i_bin):
      rbin[i_seq] = i_bin-1 # i_bin starts with 1, not 0 !
  assert flex.min(rbin)==0
  assert flex.max(rbin)==n_bins-1
  #
  # Create NCS pair object that contains all information we will need.
  # This is C++ container implemented in cctbx_project/mmtbx/ncs/tncs.h
  #
  ncs_pair = ext.pair(
    r = ([1,0,0,0,0.998630,0.052336,0,-0.052336,0.998630]), # I have code to get this
    t = ([-9,0,0]),            # I have code to get this
    radius=4.24,                     # XXX: need meaningful numbers here!
    weight=20,                     # XXX: need meaningful numbers here!
    fracscat=0.5,                  # XXX: need meaningful numbers here!
    rho_mn=flex.double(n_bins,1) ) # XXX: need meaningful numbers here!
  #
  # Prepare data that we need in order to call tncs_eps_factor_refinery
  #
  # List of symmetry rotation matrices
  sym_matrices = []
  for m_as_string in f_obs.space_group().smx():
    o = sgtbx.rt_mx(symbol=str(m_as_string), t_den=f_obs.space_group().t_den())
    m_as_double = o.r().as_double()
    print m_as_string, m_as_double
    sym_matrices.append(m_as_double)
  #
  obj = ext.tncs_eps_factor_refinery(
    tncs_pairs               = [ncs_pair],
    f_obs                    = f_obs.data(),
    sigma_f_obs              = f_obs.sigmas(),
    rbin                     = rbin,
    SigmaN                   = flex.double(f_obs.data().size(),1), # XXX: need meaningful numbers here!
    space_group              = f_obs.space_group(),
    miller_indices           = f_obs.indices(),
    fractionalization_matrix = f_obs.unit_cell().fractionalization_matrix(),
    sym_matrices             = sym_matrices)
  print "target:",obj.target()
  print "gradient:", list(obj.gradient())

if (__name__ == "__main__"):
  run()
