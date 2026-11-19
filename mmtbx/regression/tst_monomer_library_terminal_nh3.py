"""
mmtbx.monomer_library.terminal_nh3 through pdb_interpretation, without
reduce2: an N-terminal NH3 of a peptide monomer whose dictionary does not
fit mod_NH3 (no atom H, or H1-H3 on another atom than N) is restrained with
the dictionary's N-H names plus computed names (generated modification
NH3-<comp>); mod_NH3 is not applied where it would take a carbon H's name;
SCY (H2 on N) and standard residues keep mod_NH3.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.cif
import iotbx.pdb
import libtbx.load_env
import mmtbx.monomer_library.server
from libtbx.utils import null_out
from mmtbx.monomer_library import pdb_interpretation, terminal_nh3

# 0A0-GLY: HN1, HN2 (dictionary) and HN3 (computed) on N
pdb_0a0 = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   0A0 A   1      24.101   8.113 -14.105  1.00 20.00           N
HETATM    2  CA  0A0 A   1      23.417   9.374 -13.816  1.00 20.00           C
HETATM    3  C   0A0 A   1      24.034  10.606 -14.545  1.00 20.00           C
HETATM    4  O   0A0 A   1      23.726  11.736 -14.108  1.00 20.00           O
HETATM    5  CB  0A0 A   1      23.377   9.577 -12.304  1.00 20.00           C
HETATM    6  CG  0A0 A   1      24.721   9.728 -11.566  1.00 20.00           C
HETATM    7  OD1 0A0 A   1      25.780   9.685 -12.231  1.00 20.00           O
HETATM    8  OD2 0A0 A   1      24.641   9.879 -10.329  1.00 20.00           O
HETATM    9  CM  0A0 A   1      21.979   9.248 -14.335  1.00 20.00           C
HETATM   10  HN1 0A0 A   1      24.064   7.949 -14.979  1.00 20.00           H
HETATM   11  HN2 0A0 A   1      23.704   7.452 -13.662  1.00 20.00           H
HETATM   12  HN3 0A0 A   1      24.952   8.172 -13.850  1.00 20.00           H
HETATM   13 HBC1 0A0 A   1      22.843  10.362 -12.106  1.00 20.00           H
HETATM   14 HBC2 0A0 A   1      22.901   8.835 -11.899  1.00 20.00           H
HETATM   15 HMC1 0A0 A   1      21.513   8.517 -13.900  1.00 20.00           H
HETATM   16 HMC2 0A0 A   1      21.475  10.060 -14.170  1.00 20.00           H
HETATM   17 HMC3 0A0 A   1      21.965   9.081 -15.290  1.00 20.00           H
ATOM     18  N   GLY A   2      23.223  11.325 -15.249  1.00  5.27           N
ATOM     19  CA  GLY A   2      23.645  12.533 -15.945  1.00  5.73           C
ATOM     20  C   GLY A   2      23.253  13.729 -15.098  1.00  5.60           C
ATOM     21  O   GLY A   2      22.111  13.831 -14.628  1.00  5.67           O
ATOM     22  H   GLY A   2      22.389  11.143 -15.355  1.00  5.27           H
ATOM     23  HA2 GLY A   2      24.606  12.529 -16.074  1.00  5.73           H
ATOM     24  HA3 GLY A   2      23.210  12.593 -16.810  1.00  5.73           H
END
"""

# 2KZ-GLY: H1, H3 on N (names 2KZ lacks), H2 on CB (the dictionary's)
pdb_2kz = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   2KZ A   1      21.584  21.923  19.934  1.00 20.00           N
HETATM    2  CA  2KZ A   1      20.688  20.877  19.924  1.00 20.00           C
HETATM    3  C   2KZ A   1      21.181  19.479  19.862  1.00 20.00           C
HETATM    4  O   2KZ A   1      20.664  18.687  19.102  1.00 20.00           O
HETATM    5  CB  2KZ A   1      19.356  21.125  19.970  1.00 20.00           C
HETATM    6  CG  2KZ A   1      18.411  20.018  20.124  1.00 20.00           C
HETATM    7  OD1 2KZ A   1      17.085  20.247  20.036  1.00 20.00           O
HETATM    8  OD2 2KZ A   1      18.822  18.892  20.330  1.00 20.00           O
HETATM    9  H1  2KZ A   1      21.225  22.604  20.317  1.00 20.00           H
HETATM   10  H2  2KZ A   1      19.023  21.991  19.905  1.00 20.00           H
HETATM   12  H3  2KZ A   1      22.291  21.691  20.366  1.00 20.00           H
ATOM     13  N   GLY A   2      20.825  18.621  20.776  1.00  5.27           N
ATOM     14  CA  GLY A   2      21.259  17.231  20.813  1.00  5.73           C
ATOM     15  C   GLY A   2      20.158  16.376  20.213  1.00  5.60           C
ATOM     16  O   GLY A   2      18.979  16.512  20.569  1.00  5.67           O
ATOM     17  H   GLY A   2      20.302  18.818  21.430  1.00  5.27           H
ATOM     18  HA2 GLY A   2      22.073  17.118  20.298  1.00  5.73           H
ATOM     19  HA3 GLY A   2      21.423  16.953  21.728  1.00  5.73           H
END
"""

# SCY (GeoStd ideal coordinates): H1 (the dictionary's H), H2 (on N in the
# dictionary), H3
pdb_scy = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   SCY A   1      23.056 -11.316   8.334  1.00 20.00           N
HETATM    2  CA  SCY A   1      24.398 -11.863   8.288  1.00 20.00           C
HETATM    3  CB  SCY A   1      24.407 -13.293   8.816  1.00 20.00           C
HETATM    4  SG  SCY A   1      23.858 -13.382  10.547  1.00 20.00           S
HETATM    5  CD  SCY A   1      22.196 -13.986  10.444  1.00 20.00           C
HETATM    6  OCD SCY A   1      21.678 -14.373   9.429  1.00 20.00           O
HETATM    7  CE  SCY A   1      21.503 -13.943  11.778  1.00 20.00           C
HETATM    8  C   SCY A   1      24.985 -11.809   6.848  1.00 20.00           C
HETATM    9  O   SCY A   1      25.985 -12.521   6.638  1.00 20.00           O
HETATM   10  OXT SCY A   1      24.438 -11.027   6.048  1.00 20.00           O
HETATM   11  H1  SCY A   1      23.014 -10.536   7.691  1.00 20.00           H
HETATM   12  H2  SCY A   1      22.858 -10.953   9.257  1.00 20.00           H
HETATM   13  HA  SCY A   1      25.105 -11.286   8.907  1.00 20.00           H
HETATM   14  HB2 SCY A   1      23.756 -13.940   8.227  1.00 20.00           H
HETATM   15  HB3 SCY A   1      25.413 -13.707   8.795  1.00 20.00           H
HETATM   16  HE1 SCY A   1      20.730 -14.708  11.828  1.00 20.00           H
HETATM   17  HE2 SCY A   1      22.194 -14.069  12.610  1.00 20.00           H
HETATM   18  HE3 SCY A   1      21.024 -12.968  11.890  1.00 20.00           H
HETATM   19  H3  SCY A   1      22.472 -11.954   8.125  1.00 20.00           H
END
"""

# ALA-GLY: H1, H2, H3
pdb_ala = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      -1.405   9.032  -6.919  1.00  5.08           N
ATOM      2  CA  ALA A   1      -2.526   8.411  -7.691  1.00  4.71           C
ATOM      3  C   ALA A   1      -2.739   9.202  -8.960  1.00  5.06           C
ATOM      4  O   ALA A   1      -2.863  10.449  -8.918  1.00  5.50           O
ATOM      5  CB  ALA A   1      -3.787   8.434  -6.835  1.00  4.85           C
ATOM      6  H1  ALA A   1      -0.664   9.023  -7.412  1.00  5.08           H
ATOM      7  H2  ALA A   1      -1.615   9.873  -6.715  1.00  5.08           H
ATOM      8  H3  ALA A   1      -1.269   8.570  -6.170  1.00  5.08           H
ATOM      9  HA  ALA A   1      -2.328   7.490  -7.920  1.00  4.71           H
ATOM     10  HB1 ALA A   1      -4.004   9.354  -6.617  1.00  4.85           H
ATOM     11  HB2 ALA A   1      -4.515   8.032  -7.334  1.00  4.85           H
ATOM     12  HB3 ALA A   1      -3.626   7.931  -6.022  1.00  4.85           H
ATOM     13  N   GLY A   2      -2.773   8.494 -10.089  1.00  5.27           N
ATOM     14  CA  GLY A   2      -3.068   9.061 -11.398  1.00  5.73           C
ATOM     15  C   GLY A   2      -4.514   8.748 -11.734  1.00  5.60           C
ATOM     16  O   GLY A   2      -4.966   7.602 -11.604  1.00  5.67           O
ATOM     17  OXT GLY A   2      -5.200   9.700 -12.100  1.00  5.60           O
ATOM     18  H   GLY A   2      -2.622   7.648 -10.122  1.00  5.27           H
ATOM     19  HA2 GLY A   2      -2.939  10.022 -11.386  1.00  5.73           H
ATOM     20  HA3 GLY A   2      -2.488   8.673 -12.072  1.00  5.73           H
END
"""

# NPH-VAL as ready_set (old reduce) writes it: H and H2 on N; NPH's H3
# is on C3
pdb_nph_old_reduce = """\
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  N   NPH A   1      23.870  15.268 -50.490  1.00 25.06           N
HETATM    2  CA  NPH A   1      23.515  14.664 -49.210  1.00 22.94           C
HETATM    3  C   NPH A   1      22.061  14.206 -49.245  1.00 23.82           C
HETATM    4  O   NPH A   1      21.158  15.023 -49.423  1.00 24.62           O
HETATM    5  CB  NPH A   1      23.658  15.702 -48.113  1.00 21.44           C
HETATM    6  SG  NPH A   1      25.281  16.410 -47.839  1.00 31.29           S
HETATM    7  CD  NPH A   1      25.498  16.310 -46.059  0.50 33.70           C
HETATM    8  CE  NPH A   1      26.971  16.492 -45.820  0.50 39.92           C
HETATM    9  NZ  NPH A   1      27.411  17.709 -45.649  0.50 43.87           N
HETATM   10  OZ  NPH A   1      27.815  15.348 -45.806  0.50 43.26           O
HETATM   11  C10 NPH A   1      30.543  19.684 -45.064  0.50 54.23           C
HETATM   12  C1A NPH A   1      31.522  18.669 -45.046  0.50 55.19           C
HETATM   13  C2  NPH A   1      33.832  17.924 -44.852  0.50 55.79           C
HETATM   14  C3  NPH A   1      33.439  16.599 -45.020  0.50 53.96           C
HETATM   15  C4  NPH A   1      32.087  16.314 -45.202  0.50 53.40           C
HETATM   16  C4A NPH A   1      31.137  17.332 -45.219  0.50 52.42           C
HETATM   17  C5  NPH A   1      29.800  17.009 -45.411  0.50 51.26           C
HETATM   18  C6  NPH A   1      28.830  18.015 -45.436  0.50 51.14           C
HETATM   19  C6A NPH A   1      29.195  19.344 -45.262  0.50 52.54           C
HETATM   20  C7  NPH A   1      28.201  20.331 -45.288  0.50 53.28           C
HETATM   21  C8  NPH A   1      28.549  21.670 -45.111  0.50 54.83           C
HETATM   22  C9  NPH A   1      29.889  22.011 -44.911  0.50 55.29           C
HETATM   23  N1  NPH A   1      32.879  18.963 -44.861  0.50 57.78           N
HETATM   24  N10 NPH A   1      30.892  21.027 -44.889  0.50 53.87           N
HETATM   25  H   NPH A   1      24.112  16.093 -50.524  1.00 25.06           H
HETATM   26  H2  NPH A   1      23.837  14.799 -51.210  1.00 25.06           H
HETATM   27  H3  NPH A   1      34.108  15.897 -45.009  0.50 53.96           H
HETATM   28  HA  NPH A   1      24.099  13.908 -49.042  1.00 22.94           H
HETATM   29  HB2 NPH A   1      23.355  15.287 -47.290  1.00 21.44           H
HETATM   30  HB3 NPH A   1      23.062  16.441 -48.312  1.00 21.44           H
HETATM   31  HD2 NPH A   1      25.226  15.451 -45.699  0.50 33.70           H
HETATM   32  HD3 NPH A   1      25.032  17.028 -45.603  0.50 33.70           H
HETATM   33  HZ  NPH A   1      26.833  18.345 -45.684  0.50 43.87           H
HETATM   34  H4  NPH A   1      31.747  15.412 -45.311  0.50 53.40           H
HETATM   35  H5  NPH A   1      29.643  16.059 -45.528  0.50 51.26           H
HETATM   36  H7  NPH A   1      27.261  20.125 -45.413  0.50 53.28           H
HETATM   37  H8  NPH A   1      27.874  22.367 -45.116  0.50 54.83           H
HETATM   38  H9  NPH A   1      30.089  22.949 -44.766  0.50 55.29           H
HETATM   39  HC2 NPH A   1      34.778  18.099 -44.728  0.50 55.79           H
ATOM     40  N   VAL A   2      21.831  12.913 -49.037  1.00 28.67           N
ATOM     41  CA  VAL A   2      20.474  12.371 -49.062  1.00 33.96           C
ATOM     42  C   VAL A   2      20.069  11.729 -47.737  1.00 36.17           C
ATOM     43  O   VAL A   2      20.854  11.018 -47.109  1.00 40.12           O
ATOM     44  CB  VAL A   2      20.297  11.336 -50.209  1.00 39.91           C
ATOM     45  CG1 VAL A   2      18.831  10.915 -50.336  1.00 38.68           C
ATOM     46  CG2 VAL A   2      20.791  11.918 -51.524  1.00 43.67           C
ATOM     47  H   VAL A   2      22.442  12.328 -48.879  1.00 28.67           H
ATOM     48  HA  VAL A   2      19.848  13.087 -49.232  1.00 33.96           H
ATOM     49  HB  VAL A   2      20.823  10.546 -50.010  1.00 39.91           H
ATOM     50 HG11 VAL A   2      18.723  10.392 -51.146  1.00 38.68           H
ATOM     51 HG12 VAL A   2      18.585  10.381 -49.564  1.00 38.68           H
ATOM     52 HG13 VAL A   2      18.277  11.710 -50.381  1.00 38.68           H
ATOM     53 HG21 VAL A   2      20.576  11.302 -52.242  1.00 43.67           H
ATOM     54 HG22 VAL A   2      20.352  12.770 -51.676  1.00 43.67           H
ATOM     55 HG23 VAL A   2      21.751  12.045 -51.474  1.00 43.67           H
END
"""

_srv = []
def interpret(pdb_str):
  """residue names with their modifications; {atom name of residue 1:
  sorted bonded atom names}; {(name, name): N-CA ideal}"""
  if not _srv:
    _srv.append(mmtbx.monomer_library.server.server())
    _srv.append(mmtbx.monomer_library.server.ener_lib())
  p = pdb_interpretation.process(mon_lib_srv=_srv[0], ener_lib=_srv[1],
    pdb_inp=iotbx.pdb.input(source_info=None, lines=pdb_str),
    log=null_out(),
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  acp = p.all_chain_proxies
  names = set(mm.residue_name for mm in acp.all_monomer_mappings)
  atoms = acp.pdb_hierarchy.atoms()
  grm = p.geometry_restraints_manager(show_energies=False)
  bonded, n_ca = {}, []
  for prx in grm.pair_proxies(sites_cart=atoms.extract_xyz()
      ).bond_proxies.simple:
    a, b = [atoms[i] for i in prx.i_seqs]
    if a.parent().parent().resseq.strip() != '1': continue
    if b.parent().parent().resseq.strip() != '1': continue
    bonded.setdefault(a.name.strip(), []).append(b.name.strip())
    bonded.setdefault(b.name.strip(), []).append(a.name.strip())
    if sorted([a.name.strip(), b.name.strip()]) == ['CA', 'N']:
      n_ca.append(prx.distance_ideal)
  first = [a for a in atoms if a.parent().parent().resseq.strip() == '1']
  h = [a.name.strip() for a in first if a.element.strip() == 'H']
  return names, dict((k, sorted(v)) for k, v in bonded.items()), h, n_ca

def test_0a0():
  # the generated modification: every H on N restrained, N-CA as mod_NH3
  # (CA sp3)
  names, bonded, h, n_ca = interpret(pdb_0a0)
  assert '0A0%NH3-0A0' in names, names
  for x in h:
    assert x in bonded, (x, bonded)
  for x in ('HN1', 'HN2', 'HN3'):
    assert bonded[x] == ['N'], (x, bonded[x])
  assert len(n_ca) == 1 and abs(n_ca[0] - 1.491) < 1.e-6, n_ca

def test_2kz_carbon_h():
  # H1 and H3 on N: mod_NH3 is not applied (it would add H2 to N and take
  # CB's H2); H2 stays restrained to CB; H1 and H3 are unknown
  names, bonded, h, _ = interpret(pdb_2kz)
  assert '2KZ' in names and not [x for x in names if '%' in x and
    x.startswith('2KZ')], names
  assert bonded['H2'] == ['CB'], bonded['H2']
  assert 'H1' not in bonded and 'H3' not in bonded, bonded

def test_scy():
  # H2 on N in the dictionary: mod_NH3 still applies and replaces it
  names, bonded, h, _ = interpret(pdb_scy)
  assert 'SCY%NH3' in names, names
  for x in ('H1', 'H2', 'H3'):
    assert bonded.get(x) == ['N'], (x, bonded.get(x))

def test_standard_terminus():
  names, bonded, h, n_ca = interpret(pdb_ala)
  assert 'ALA%NH3' in names, names
  for x in h:
    assert x in bonded, (x, bonded)
  assert len(n_ca) == 1 and abs(n_ca[0] - 1.491) < 1.e-6, n_ca

def test_partial_and_sp3():
  # NPH needs two computed names (H2, HN3): one of them present is enough,
  # as for mod_NH3's H1-H3; CA sp3 by neighbour count, also for
  # mon_lib-format entries whose bonds are all "coval" (MSE)
  from mmtbx.monomer_library import terminal_nh3
  interpret(pdb_ala)
  srv = _srv[0]
  nph = srv.get_comp_comp_id_direct('NPH')
  assert terminal_nh3.extra_h_names(nph) == ['H2', 'HN3']
  for unexpected in (['H2'], ['HN3'], ['H2', 'HN3']):
    assert terminal_nh3.terminal_mod(nph, unexpected) is not None, \
      unexpected
  assert terminal_nh3.terminal_mod(nph, ['HXX']) is None
  for comp, sp3 in (('0A0', True), ('MSE', True), ('2KZ', False)):
    assert terminal_nh3.ca_sp3(srv.get_comp_comp_id_direct(comp)) == sp3, \
      comp

def n_h_values(srv, ener_lib):
  """{H name: (neutron ideal, esd)} of the N-H bonds of residue 1"""
  params = pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  p = pdb_interpretation.process(mon_lib_srv=srv, ener_lib=ener_lib,
    params=params, pdb_inp=iotbx.pdb.input(source_info=None, lines=pdb_0a0),
    log=null_out(),
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  atoms = p.all_chain_proxies.pdb_hierarchy.atoms()
  grm = p.geometry_restraints_manager(show_energies=False)
  result = {}
  for prx in grm.pair_proxies(sites_cart=atoms.extract_xyz()
      ).bond_proxies.simple:
    a, b = sorted([atoms[i] for i in prx.i_seqs],
      key=lambda x: x.name.strip() != 'N')
    if (a.parent().parent().resseq.strip() == '1' and a.name.strip() == 'N'
        and b.name.strip().startswith('HN')):
      result[b.name.strip()] = (round(prx.distance_ideal, 3),
        round(prx.weight**-0.5, 3))
  return result

def test_replacement_cif():
  # a CIF for 0A0 loaded into the same server after a first interpretation
  # (N-H esd 0.03, neutron 1.08): the generated HN3 takes the new values
  # with the dictionary's HN1 and HN2, as in a fresh server
  path = libtbx.env.find_in_repositories(
    relative_path='chem_data/geostd/0/data_0A0.cif', test=os.path.isfile)
  if path is None:
    print('Skipping test_replacement_cif: chem_data not available')
    return
  srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  assert n_h_values(srv, ener_lib) == dict((h, (1.04, 0.02))
    for h in ('HN1', 'HN2', 'HN3')), n_h_values(srv, ener_lib)
  text = open(path).read()
  old = 'single  0.890  0.020  1.040'
  assert text.count(old) == 2, text.count(old)
  cif = iotbx.cif.reader(input_string=text.replace(old,
    'single  0.890  0.030  1.080')).model()
  srv.process_cif_object(cif_object=cif, file_name='0A0_replacement.cif')
  values = n_h_values(srv, ener_lib)
  assert values == dict((h, (1.08, 0.03)) for h in ('HN1', 'HN2', 'HN3')), \
    values

def mod_n_h(srv, comp):
  """{N-H name: (x-ray, esd, neutron)} of comp with its generated
  modification applied"""
  mon = srv.get_comp_comp_id_direct(comp)
  mod = terminal_nh3.nh3_mod(mon, terminal_nh3.extra_h_names(mon),
    srv.mod_mod_id_dict['NH3'])
  new = mon.apply_mod(mod)
  result = {}
  for b in new.bond_list:
    ids = set([b.atom_id_1, b.atom_id_2])
    if 'N' in ids and (ids - set(['N'])) & set(
        terminal_nh3.n_h_names(mon) + terminal_nh3.extra_h_names(mon)):
      result[list(ids - set(['N']))[0]] = (b.value_dist, b.value_dist_esd,
        b.value_dist_neutron)
  return result

def test_n_h_values():
  # one value set for the three N-H: the dictionary's own where its N-H
  # agree (0A0), else mod_NH3's (TSY and 2KZ: 0.86 and 0.89; NPH: no
  # neutron value)
  srv = mmtbx.monomer_library.server.server()
  own = mod_n_h(srv, '0A0')
  assert sorted(own) == ['HN1', 'HN2', 'HN3'], own
  assert set(own.values()) == set([(0.89, 0.02, 1.04)]), own
  assert terminal_nh3.n_h_values(srv.get_comp_comp_id_direct('0A0'))[1]
  for comp in ('TSY', 'NPH', '2KZ'):
    values = mod_n_h(srv, comp)
    assert len(values) == 3, (comp, values)
    assert set(values.values()) == set([(0.89, 0.02, 1.04)]), (comp, values)
  mon = srv.get_comp_comp_id_direct('TSY')
  assert terminal_nh3.n_h_values(mon)[1] is False
  # one nuclear length in both columns is the dictionary's own; an x-ray
  # length repeated (NPH 0.86/0.86) is not
  mon = srv.get_comp_comp_id_direct('0A0')
  for b in mon.bond_list:
    if 'N' in (b.atom_id_1, b.atom_id_2) and set([b.atom_id_1,
        b.atom_id_2]) & set(terminal_nh3.n_h_names(mon)):
      b.value_dist = b.value_dist_neutron = 1.03
  assert terminal_nh3.n_h_values(mon) == ((1.03, 0.02, 1.03), True)
  for b in mon.bond_list:
    if 'N' in (b.atom_id_1, b.atom_id_2) and set([b.atom_id_1,
        b.atom_id_2]) & set(terminal_nh3.n_h_names(mon)):
      b.value_dist = b.value_dist_neutron = 0.86
  assert terminal_nh3.n_h_values(mon) == ((0.89, 0.02, 1.04), False)

def test_imine_and_bond_order():
  # an N=CA imine is not an amine terminus (HHI, ITR); 2KZ's N-CA (single,
  # 1.402, sp2 CA) is. Without a bond order (coval) the ideal length
  # decides: N_CA_SINGLE_MIN
  srv = mmtbx.monomer_library.server.server()
  for comp in ('HHI', 'ITR'):
    mon = srv.get_comp_comp_id_direct(comp)
    assert not terminal_nh3.n_ca_single(mon), comp
    assert not terminal_nh3.applies(mon), comp
  mon = srv.get_comp_comp_id_direct('2KZ')
  assert terminal_nh3.n_ca_single(mon) and terminal_nh3.applies(mon)
  mon = srv.get_comp_comp_id_direct('0A0')
  n_ca = [b for b in mon.bond_list if set([b.atom_id_1, b.atom_id_2]) ==
    set(['N', 'CA'])][0]
  n_ca.type = 'coval'
  for d, expected in ((1.30, False), (terminal_nh3.N_CA_SINGLE_MIN, True),
                      (1.458, True)):
    n_ca.value_dist = d
    assert terminal_nh3.n_ca_single(mon) == expected, (d, expected)
  n_ca.value_dist = None
  assert not terminal_nh3.n_ca_single(mon)

def test_computed_rule_name():
  # 0AA: the computed HN2 is also a name of the HN1/HN2 -> mod_NH2 rule;
  # NH3-0AA restrains it, never mod_NH2
  srv = mmtbx.monomer_library.server.server()
  mon = srv.get_comp_comp_id_direct('0AA')
  assert terminal_nh3.extra_h_names(mon) == ['HN2']
  nh3 = srv.mod_mod_id_dict['NH3']
  assert terminal_nh3.computed_names_mod(mon, ['HN2'], nh3) is not None
  assert terminal_nh3.computed_names_mod(mon, ['DN2'], nh3) is not None
  # a supplied rule name that is not computed: the rules as before
  assert terminal_nh3.computed_names_mod(mon, ['HN2', 'H1'], nh3) is None
  assert terminal_nh3.computed_names_mod(mon, ['OXT'], nh3) is None

def test_h_series():
  # dictionary N-H named H: the missing H continue the series H, H2, H3
  # (CCD, old reduce) where free, then the fallback names. NPH H3 is on
  # C3: H, H2, HN3. Old reduce's H2 is restrained to N, never by mod_NH3
  # (it would take C3's H3); as D (all_d), D and D2. 0G6 (H1 on another
  # atom): H2, H3
  interpret(pdb_ala)
  srv = _srv[0]
  assert terminal_nh3.extra_h_names(
    srv.get_comp_comp_id_direct('NPH')) == ['H2', 'HN3']
  assert terminal_nh3.extra_h_names(
    srv.get_comp_comp_id_direct('0G6')) == ['H2', 'H3']
  # 2KZ (H and H4 on N) is not a lone H: H5, as before
  assert terminal_nh3.extra_h_names(
    srv.get_comp_comp_id_direct('2KZ')) == ['H5']
  d = pdb_nph_old_reduce.replace(' H   NPH', ' D   NPH').replace(
    ' H2  NPH', ' D2  NPH')
  d = '\n'.join(l[:76] + ' D' + l[78:] if l[12:16] in (' D  ', ' D2 ')
    else l for l in d.split('\n'))
  for pdb_str, h, h2 in ((pdb_nph_old_reduce, 'H', 'H2'), (d, 'D', 'D2')):
    names, bonded, hs, _ = interpret(pdb_str)
    # every atom typed (old reduce's H2 stopped model.process())
    p = pdb_interpretation.process(mon_lib_srv=_srv[0], ener_lib=_srv[1],
      pdb_inp=iotbx.pdb.input(source_info=None, lines=pdb_str),
      log=null_out(),
      substitute_non_crystallographic_unit_cell_if_necessary=True)
    assert p.all_chain_proxies.fatal_problems_message() is None, \
      p.all_chain_proxies.fatal_problems_message()
    assert 'NPH%NH3-NPH' in names, names
    assert not [x for x in names if x.startswith('NPH%NH3') and
      x != 'NPH%NH3-NPH'], names
    for x in hs:
      assert x in bonded, (x, bonded)
    assert bonded[h] == ['N'] and bonded[h2] == ['N'], bonded
    assert bonded['H3'] == ['C3'], bonded['H3']

if __name__ == '__main__':
  test_0a0()
  test_2kz_carbon_h()
  test_scy()
  test_standard_terminus()
  test_partial_and_sp3()
  test_replacement_cif()
  test_n_h_values()
  test_imine_and_bond_order()
  test_computed_rule_name()
  test_h_series()
  print('OK')
