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
from mmtbx.monomer_library import pdb_interpretation

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
  # NPH needs two computed names (HN3, HNT): one of them present is enough,
  # as for mod_NH3's H1-H3; CA sp3 by neighbour count, also for
  # mon_lib-format entries whose bonds are all "coval" (MSE)
  from mmtbx.monomer_library import terminal_nh3
  interpret(pdb_ala)
  srv = _srv[0]
  nph = srv.get_comp_comp_id_direct('NPH')
  assert terminal_nh3.extra_h_names(nph) == ['HN3', 'HNT']
  for unexpected in (['HN3'], ['HNT'], ['HN3', 'HNT']):
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

if __name__ == '__main__':
  test_0a0()
  test_2kz_carbon_h()
  test_scy()
  test_standard_terminus()
  test_partial_and_sp3()
  test_replacement_cif()
  print('OK')
