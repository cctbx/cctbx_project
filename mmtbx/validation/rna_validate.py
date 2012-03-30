import os, sys
from iotbx import pdb
from mmtbx.monomer_library import rna_sugar_pucker_analysis
from iotbx.pdb import common_residue_names_get_class
import iotbx.phil
from mmtbx.monomer_library import pdb_interpretation
from mmtbx import monomer_library
from cctbx import geometry_restraints
from libtbx import easy_run
from libtbx.utils import Usage

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
  rna_validate {
    pdb = None
      .type = path
      .help = '''Enter a PDB file name'''
    outliers_only = False
      .type = bool
      .help = '''Only print rotamer outliers'''
    rna_sugar_pucker_analysis
      .short_caption = RNA sugar pucker analysis
      .style = box noauto auto_align menu_item parent_submenu:advanced
    {
      include scope mmtbx.monomer_library.rna_sugar_pucker_analysis.master_phil
    }
  }
  """,process_includes=True)

class rna_validate(object):

  def __init__(self):
    self.params = None
    self.rna_backbone_atoms = ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C1'",
                               "C3'", "O3'", "C2'", "O2'", "N1", "N9" ] #version 3.x naming

  def usage(self):
    return """
phenix.rna_validate file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  outliers_only=False   only print outliers

Example:

  phenix.rna_validate pdb=1u8d.pdb outliers_only=True

"""

  def run(self, args):
    if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
      raise Usage(self.usage())
    master_phil = get_master_phil()
    import iotbx.utils
    input_objects = iotbx.utils.process_command_line_inputs(
      args=args,
      master_phil=master_phil,
      input_types=("pdb",))
    work_phil = master_phil.fetch(sources=input_objects["phil"])
    work_params = work_phil.extract()
    if len(input_objects["pdb"]) != 1:
      raise Usage("phenix.rna_validate mypdb.pdb")
    file_obj = input_objects["pdb"][0]
    filename = file_obj.file_name
    self.params=work_params
    if filename and os.path.exists(filename):
      try:
        import iotbx.pdb
      except ImportError:
        print "iotbx not loaded"
        return
      pdb_io = iotbx.pdb.input(filename)
    else:
      print "Please enter a file name"
      return
    self.analyze_pdb(pdb_io=pdb_io)
    self.print_results()

  def analyze_pdb(self,
                  params=None,
                  pdb_io=None,
                  processed_pdb_file=None,
                  outliers_only=False,
                  show_errors = False,
                  out=sys.stdout):
    if processed_pdb_file is None:
      mon_lib_srv = monomer_library.server.server()
      ener_lib = monomer_library.server.ener_lib()
      self.processed_pdb_file = pdb_interpretation.process(
        mon_lib_srv=mon_lib_srv,
        ener_lib=ener_lib,
        pdb_inp=pdb_io,
        substitute_non_crystallographic_unit_cell_if_necessary=True)
    else:
      self.processed_pdb_file=processed_pdb_file
    self.pdb_hierarchy = self.processed_pdb_file.all_chain_proxies.pdb_hierarchy
    self.pucker_outliers = self.pucker_evaluate(
                           hierarchy=self.pdb_hierarchy)
    self.bond_outliers, self.angle_outliers = self.bond_and_angle_evaluate()
    self.suite_outliers = self.run_suitename()

  def print_results(self):
    print "RNA Validation"
    print "-----------------------------------------------"
    print "Pucker Outliers:"
    print "#residue:delta_angle:is_delta_outlier:epsilon_angle:is_epsilon_outler"
    for outlier in self.pucker_outliers:
      if outlier[1][1]:
        is_delta_outlier="yes"
      else:
        is_delta_outlier="no"
      if outlier[1][3]:
        is_epsilon_outlier="yes"
      else:
        is_epsilon_outlier="no"
      if outlier[1][0] == None:
        delta_angle = -1.0
      else:
        delta_angle = outlier[1][0]
      if outlier[1][2] == None:
        ep_angle = -1.0
      else:
        ep_angle = outlier[1][2]
      print "%s:%.3f:%s:%.3f:%s" % (outlier[0],
                                      delta_angle,
                                      is_delta_outlier,
                                      ep_angle,
                                      is_epsilon_outlier)
    print "\n-----------------------------------------------"
    print "Bond Length Outliers:"
    print "#residue:atom_1:atom_2:num_sigmas"
    for outlier in self.bond_outliers:
      for pair in outlier[1]:
        print "%s:%s:%s:%.3f" % (outlier[0], pair[0], pair[1], pair[2])
    print "\n-----------------------------------------------"
    print "Angle Outliers:"
    print "#residue:atom_1:atom_2:atom_3:num_sigmas"
    for outlier in self.angle_outliers:
      for pair in outlier[1]:
        print "%s:%s:%s:%s:%.3f" % (outlier[0], pair[0], pair[1], pair[2], pair[3])
    print "\n-----------------------------------------------"
    print "Suite Outliers:"
    print "#suiteID:triaged_angle"
    for outlier in self.suite_outliers:
      print "%s:%s" % (outlier[0], outlier[1])

  def run_suitename(self):
    suite_outliers = []
    suitename = "phenix.suitename -report -"
    backbone_dihedrals = self.get_rna_backbone_dihedrals()
    suitename_out = easy_run.fully_buffered(suitename,
                         stdin_lines=backbone_dihedrals).stdout_lines
    for line in suitename_out:
      if '!!' in line:
        temp = line.split(":")
        key = ' '+temp[5][0:3]+temp[2]+temp[3]+temp[4]
        temp2 = temp[5].split(" ")
        suite_outliers.append([key,temp2[len(temp2)-1]])
    return suite_outliers


  def match_dihedral_to_name(self, atoms):
    name = None
    alpha = ["O3'","P","O5'","C5'"]
    beta = ["P","O5'","C5'","C4'"]
    gamma = ["O5'","C5'","C4'","C3'"]
    delta = ["C5'","C4'","C3'","O3'"]
    epsilon = ["C4'","C3'","O3'","P"]
    zeta = ["C3'","O3'","P","O5'"]
    if atoms == alpha:
      name = "alpha"
    elif atoms == beta:
      name = "beta"
    elif atoms == gamma:
      name = "gamma"
    elif atoms == delta:
      name = "delta"
    elif atoms == epsilon:
      name = "epsilon"
    elif atoms == zeta:
      name = "zeta"
    return name

  def get_rna_backbone_dihedrals(self):
    bb_dihedrals = {}
    formatted_out = []
    sites_cart = self.processed_pdb_file.all_chain_proxies.sites_cart
    i_seq_name_hash = self.build_name_hash(
                      pdb_hierarchy=self.pdb_hierarchy)
    geometry = self.processed_pdb_file.geometry_restraints_manager()
    dihedral_proxies = geometry.dihedral_proxies
    for dp in dihedral_proxies:
      atoms = []
      for i in dp.i_seqs:
        atoms.append(i_seq_name_hash[i][0:4].strip())
      name = self.match_dihedral_to_name(atoms=atoms)
      if name is not None:
        restraint = geometry_restraints.dihedral(
                                                 sites_cart=sites_cart,
                                                 proxy=dp)
        key = i_seq_name_hash[dp.i_seqs[1]][5:]
        try:
          bb_dihedrals[key][name]=restraint.angle_model
        except Exception:
          bb_dihedrals[key] = {}
          bb_dihedrals[key][name]=restraint.angle_model
    for key in bb_dihedrals.keys():
      resname = key[0:3]
      chainID = key[3:5]
      resnum = key[5:9]
      i_code = key[9:]
      try:
        alpha = "%.3f" % bb_dihedrals[key]['alpha']
      except Exception:
        alpha = '__?__'
      try:
        beta = "%.3f" % bb_dihedrals[key]['beta']
      except Exception:
        beta = '__?__'
      try:
        gamma = "%.3f" % bb_dihedrals[key]['gamma']
      except Exception:
        gamma = '__?__'
      try:
        delta = "%.3f" % bb_dihedrals[key]['delta']
      except Exception:
        delta = '__?__'
      try:
        epsilon = "%.3f" % bb_dihedrals[key]['epsilon']
      except Exception:
        epsilon = '__?__'
      try:
        zeta = "%.3f" % bb_dihedrals[key]['zeta']
      except Exception:
        zeta = '__?__'
      eval = "%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s" \
             % (" ",
                "1",
                chainID,
                resnum,
                i_code,
                resname,
                alpha,
                beta,
                gamma,
                delta,
                epsilon,
                zeta)
      formatted_out.append(eval)
    formatted_out.sort()
    backbone_dihedrals = ""
    for line in formatted_out:
      backbone_dihedrals += line+'\n'
    return backbone_dihedrals

  def bond_and_angle_evaluate(self):
    bond_hash = {}
    angle_hash = {}
    bond_outliers = []
    angle_outliers = []
    geometry = self.processed_pdb_file.geometry_restraints_manager()
    assert (geometry is not None)
    flags = geometry_restraints.flags.flags(default=True)
    i_seq_name_hash = self.build_name_hash(
                      pdb_hierarchy=self.processed_pdb_file.all_chain_proxies.pdb_hierarchy)
    pair_proxies = geometry.pair_proxies(flags=flags,
                            sites_cart=self.processed_pdb_file.all_chain_proxies.sites_cart)
    bond_proxies = pair_proxies.bond_proxies
    sites_cart=self.processed_pdb_file.all_chain_proxies.sites_cart
    for proxy in bond_proxies.simple:
      restraint = geometry_restraints.bond(
                                           sites_cart=sites_cart,
                                           proxy=proxy)
      atom1 = i_seq_name_hash[proxy.i_seqs[0]][0:4]
      atom2 = i_seq_name_hash[proxy.i_seqs[1]][0:4]
      residue = i_seq_name_hash[proxy.i_seqs[0]][4:]
      if (atom1.strip() not in self.rna_backbone_atoms or
          atom2.strip() not in self.rna_backbone_atoms):
        continue
      sigma = (1/restraint.weight)**(.5)
      num_sigmas = restraint.delta / sigma
      if abs(num_sigmas) >= 4:
        try:
          bond_hash[residue].append((atom1,atom2,num_sigmas))
        except Exception:
          bond_hash[residue] = []
          bond_hash[residue].append((atom1,atom2,num_sigmas))
    for proxy in geometry.angle_proxies:
      restraint = geometry_restraints.angle(
                                           sites_cart=sites_cart,
                                           proxy=proxy)
      atom1 = i_seq_name_hash[proxy.i_seqs[0]][0:4]
      atom2 = i_seq_name_hash[proxy.i_seqs[1]][0:4]
      atom3 = i_seq_name_hash[proxy.i_seqs[2]][0:4]
      residue = i_seq_name_hash[proxy.i_seqs[0]][4:]
      if (atom1.strip() not in self.rna_backbone_atoms or
          atom2.strip() not in self.rna_backbone_atoms or
          atom3.strip() not in self.rna_backbone_atoms):
        continue
      # double sigma to account for artificially small sigmas used for
      # weighting
      sigma = ((1/restraint.weight)**(.5))*2
      num_sigmas = restraint.delta / sigma
      if abs(num_sigmas) >= 4:
        try:
          angle_hash[residue].append((atom1,atom2,atom3,num_sigmas))
        except Exception:
          angle_hash[residue] = []
          angle_hash[residue].append((atom1,atom2,atom3,num_sigmas))
    for key in bond_hash.keys():
      bond_outliers.append([key, bond_hash[key]])
    for key in angle_hash.keys():
      angle_outliers.append([key, angle_hash[key]])
    return bond_outliers, angle_outliers

  def pucker_evaluate(self, hierarchy):
    if self.params is None:
      master_phil = get_master_phil()
      self.params = master_phil.extract()
    params = self.params.rna_validate.rna_sugar_pucker_analysis
    self.pucker_states = []
    self.pucker_perp_xyz = {}
    self.pucker_dist = {}
    outliers = []
    from iotbx.pdb.rna_dna_detection import residue_analysis
    for model in hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          altloc = conformer.altloc
          if altloc == "":
            altloc = " "
          residues = conformer.residues()
          for i_residue,residue in enumerate(residues):
            def _get_next_residue():
              j = i_residue + 1
              if (j == len(residues)): return None
              return residues[j]
            ra1 = residue_analysis(
              residue_atoms=residue.atoms(),
              distance_tolerance=params.bond_detection_distance_tolerance)
            if (ra1.problems is not None): continue
            if (not ra1.is_rna): continue
            residue_2_p_atom = None
            next_pdb_residue = _get_next_residue()
            if (next_pdb_residue is not None):
              residue_2_p_atom = next_pdb_residue.find_atom_by(name=" P  ")
            if common_residue_names_get_class(residue.resname) != "common_rna_dna":
              continue
            ana = rna_sugar_pucker_analysis.evaluate(
              params=params,
              residue_1_deoxy_ribo_atom_dict=ra1.deoxy_ribo_atom_dict,
              residue_1_c1p_outbound_atom=ra1.c1p_outbound_atom,
              residue_2_p_atom=residue_2_p_atom)
            self.pucker_states.append(ana)
            if ana.is_delta_outlier or ana.is_epsilon_outlier:
              key = residue.id_str()[8:-1]
              key = altloc+key
              outliers.append([key,
                               [ana.delta,
                               ana.is_delta_outlier,
                               ana.epsilon,
                               ana.is_epsilon_outlier]])
              self.pucker_perp_xyz[key] = [ana.p_perp_xyz, ana.o3p_perp_xyz]
              self.pucker_dist[key] = [ana.p_distance_c1p_outbound_line,
                                       ana.o3p_distance_c1p_outbound_line]
    return outliers

  def build_name_hash(self, pdb_hierarchy):
    i_seq_name_hash = dict()
    for atom in pdb_hierarchy.atoms():
      i_seq_name_hash[atom.i_seq]=atom.pdb_label_columns()
    return i_seq_name_hash

class rna_pucker_ref(object):

  def __init__(self):
    #residue 21 is 3', residue 22 is 2', residue 3 included for P only
    self.sample_bases = """
CRYST1  132.298   35.250   42.225  90.00  90.95  90.00 C 1 2 1       4
ATOM    132  P     A A  21       8.804   4.707  -0.950  1.00 21.07           P
ATOM    133  OP1   A A  21       8.958   5.041  -2.388  1.00 19.56           O
ATOM    134  OP2   A A  21       9.726   3.742  -0.339  1.00 17.87           O
ATOM    135  O5'   A A  21       8.914   6.057  -0.104  1.00 18.44           O
ATOM    136  C5'   A A  21       8.311   7.240  -0.587  1.00 16.49           C
ATOM    137  C4'   A A  21       8.458   8.398   0.384  1.00 15.52           C
ATOM    138  O4'   A A  21       7.781   8.097   1.635  1.00 16.13           O
ATOM    139  C3'   A A  21       9.885   8.660   0.849  1.00 16.12           C
ATOM    140  O3'   A A  21      10.578   9.415  -0.119  1.00 15.84           O
ATOM    141  C2'   A A  21       9.670   9.462   2.123  1.00 14.81           C
ATOM    142  O2'   A A  21       9.390  10.824   1.853  1.00 16.92           O
ATOM    143  C1'   A A  21       8.420   8.794   2.696  1.00 16.03           C
ATOM    144  N9    A A  21       8.745   7.898   3.805  1.00 14.49           N
ATOM    145  C8    A A  21       9.200   6.611   3.807  1.00 14.35           C
ATOM    146  N7    A A  21       9.500   6.161   5.002  1.00 14.14           N
ATOM    147  C5    A A  21       9.192   7.227   5.846  1.00 13.55           C
ATOM    148  C6    A A  21       9.288   7.400   7.255  1.00 13.46           C
ATOM    149  N6    A A  21       9.760   6.467   8.083  1.00 13.17           N
ATOM    150  N1    A A  21       8.881   8.588   7.778  1.00 14.41           N
ATOM    151  C2    A A  21       8.419   9.532   6.935  1.00 14.69           C
ATOM    152  N3    A A  21       8.294   9.483   5.600  1.00 16.23           N
ATOM    153  C4    A A  21       8.702   8.288   5.122  1.00 13.95           C
ATOM    154  P     U A  22      11.833   8.778  -0.860  1.00 17.19           P
ATOM    155  OP1   U A  22      12.143   9.766  -1.930  1.00 14.62           O
ATOM    156  OP2   U A  22      11.572   7.360  -1.223  1.00 16.11           O
ATOM    157  O5'   U A  22      12.955   8.795   0.274  1.00 17.34           O
ATOM    158  C5'   U A  22      14.086   7.879   0.260  1.00 15.80           C
ATOM    159  C4'   U A  22      14.668   7.761   1.654  1.00 13.76           C
ATOM    160  O4'   U A  22      14.980   9.104   2.130  1.00 13.53           O
ATOM    161  C3'   U A  22      13.575   7.229   2.579  1.00 14.34           C
ATOM    162  O3'   U A  22      14.069   6.400   3.635  1.00 16.58           O
ATOM    163  C2'   U A  22      12.950   8.479   3.204  1.00 14.31           C
ATOM    164  O2'   U A  22      12.419   8.301   4.510  1.00 12.30           O
ATOM    165  C1'   U A  22      14.140   9.444   3.207  1.00 13.55           C
ATOM    166  N1    U A  22      13.832  10.876   3.147  1.00 14.29           N
ATOM    167  C2    U A  22      14.055  11.619   4.286  1.00 13.16           C
ATOM    168  O2    U A  22      14.534  11.153   5.327  1.00 15.12           O
ATOM    169  N3    U A  22      13.706  12.925   4.186  1.00 13.49           N
ATOM    170  C4    U A  22      13.172  13.569   3.099  1.00 13.63           C
ATOM    171  O4    U A  22      12.789  14.728   3.228  1.00 12.03           O
ATOM    172  C5    U A  22      12.990  12.747   1.948  1.00 12.63           C
ATOM    173  C6    U A  22      13.320  11.453   2.009  1.00 13.42           C
ATOM    174  P     A A  23      14.502   4.869   3.397  1.00 16.90           P
ATOM    175  OP1   A A  23      13.773   4.291   2.214  1.00 17.04           O
ATOM    176  OP2   A A  23      14.351   4.224   4.726  1.00 16.03           O
ATOM    177  O5'   A A  23      16.057   4.949   3.092  1.00 16.99           O
ATOM    178  C5'   A A  23      16.937   5.684   3.940  1.00 16.50           C
ATOM    179  C4'   A A  23      18.220   5.970   3.200  1.00 16.03           C
ATOM    180  O4'   A A  23      17.896   6.559   1.914  1.00 15.15           O
ATOM    181  C3'   A A  23      19.123   7.000   3.865  1.00 15.88           C
ATOM    182  O3'   A A  23      20.007   6.350   4.779  1.00 16.04           O
ATOM    183  C2'   A A  23      19.932   7.535   2.689  1.00 14.37           C
ATOM    184  O2'   A A  23      21.038   6.689   2.437  1.00 14.07           O
ATOM    185  C1'   A A  23      18.924   7.447   1.539  1.00 15.13           C
ATOM    186  N9    A A  23      18.313   8.717   1.157  1.00 14.71           N
ATOM    187  C8    A A  23      18.027   9.134  -0.113  1.00 14.98           C
ATOM    188  N7    A A  23      17.513  10.334  -0.171  1.00 15.52           N
ATOM    189  C5    A A  23      17.440  10.727   1.159  1.00 14.12           C
ATOM    190  C6    A A  23      16.984  11.892   1.762  1.00 13.27           C
ATOM    191  N6    A A  23      16.517  12.944   1.074  1.00 13.02           N
ATOM    192  N1    A A  23      17.029  11.962   3.109  1.00 13.36           N
ATOM    193  C2    A A  23      17.522  10.921   3.791  1.00 12.32           C
ATOM    194  N3    A A  23      17.996   9.768   3.328  1.00 12.80           N
ATOM    195  C4    A A  23      17.924   9.737   1.988  1.00 13.37           C
"""
    self.rna_backbone_atoms = ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C1'",
                               "C3'", "O3'", "C2'", "O2'"] #version 3.x naming
    self.get_rna_pucker_ref()

  def get_rna_pucker_ref(self, test_pucker=False):
    flags = geometry_restraints.flags.flags(default=True)
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    self.bond_dict = {}
    self.angle_dict = {}
    self.processed_pdb_file = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=self.sample_bases)
    self.geometry = self.processed_pdb_file.geometry_restraints_manager()
    #confirm puckers are correct
    if test_pucker:
      r = rna_validate()
      r.pucker_evaluate(hierarchy=self.processed_pdb_file.all_chain_proxies.pdb_hierarchy)
      assert not r.pucker_states[0].is_2p
      assert r.pucker_states[1].is_2p
    i_seq_name_hash = self.build_name_hash(
                      pdb_hierarchy=self.processed_pdb_file.all_chain_proxies.pdb_hierarchy)
    pair_proxies = self.geometry.pair_proxies(flags=flags,
                                 sites_cart=self.processed_pdb_file.all_chain_proxies.sites_cart)
    bond_proxies = pair_proxies.bond_proxies
    for bond in bond_proxies.simple:
      atom1 = i_seq_name_hash[bond.i_seqs[0]][0:4]
      atom2 = i_seq_name_hash[bond.i_seqs[1]][0:4]
      if (atom1.strip() not in self.rna_backbone_atoms and
          atom2.strip() not in self.rna_backbone_atoms):
        continue
      key = atom1+atom2
      sigma = (1/bond.weight)**(.5)
      try:
        self.bond_dict[key].append((bond.distance_ideal, sigma))
      except Exception:
        self.bond_dict[key] = []
        self.bond_dict[key].append((bond.distance_ideal, sigma))
    for angle in self.geometry.angle_proxies:
      atom1 = i_seq_name_hash[angle.i_seqs[0]][0:4]
      atom2 = i_seq_name_hash[angle.i_seqs[1]][0:4]
      atom3 = i_seq_name_hash[angle.i_seqs[2]][0:4]
      if (atom1.strip() not in self.rna_backbone_atoms and
          atom2.strip() not in self.rna_backbone_atoms and
          atom3.strip() not in self.rna_backbone_atoms):
        continue
      key = atom1+atom2+atom3
      sigma = (1/angle.weight)**(.5)
      try:
        self.angle_dict[key].append((angle.angle_ideal, sigma))
      except Exception:
        self.angle_dict[key] = []
        self.angle_dict[key].append((angle.angle_ideal, sigma))

  def build_name_hash(self, pdb_hierarchy):
    i_seq_name_hash = dict()
    for atom in pdb_hierarchy.atoms():
      i_seq_name_hash[atom.i_seq]=atom.pdb_label_columns()
    return i_seq_name_hash
