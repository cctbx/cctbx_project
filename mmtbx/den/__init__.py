import iotbx.phil
from cctbx.array_family import flex
from libtbx import easy_pickle
import sys

den_params = iotbx.phil.parse("""
 reference_file = None
   .type = path
   .optional = true
   .short_caption = DEN reference model
   .style = file_type:pdb input_file
 gamma = 0.5
   .type = float
 kappa = 0.1
   .type = float
 weight = 1.0
   .type = float
 optimize = False
   .type = bool
   .short_caption = Optimize DEN parameters
   .help = If selected, will run a grid search to determine optimal values \
    of the weight and gamma parameters for DEN refinement.  This is very \
    slow, but highly recommended, and can be parallelized across multiple \
    CPU cores.
 opt_gamma_values = 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
   .type = floats
   .short_caption = Gamma values for optimization
 opt_weight_values = 0.5, 1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0
   .type = floats
   .short_caption = Weight values for optimization
 num_cycles = 8
   .type = int
   .short_caption = Number of cycles
 refine_adp = True
   .type = bool
   .short_caption = Refine B-factors
 final_refinement_cycle = False
   .type = bool
 verbose = False
   .type = bool
 annealing_type = *torsion \
                   cartesian
   .type = choice(multi=False)
   .help = select strategy to apply DEN restraints
 output_kinemage = False
   .type = bool
   .help = output kinemage representation of starting DEN restraints
 restraint_network {
  lower_distance_cutoff = 3.0
    .type = float
  upper_distance_cutoff = 15.0
    .type = float
  sequence_separation_low = 0
    .type = int
  sequence_separation_limit = 10
    .type = int
  ndistance_ratio = 1.0
    .type = float
  export_den_pairs = False
    .type = bool
    .expert_level = 3
  den_network_file = None
    .type = path
    .optional = True
    .expert_level = 3
 }
""")

den_params_development = iotbx.phil.parse("""
 scale= 150.0
   .type = float
 relax_ncycle = 0
   .type = int
 post_ncycle = 0
   .type = int
 minimum_start = True
   .type = bool
 exponent = *2 4
   .type = choice(multi=False)
 atom_select = None
   .type = atom_selection
   .multiple = True
   .optional = True
""")

class den_restraints(object):

  def __init__(self,
               pdb_hierarchy,
               params,
               pdb_hierarchy_ref=None,
               log=None):
    if(log is None): log = sys.stdout
    self.log = log
    if len(pdb_hierarchy.models()) > 1:
      raise Sorry("More than one model in input model. DEN refinement "+
                  "is only available for a single model.")
    if pdb_hierarchy_ref is not None:
      if len(pdb_hierarchy_ref.models()) > 1:
        raise Sorry("More than one model in reference model. "+
                    "DEN refinement "+
                    "is only available for a single model.")
    self.pdb_hierarchy = pdb_hierarchy
    if pdb_hierarchy_ref is None:
      print >> self.log, "No input DEN reference model...restraining model "+ \
        "to starting structure"
      self.pdb_hierarchy_ref = pdb_hierarchy
    else:
      self.pdb_hierarchy_ref = pdb_hierarchy_ref
    from mmtbx.torsion_restraints import utils
    import boost.python
    self.ext = boost.python.import_ext("mmtbx_den_restraints_ext")
    self.params = params
    self.kappa = params.kappa
    self.gamma = params.gamma
    self.weight = params.weight
    self.num_cycles = params.num_cycles
    self.annealing_type = params.annealing_type
    self.ndistance_ratio = \
      params.restraint_network.ndistance_ratio
    self.lower_distance_cutoff = \
      params.restraint_network.lower_distance_cutoff
    self.upper_distance_cutoff = \
      params.restraint_network.upper_distance_cutoff
    self.sequence_separation_low = \
      params.restraint_network.sequence_separation_low
    self.sequence_separation_limit = \
      params.restraint_network.sequence_separation_limit
    self.den_network_file = \
      params.restraint_network.den_network_file
    self.export_den_pairs = \
      params.restraint_network.export_den_pairs
    self.den_proxies = None

    self.atoms_per_chain = \
      self.count_atoms_per_chain(pdb_hierarchy=pdb_hierarchy)
    self.atoms_per_chain_ref = \
      self.count_atoms_per_chain(pdb_hierarchy=self.pdb_hierarchy_ref)
    self.resid_hash_ref = \
      utils.build_resid_hash(pdb_hierarchy=self.pdb_hierarchy_ref)
    self.i_seq_hash = \
      utils.build_i_seq_hash(pdb_hierarchy=pdb_hierarchy)
    self.i_seq_hash_ref = \
      utils.build_i_seq_hash(pdb_hierarchy=self.pdb_hierarchy_ref)
    self.name_hash = \
      utils.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    self.name_hash_ref = \
      utils.build_name_hash(pdb_hierarchy=self.pdb_hierarchy_ref)
    self.ref_atom_pairs, self.ref_distance_hash = \
      self.find_atom_pairs(pdb_hierarchy=self.pdb_hierarchy_ref,
                           resid_hash=self.resid_hash_ref)
    self.remove_non_matching_pairs()
    if self.den_network_file is not None:
      self.den_atom_pairs = self.load_den_network()
    else:
      self.random_ref_atom_pairs = \
        self.select_random_den_restraints()
      self.den_atom_pairs = self.get_den_atom_pairs()
    if self.export_den_pairs:
      self.dump_den_network()
    self.build_den_restraints()

  def find_atom_pairs(self, pdb_hierarchy, resid_hash):
    print >> self.log, "finding DEN atom pairs..."
    atom_pairs = {}
    distance_hash = {}
    atom_pairs_test = {}
    distance_hash_test = {}
    #only supports first model
    low_dist_sq = self.lower_distance_cutoff**2
    high_dist_sq = self.upper_distance_cutoff**2
    residue_range = \
      self.sequence_separation_limit - self.sequence_separation_low
    for chain in pdb_hierarchy.models()[0].chains():
      found_conformer = None
      for conformer in chain.conformers():
        if not conformer.is_protein() and not conformer.is_na():
          continue
        elif found_conformer is None:
          found_conformer = conformer
        else:
          print >> self.log, "warning, multiple conformers found, using first"
      if found_conformer is not None:
        atom_pairs[chain.id] = []
        atom_pairs_test[chain.id] = []
        for i, res1 in enumerate(found_conformer.residues()):
          for res2 in found_conformer.residues()[i:i+residue_range+1]:
            separation = res2.resseq_as_int() - res1.resseq_as_int()
            if separation < self.sequence_separation_low or \
               separation > self.sequence_separation_limit:
              continue
            for j, atom1 in enumerate(res1.atoms()):
              for atom2 in res2.atoms():
                if atom2.i_seq <= atom1.i_seq:
                  continue
                dist = distance_squared(atom1.xyz, atom2.xyz)
                if dist >= low_dist_sq and \
                   dist <= high_dist_sq:
                  atom_pairs[chain.id].append( (atom1.i_seq, atom2.i_seq) )
                  distance_hash[(atom1.i_seq, atom2.i_seq)] = \
                    (dist**(0.5))
    return atom_pairs, distance_hash

  # remove any pairs of reference model atoms that do not
  # have matching atom pairs in the working model
  def remove_non_matching_pairs(self):
    print >> self.log, "removing non-matching pairs..."
    temp_atom_pairs = {}
    #temp_distance_hash = {}
    for chain in self.ref_atom_pairs.keys():
      temp_atom_pairs[chain] = []
      for i, pair in enumerate(self.ref_atom_pairs[chain]):
        ref_atom1 = self.name_hash_ref[pair[0]]
        ref_atom2 = self.name_hash_ref[pair[1]]
        model_atom1 = self.i_seq_hash.get(ref_atom1)
        model_atom2 = self.i_seq_hash.get(ref_atom2)
        if model_atom1 != None and model_atom2 != None:
          temp_atom_pairs[chain].append(pair)
    self.ref_atom_pairs = temp_atom_pairs

  def count_atoms_per_chain(self, pdb_hierarchy):
    atoms_per_chain = {}
    for chain in pdb_hierarchy.models()[0].chains():
      found_conformer = False
      for conformer in chain.conformers():
        if not conformer.is_protein() and not conformer.is_na():
          continue
        else:
          found_conformer = True
      if found_conformer:
        atoms_per_chain[chain.id] = chain.atoms_size()
    return atoms_per_chain

  def select_random_den_restraints(self):
    from cctbx.array_family import flex
    print >> self.log, "selecting random DEN restraints..."
    random_pairs = {}
    for chain in self.ref_atom_pairs.keys():
      random_pairs[chain] = []
      pair_list_size = len(self.ref_atom_pairs[chain])
      num_restraints = round(self.atoms_per_chain_ref[chain] *
                             self.ndistance_ratio)
      if num_restraints > pair_list_size:
        num_restraints = pair_list_size
      random_selection = \
        flex.random_selection(pair_list_size, int(num_restraints))
      for i in random_selection:
          random_pairs[chain].append(self.ref_atom_pairs[chain][i])
    return random_pairs

  def dump_den_network(self):
    den_dump = {}
    self.get_selection_strings()
    #print self.selection_string_hash
    #print self.random_ref_atom_pairs
    for chain in self.den_atom_pairs.keys():
      den_dump[chain] = []
      for pair in self.den_atom_pairs[chain]:
        i_seq_1 = pair[0]
        i_seq_2 = pair[1]
        select_1 = self.selection_string_hash[i_seq_1]
        select_2 = self.selection_string_hash[i_seq_2]
        dump_pair = (select_1, select_2)
        den_dump[chain].append(dump_pair)
    output_prefix = "den"
    easy_pickle.dump(
      "%s.pkl"%output_prefix,
      den_dump)

  def load_den_network(self):
    den_atom_pairs = {}
    network_pairs = easy_pickle.load(
      self.den_network_file)
    #check for current model compatibility
    sel_cache = self.pdb_hierarchy.atom_selection_cache()
    for chain in network_pairs.keys():
      den_atom_pairs[chain] = []
      for pair in network_pairs[chain]:
        string = "(%s) or (%s)" % (pair[0], pair[1])
        iselection = sel_cache.selection(string=string).iselection()
        if iselection.size() != 2:
          raise Sorry(
            "input DEN network does not match current model")
        den_atom_pairs[chain].append(iselection)
    return den_atom_pairs

  def get_selection_strings(self):
    selection_string_hash = {}
    atom_labels = list(self.pdb_hierarchy.atoms_with_labels())
    segids = flex.std_string([ a.segid for a in atom_labels ])
    for a in atom_labels:
      chain = a.chain_id
      resid = a.resid()
      resname = a.resname
      atomname = a.name
      altloc = a.altloc
      segid = a.segid
      selection_string = \
        "name '%s' and resname '%s' and chain '%s' and resid '%s'" % \
        (atomname, resname, chain, resid) + \
        " and segid '%s'" % (segid)
      if altloc != "":
        selection_string += " and altid '%s'" % altloc
      selection_string_hash[a.i_seq] = selection_string
    self.selection_string_hash = selection_string_hash

  def get_den_atom_pairs(self):
    den_atom_pairs = {}
    for chain in self.random_ref_atom_pairs.keys():
      den_atom_pairs[chain] = []
      for pair in self.random_ref_atom_pairs[chain]:
        i_seq_a = self.i_seq_hash[self.name_hash_ref[pair[0]]]
        i_seq_b = self.i_seq_hash[self.name_hash_ref[pair[1]]]
        i_seqs = flex.size_t([i_seq_a, i_seq_b])
        den_atom_pairs[chain].append(i_seqs)
    return den_atom_pairs

  def build_den_restraints(self):
    print >> self.log, "building DEN restraints..."
    den_proxies = self.ext.shared_den_simple_proxy()
    for chain in self.den_atom_pairs.keys():
      for pair in self.den_atom_pairs[chain]:
        i_seq_a = self.i_seq_hash_ref[self.name_hash[pair[0]]]
        i_seq_b = self.i_seq_hash_ref[self.name_hash[pair[1]]]
        distance_ideal = self.ref_distance_hash[ (i_seq_a, i_seq_b) ]
        i_seqs = tuple(pair)
        proxy = self.ext.den_simple_proxy(
          i_seqs=i_seqs,
          eq_distance=distance_ideal,
          eq_distance_start=distance_ideal,
          weight=self.weight)
        den_proxies.append(proxy)
    self.den_proxies = den_proxies

  def target_and_gradients(self,
                           sites_cart,
                           gradient_array):
    return self.ext.den_simple_residual_sum(
      sites_cart,
      self.den_proxies,
      gradient_array,
      self.weight)

  def update_eq_distances(self,
                          sites_cart):
    self.ext.den_update_eq_distances(sites_cart,
                            self.den_proxies,
                            self.gamma,
                            self.kappa)

  def get_optimization_grid(self):
    # defaults adapted from DEN Nature paper Fig. 1
    gamma_array = self.params.opt_gamma_values
    weight_array = self.params.opt_weight_values
    grid = []
    for g in gamma_array:
      for w in weight_array:
        grid.append( (g, w) )
    return grid

  def show_den_summary(self, sites_cart):
    print >> self.log, "DEN restraints summary:"
    print >> self.log, "%s | %s | %s | %s | %s " % \
      ("    atom 1     ",
       "    atom 2     ",
       "model dist",
       "  eq dist ",
       "eq dist start")
    for dp in self.den_proxies:
      i_seqs = dp.i_seqs
      a_xyz = sites_cart[i_seqs[0]]
      b_xyz = sites_cart[i_seqs[1]]
      distance_sq = distance_squared(a_xyz, b_xyz)
      distance = distance_sq**(0.5)
      print >> self.log, \
        "%s | %s |   %6.3f   |   %6.3f   |   %6.3f  " % \
        (self.name_hash[i_seqs[0]],
         self.name_hash[i_seqs[1]],
         distance,
         dp.eq_distance,
         dp.eq_distance_start)

  def output_kinemage(self, sites_cart):
    from mmtbx.kinemage import validation
    f = file("den_restraints.kin", "w")
    vec_header = "@kinemage\n"
    vec_header += "@vectorlist {DEN} color= magenta master= {DEN}\n"
    f.write(vec_header)
    for dp in self.den_proxies:
      i_seqs = dp.i_seqs
      eq_distance = dp.eq_distance
      site_a = sites_cart[i_seqs[0]]
      site_b = sites_cart[i_seqs[1]]
      sites = [site_a, site_b]
      #distance_sq = distance_squared(site_a, site_b)
      #distance = distance_sq**(0.5)
      #diff = distance - eq_distance
      #spring = validation.add_spring(sites, diff, "DEN")
      vec = validation.kin_vec(start_key="A",
                               start_xyz=site_a,
                               end_key="B",
                               end_xyz=site_b,
                               width=None)
      f.write(vec)
      #STOP()
    #for chain in self.random_ref_atom_pairs.keys():
    #  for pair in self.random_ref_atom_pairs[chain]:
    #    start_xyz = sites_cart[pair[0]]
    #    end_xyz = sites_cart[pair[1]]
    #    vec = validation.kin_vec(start_key="A",
    #                             start_xyz=start_xyz,
    #                             end_key="B",
    #                             end_xyz=end_xyz,
    #                             width=None)
    #   f.write(vec)
    f.close()
    STOP()

def distance_squared(a, b):
  return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)

