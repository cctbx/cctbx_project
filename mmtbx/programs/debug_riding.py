from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
from mmtbx.hydrogens import connectivity
import mmtbx.model
import iotbx.pdb


master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
iseq = None
  .type = int(value_min=0)
'''

# =============================================================================

class Program(ProgramTemplate):
  description = '''

Print parameterization

'''

  datatypes = ['model', 'restraint', 'phil']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

  # ---------------------------------------------------------------------------

  def run(self):

    model_fn = self.data_manager.get_default_model_name()
    print('\nUsing model:', model_fn, file=self.logger)

    model = self.data_manager.get_model()
    model.set_log(log = null_out())
    model.set_stop_for_unknowns(False)
    model.process(make_restraints=True)

    #model.setup_riding_h_manager()

    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    geometry = model.restraints_manager.geometry
    sites_cart = model.get_sites_cart()

    connectivity_manager = connectivity.determine_connectivity(
      pdb_hierarchy       = hierarchy,
      geometry_restraints = geometry)
    h_connectivity = connectivity_manager.h_connectivity

    self.print_connectivity(
      h_connectivity = h_connectivity,
      atoms = atoms)

    bond_proxies_simple, asu = geometry.get_all_bond_proxies(sites_cart=sites_cart)
    angle_proxies = geometry.get_all_angle_proxies()
    dihedral_proxies = geometry.dihedral_proxies
    planarity_proxies = geometry.planarity_proxies

    make_sub_header('Proxies', out=self.logger)

    if self.params.iseq:
      print('Showing proxies for particular atom:', file=self.logger)
      print(self.clean(atoms[self.params.iseq].id_str()), file=self.logger)

    self.print_bond_proxies(
      bond_proxies = bond_proxies_simple,
      atoms = atoms,
      iseq = self.params.iseq)
    self.print_angle_proxies(
      angle_proxies = angle_proxies,
      atoms = atoms,
      iseq = self.params.iseq)
    self.print_dihedral_proxies(
      dihedral_proxies = dihedral_proxies,
      atoms = atoms,
      iseq = self.params.iseq)
    self.print_planarity_proxies(
      planarity_proxies = planarity_proxies,
      atoms = atoms,
      iseq = self.params.iseq)

  # ----------------------------------------------------------------------------

  def clean(self, s):
    return s.replace('pdb=', '').replace('"', '')#.strip()

  def print_bond_proxies(self, bond_proxies, atoms, iseq = None):
    print('\nbond proxies', file=self.logger)
    for bp in bond_proxies:
      if iseq and iseq not in bp.i_seqs: continue
      i1, i2 = bp.i_seqs
      s1 = self.clean(atoms[i1].id_str())
      s2 = self.clean(atoms[i2].id_str())
      print("\t%s --- %s " %  (s1, s2), file=self.logger)

  def print_angle_proxies(self, angle_proxies, atoms, iseq = None):
    print('\nangle proxies', file=self.logger)
    for ap in angle_proxies:
      if iseq  and iseq not in ap.i_seqs: continue
      i1, i2, i3 = ap.i_seqs
      s1 = self.clean(atoms[i1].id_str())
      s2 = self.clean(atoms[i2].id_str())
      s3 = self.clean(atoms[i3].id_str())
      print('\t%s - %s - %s  (angle ideal: %s)' % (s1, s2, s3, ap.angle_ideal))

  def print_dihedral_proxies(self, dihedral_proxies, atoms, iseq = None):
    print('\ndihedral proxies', file=self.logger)
    for ap in dihedral_proxies:
      if iseq and iseq not in ap.i_seqs: continue
      i1, i2, i3, i4 = ap.i_seqs
      s1 = self.clean(atoms[i1].id_str())
      s2 = self.clean(atoms[i2].id_str())
      s3 = self.clean(atoms[i3].id_str())
      s4 = self.clean(atoms[i4].id_str())
      print('\t%s - %s - %s - %s  (angle ideal: %s)' % (s1, s2, s3, s4,
        ap.angle_ideal), file=self.logger)

  def print_planarity_proxies(self, planarity_proxies, atoms, iseq = None):
    print('\nplanarity proxies', file=self.logger)
    for pp in planarity_proxies:
      if iseq and iseq not in pp.i_seqs: continue
      #print(*pp.i_seqs, sep = ", ")
      for iseq in pp.i_seqs:
        print('\t', self.clean(atoms[iseq].id_str()), file=self.logger)

  # ----------------------------------------------------------------------------

  def print_connectivity(self, h_connectivity, atoms):
    """Print information of connectivity for each H atom."""
    make_sub_header(' H-atom connectivity ', out=self.logger)
    names = list(atoms.extract_name())
    for neighbors in h_connectivity:
      if neighbors is None: continue
      ih = neighbors.ih
      #if (ih != 1620): continue
      labels = atoms[ih].fetch_labels()
      #if (labels.resseq.strip() != '2'): continue
      #print labels.resseq.strip()
      if (neighbors.number_non_h_neighbors == 0):
        print('%s %s (%s)' % (names[ih], ih, labels.resseq.strip()))
      else:
        i_a0 = neighbors.a0['iseq']
        i_a1 = neighbors.a1['iseq']
        string = names[i_a1]+'('+str(i_a1)+')'
        if 'iseq' in neighbors.a2:
          i_a2 = neighbors.a2['iseq']
          string = string + names[i_a2]
        if 'iseq' in neighbors.a3:
          string = string + names[neighbors.a3['iseq']]
        output = (names[ih], ih, labels.resseq.strip(), names[i_a0], i_a0, string)
        stringh = ''
        if 'iseq' in neighbors.h1:
          stringh = stringh + names[neighbors.h1['iseq']]
        if 'iseq' in neighbors.h2:
          stringh = stringh + names[neighbors.h2['iseq']]
        if 'iseq' in neighbors.b1:
          stringb1 = names[neighbors.b1['iseq']]
        else:
          stringb1 = 'n/a'
        output = output + (stringh,) + (stringb1,)#+ (self.third_neighbors_raw[i_a0],)
        print('%s %s (%s) , %s (%s) , %s , %s,%s' % output)
