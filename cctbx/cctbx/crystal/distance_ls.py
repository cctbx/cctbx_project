from mmtbx import stereochemistry
from cctbx.crystal import minimization
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx.crystal.neighbors import is_sym_equiv_interaction, show_distances
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx import matrix
from libtbx.itertbx import count
from libtbx.test_utils import approx_equal

class create_bond_proxies:

  def __init__(self, structure, distance_cutoff=3.5, distance_ideal=3.1,
                     diff_vec_frac_tolerance=1.e-8):
    distance_cutoff_plus = distance_cutoff * (1+1.e-4)
    distance_cutoff_minus = distance_cutoff * (1-1.e-4)
    asu_mappings = structure.asu_mappings(
      buffer_thickness=distance_cutoff_plus)
    pair_generator = crystal.neighbors_simple_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=distance_cutoff_plus)
    distances_list = [flex.double()
      for i in xrange(structure.scatterers().size())]
    pairs_list = [[]
      for i in xrange(structure.scatterers().size())]
    cb_op_to_niggli_cell = structure.change_of_basis_op_to_niggli_cell()
    for pair in pair_generator:
      diff_vec_frac = structure.unit_cell().fractionalize(pair.diff_vec)
      diff_vec_frac_niggli = matrix.col(cb_op_to_niggli_cell.c()*diff_vec_frac)
      if (diff_vec_frac_niggli.each_abs().max() < 1+diff_vec_frac_tolerance):
        distances_list[pair.i_seq].append(pair.dist_sq)
        pairs_list[pair.i_seq].append(pair)
    for i,distances,pairs in zip(count(),distances_list,pairs_list):
      perm = flex.sort_permutation(distances)
      pairs_list[i] = flex.select(pairs, perm)
    site_symmetries = []
    scatterers = structure.scatterers()
    for scatterer in scatterers:
      site_symmetries.append(structure.site_symmetry(scatterer.site))
    proxies = stereochemistry.restraints_shared_bond_sym_proxy()
    bond_counts = flex.size_t(structure.scatterers().size(), 0)
    for pairs in pairs_list:
      if (len(pairs) > 0):
        i_seq = pairs[0].i_seq
        scatterer_i = scatterers[i_seq]
        site_i = scatterer_i.site
        site_symmetry_i = site_symmetries[i_seq]
        rt_mx_i_inverse = asu_mappings.get_rt_mx(i_seq=i_seq,i_sym=0).inverse()
      for i_pair,pair in zip(count(), pairs):
        assert pair.i_seq == i_seq
        if (pair.dist_sq**.5 > distance_cutoff):
          break
        scatterer_j = scatterers[pair.j_seq]
        site_symmetry_j = site_symmetries[pair.j_seq]
        rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
        rt_mx_ji = rt_mx_i_inverse.multiply(rt_mx_j)
        bond_counts[pair.i_seq] += 1
        if (pair.j_sym == 0):
          bond_counts[pair.j_seq] += 1
        proxies.append(stereochemistry.restraints_bond_sym_proxy(
          pair=pair,
          distance_ideal=distance_ideal,
          weight=100))
        if (pair.dist_sq**.5 > distance_cutoff_minus):
          for pair_fwd in pairs[i_pair+1:]:
            if (    pair_fwd.j_seq == pair.j_seq
                and abs(pair_fwd.dist_sq**.5-pair.dist_sq**.5) < 1.e-4):
              rt_mx_j_fwd = asu_mappings.get_rt_mx(
                i_seq=pair_fwd.j_seq, i_sym=pair_fwd.j_sym)
              rt_mx_ji_fwd = rt_mx_i_inverse.multiply(rt_mx_j_fwd)
              if (is_sym_equiv_interaction(
                    unit_cell=structure.unit_cell(),
                    i_seq=i_seq,
                    site_i=site_i,
                    j_seq=pair.j_seq,
                    site_j=scatterer_j.site,
                    special_op_j=site_symmetries[pair.j_seq].special_op(),
                    rt_mx_ji_1=rt_mx_ji,
                    rt_mx_ji_2=rt_mx_ji_fwd)):
                proxies.append(stereochemistry.restraints_bond_sym_proxy(
                  pair=pair_fwd,
                  distance_ideal=distance_ideal,
                  weight=1))
                bond_counts[i_seq] += 1
                if (pair_fwd.j_sym == 0):
                  bond_counts[pair_fwd.j_seq] += 1
    self.site_symmetries = site_symmetries
    self.asu_mappings = asu_mappings
    self.proxies = proxies
    self.bond_counts = bond_counts

def get_bond_site_symmetry(structure, site_i, site_ji):
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry=structure,
    min_distance_sym_equiv=structure.min_distance_sym_equiv()*.5*(1-1.e-5))
  result = special_position_settings.site_symmetry(
    site=(matrix.col(site_i) + matrix.col(site_ji)) / 2)
  assert result.distance_moved() < 1.e-6
  return result

def add_oxygen(framework_structure, bond_proxies):
  asu_mappings = bond_proxies.asu_mappings
  bonds_processed = [{}
    for i in xrange(framework_structure.scatterers().size())]
  bond_centers = []
  for proxy in bond_proxies.proxies:
    pair = proxy.pair
    if (pair.i_seq > pair.j_seq): continue
    rt_mx_i_inverse=asu_mappings.get_rt_mx(i_seq=pair.i_seq,i_sym=0).inverse()
    rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
    rt_mx_ji = rt_mx_i_inverse.multiply(rt_mx_j)
    ij_rt_mx = bonds_processed[pair.i_seq].setdefault(pair.j_seq, [])
    is_sym_equiv = 00000
    for rt_mx_ji_prev in ij_rt_mx:
      if (is_sym_equiv_interaction(
            unit_cell=framework_structure.unit_cell(),
            i_seq=pair.i_seq,
            site_i=framework_structure.scatterers()[pair.i_seq].site,
            j_seq=pair.j_seq,
            site_j=framework_structure.scatterers()[pair.j_seq].site,
            special_op_j=bond_proxies.site_symmetries[pair.j_seq].special_op(),
            rt_mx_ji_1=rt_mx_ji,
            rt_mx_ji_2=rt_mx_ji_prev)):
        is_sym_equiv = 0001
        break
    if (not is_sym_equiv):
      bond_site_symmetry = get_bond_site_symmetry(
        structure=framework_structure,
        site_i=framework_structure.scatterers()[pair.i_seq].site,
        site_ji=rt_mx_ji*framework_structure.scatterers()[pair.j_seq].site)
      bond_centers.append(bond_site_symmetry.exact_site())
      ij_rt_mx.append(rt_mx_ji)
  complete_structure = framework_structure.deep_copy_scatterers()
  for i,bond_center in zip(count(1), bond_centers):
    complete_structure.add_scatterer(xray.scatterer(
      label="O%d"%i,
      site=bond_center))
  return complete_structure

def run(distance_cutoff=3.5):
  command_line = (iotbx_option_parser(
    usage="python distance_ls.py [options] studat_file [...]",
    description="Example: python distance_ls.py strudat --tag=SOD")
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
  ).process()
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  for file_name in command_line.args:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      if (    command_line.options.tag is not None
          and command_line.options.tag != entry.tag):
        continue
      print "strudat tag:", entry.tag
      structure = entry.as_xray_structure()
      structure.show_summary().show_scatterers()
      show_distances(structure, distance_cutoff=distance_cutoff)
      proxies = create_bond_proxies(
        structure=structure, distance_cutoff=distance_cutoff)
      print "number of bond proxies:", proxies.proxies.size()
      print "proxies.bond_counts:", list(proxies.bond_counts)
      if (proxies.bond_counts.count(4) != proxies.bond_counts.size()):
        print "Not fully 4-connected:", entry.tag
      print "number of bond proxies:", proxies.proxies.size()
      if (0):
        sites_cart = structure.sites_cart()
        gradients_cart = flex.vec3_double(sites_cart.size(), [0,0,0])
        residual_sum = stereochemistry.restraints_bond_residual_sum(
          sites_cart=sites_cart,
          asu_mappings=proxies.asu_mappings,
          proxies=proxies.proxies,
          gradient_array=gradients_cart)
        gradients_frac = gradients_cart \
                       * structure.unit_cell().orthogonalization_matrix()
        print "residual sum:", residual_sum
      if (0):
        for i_site,scatterer in zip(count(), structure.scatterers()):
          site = scatterer.site
          blanks = " "*len(scatterer.label)
          site_symmetry = proxies.site_symmetries[i_site]
          special_op = sgtbx.rt_mx(site_symmetry.special_op().r())
          special_op_tp = float(special_op.r().as_rational()).transpose()
          print scatterer.label, gradients_frac[i_site]
          grad_special=(special_op_tp*matrix.col(gradients_frac[i_site])).elems
          print blanks, grad_special, special_op
          if (not approx_equal(gradients_frac[i_site], grad_special)):
            print "MISMATCH GRAD", structure.space_group_info(), special_op
            raise AssertionError
          elif (0):
            print "    GOOD GRAD", structure.space_group_info(), special_op
      if (1):
        if (0):
          complete_structure = structure
          complete_proxies = proxies
        else:
          complete_structure = add_oxygen(structure, proxies)
          complete_structure.show_summary().show_scatterers()
          show_distances(complete_structure, distance_cutoff=distance_cutoff/2.)
          complete_proxies = create_bond_proxies(
            structure=complete_structure,
            distance_cutoff=distance_cutoff/2.,
            distance_ideal=1.61)
        for proxy in complete_proxies.proxies:
          print "proxy:", proxy.pair.i_seq, proxy.pair.j_seq,
          print proxy.distance_ideal
        if (1):
          sites_cart = complete_structure.sites_cart()
        else:
          sites_frac = flex.vec3_double(flex.random_double(
            size=complete_structure.scatterers().size()*3))
          sites_special = flex.vec3_double()
          for site_frac,site_symmetry in zip(sites_frac,
                                             complete_proxies.site_symmetries):
            sites_special.append(site_symmetry.special_op()*site_frac)
          sites_cart = complete_structure.unit_cell() \
            .orthogonalization_matrix() * sites_special
        minimized = minimization.lbfgs(
          sites_cart=sites_cart,
          site_symmetries=complete_proxies.site_symmetries,
          asu_mappings=complete_proxies.asu_mappings,
          bond_sym_proxies=complete_proxies.proxies,
          lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
            max_iterations=100))
        print minimized.minimizer.error
        print "first:", minimized.first_target_value
        print "final:", minimized.final_target_value
        sites_frac = complete_structure.unit_cell().fractionalization_matrix()\
                   * sites_cart
        minimized_structure = complete_structure.deep_copy_scatterers()
        for scatterer,site in zip(minimized_structure.scatterers(),sites_frac):
          scatterer.site = site
        print "minimized_structure:"
        minimized_structure.show_summary().show_scatterers()
        show_distances(minimized_structure, distance_cutoff=distance_cutoff)
      print

if (__name__ == "__main__"):
  run()
