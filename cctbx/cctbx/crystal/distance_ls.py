from mmtbx import stereochemistry
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx.crystal.neighbors import is_sym_equiv_interaction, show_distances
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
from scitbx import matrix
from scitbx.python_utils.math_utils import iround
from libtbx.itertbx import count
from libtbx.test_utils import approx_equal
import math
import sys

def get_bond_site_symmetry(structure, site_i, site_ji):
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry=structure,
    min_distance_sym_equiv=structure.min_distance_sym_equiv()*.5*(1-1.e-5))
  result = special_position_settings.site_symmetry(
    site=(matrix.col(site_i) + matrix.col(site_ji)) / 2)
  assert result.distance_moved() < 1.e-6
  return result

class create_bond_proxies:

  def __init__(self, structure, distance_cutoff=3.5, distance_ideal=3.1,
                     diff_vec_frac_tolerance=1.e-8):
    distance_cutoff_plus = distance_cutoff * (1+1.e-4)
    distance_cutoff_minus = distance_cutoff * (1-1.e-4)
    asu_mappings = structure.asu_mappings(
      buffer_thickness=distance_cutoff_plus)
    pair_generator = crystal.neighbors_fast_pair_generator(
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
    bond_counts = flex.double(structure.scatterers().size(), 0)
    for pairs in pairs_list:
      bond_centers_dict = {}
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
        bond_centers = bond_centers_dict.get(pair.j_seq, None)
        if (bond_centers is None):
          bond_centers = bond_centers_dict[pair.j_seq] = flex.vec3_double()
        scatterer_j = scatterers[pair.j_seq]
        rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
        rt_mx_ji = rt_mx_i_inverse.multiply(rt_mx_j)
        if (pair.j_seq >= i_seq):
          bond_site_symmetry = get_bond_site_symmetry(
            structure=structure,
            site_i=site_i,
            site_ji=rt_mx_ji*scatterer_j.site)
          if (   bond_centers.size() == 0
              or sgtbx.min_sym_equiv_distance_info(
                   reference_sites=sgtbx.sym_equiv_sites(bond_site_symmetry),
                   others=bond_centers).dist() > 1.e-6):
            bond_centers.append(bond_site_symmetry.exact_site())
            proxies.append(stereochemistry.restraints_bond_sym_proxy(
              pair=pair,
              distance_ideal=distance_ideal,
              weight=bond_site_symmetry.multiplicity()))
        bond_counts[i_seq] += 1
        if (pair.j_seq != i_seq):
          n_i = len(site_symmetries[i_seq].matrices())
          n_j = len(site_symmetries[pair.j_seq].matrices())
          f_j = float(n_j) / n_i
          bond_counts[pair.j_seq] += f_j
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
                bond_counts[i_seq] += 1
                if (pair.j_seq != i_seq):
                  bond_counts[pair.j_seq] += f_j
    self.site_symmetries = site_symmetries
    self.asu_mappings = asu_mappings
    self.proxies = proxies
    self.bond_counts = flex.size_t([iround(c) for c in bond_counts])
    assert flex.max(flex.abs(bond_counts-self.bond_counts.as_double())) < 1.e-8

def run(distance_cutoff=3.5):
  command_line = (iotbx_option_parser(
    usage="python distance_ls.py [options] studat_file [...]",
    description="Example: python distance_ls.py strudat --tag=SOD")
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag is it appears in the strudat file")
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
      sites_cart = structure.sites_cart()
      gradient_array = flex.vec3_double(sites_cart.size(), [0,0,0])
      residual_sum = stereochemistry.restraints_bond_residual_sum(
        sites_cart=sites_cart,
        asu_mappings=proxies.asu_mappings,
        proxies=proxies.proxies,
        gradient_array=gradient_array)
      print "residual sum:", residual_sum
      for grad, site_symmetry in zip(gradient_array, proxies.site_symmetries):
        print grad, sgtbx.rt_mx(site_symmetry.special_op().r())
      structure_p1 = structure.expand_to_p1(append_number_to_labels=0001)
      structure_p1.show_summary().show_scatterers()
      show_distances(structure_p1, distance_cutoff=distance_cutoff)
      proxies_p1 = create_bond_proxies(
        structure=structure_p1, distance_cutoff=distance_cutoff)
      print "number of bond proxies_p1:", proxies_p1.proxies.size()
      print "proxies_p1.bond_counts:", list(proxies_p1.bond_counts)
      residual_sum_p1 = stereochemistry.restraints_bond_residual_sum(
        sites_cart=structure_p1.sites_cart(),
        asu_mappings=proxies_p1.asu_mappings,
        proxies=proxies_p1.proxies,
        gradient_array=None)
      print "residual sum p1:", residual_sum_p1
      if (residual_sum != 0):
        print "residual ratio: %.6g" % (residual_sum_p1 / residual_sum),
        print structure.space_group_info()
      assert abs(residual_sum_p1-residual_sum) < 1.e-6
      print

if (__name__ == "__main__"):
  run()
