from mmtbx import stereochemistry
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx.crystal.neighbors import is_sym_equiv_interaction, show_distances
from cctbx.array_family import flex
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
          weight=1))
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
      orth_mx = matrix.sqr(structure.unit_cell().orthogonalization_matrix())
      if (1):
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
      print

if (__name__ == "__main__"):
  run()
