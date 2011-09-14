class potential_object(object):

  def __init__(O, f_obs, xray_structure):
    O.f_obs = f_obs
    O.xray_structure = xray_structure
    O.last_sites_moved = None
    O.f = None

  def e_pot(O, sites_moved):
    if (   O.last_sites_moved is None
        or O.last_sites_moved.id() != sites_moved.id()):
      O.last_sites_moved = sites_moved
      xs = O.xray_structure
      assert len(sites_moved) == xs.scatterers().size()
      xs.set_sites_cart(sites_cart=sites_moved)
      f_calc = O.f_obs.structure_factors_from_scatterers(
        xray_structure=xs, algorithm="direct", cos_sin_table=False).f_calc()
      from cctbx import xray
      tg = xray.targets_least_squares(
        compute_scale_using_all_data=True,
        obs_type="F",
        obs=O.f_obs.data(),
        weights=None,
        r_free_flags=None,
        f_calc=f_calc.data(),
        derivatives_depth=0,
        scale_factor=0)
      O.f = tg.target_work()
    return O.f

def sample_e_pot(id_code, f_obs, xray_structure, edge_list, params):
  if (edge_list is None):
    print "NO_TARDY: no edge_list"
    return
  #
  xs = xray_structure
  if (xs.special_position_indices().size() != 0):
    print "NO_TARDY: special positions"
    return
  sites_cart = xs.sites_cart()
  labels = xs.scatterers().extract_labels()
  pat = xs.pair_asu_table(distance_cutoff=2.5)
  if (0):
    pat.show_distances(sites_cart=sites_cart, site_labels=labels)
  pst = pat.extract_pair_sym_table()
  for i,j in edge_list:
    assert i <= j
    sym_dict = pst[i].get(j)
    if (sym_dict is None):
      print "NO_TARDY: large distance for edge:", labels[i], labels[j]
      return
  #
  import scitbx.graph.tardy_tree
  tt = scitbx.graph.tardy_tree.construct(
    sites=sites_cart, edge_list=edge_list)
  import scitbx.rigid_body
  tardy_model = scitbx.rigid_body.tardy_model(
    labels=labels,
    sites=sites_cart,
    masses=xs.atomic_weights(),
    tardy_tree=tt,
    potential_obj=potential_object(f_obs, xs),
    near_singular_hinges_angular_tolerance_deg=5)
  if (tardy_model.number_of_trees != 1):
    print "NO_TARDY: multiple trees"
    return
  #
  print "Single tardy tree:", \
    id_code, xs.scatterers().size(), xs.space_group_info()
  tt.show_summary(vertex_labels=labels, prefix="  ")
  print "dof each joint:", list(tardy_model.degrees_of_freedom_each_joint())
  print "q_size each joint:", list(tardy_model.q_size_each_joint())
  q_packed = tardy_model.pack_q()
  print "q_packed.size():", q_packed.size()
  print "q_packed:", numstr(q_packed)
  if (params.iq < 0):
    return
  assert params.iq < q_packed.size()
  #
  xy = []
  from libtbx.utils import xsamples
  from math import pi
  for q_deg in xsamples(params.qmin, params.qmax, params.qstep):
    q_packed[params.iq] = q_deg*pi/180
    tardy_model.unpack_q(q_packed=q_packed)
    e_pot = tardy_model.e_pot()
    xy.append((q_deg,e_pot))
  from libtbx import pyplot
  pyplot.plot_pairs(xy, "r-")
  pyplot.show()
