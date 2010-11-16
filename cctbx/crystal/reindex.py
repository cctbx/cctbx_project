from cctbx import sgtbx
from cctbx.sgtbx import cosets
from cctbx.sgtbx import lattice_symmetry


class reindexing_operators(object):
  def __init__(self,
               xs1,
               xs2,
               relative_length_tolerance=0.05,
               absolute_angle_tolerance=10,
               max_delta=3.0,
               anomalous_flag=True,
               out=None):
    # first we have to go to the niggli setting
    self.cb_op_to_n_1 = xs1.change_of_basis_op_to_niggli_cell()
    self.cb_op_to_n_2 = xs2.change_of_basis_op_to_niggli_cell()

    nxs1 = xs1.change_basis( self.cb_op_to_n_1 )
    nxs2 = xs2.change_basis( self.cb_op_to_n_2 )

    # get the lattice symmetry please
    lat_sym_1 = lattice_symmetry.group( nxs1.unit_cell(), max_delta  )
    lat_sym_2 = lattice_symmetry.group( nxs2.unit_cell(), max_delta )
    # and the intensity symmetry please
    int_sym_1 = nxs1.reflection_intensity_symmetry( anomalous_flag ).space_group()
    int_sym_2 = nxs2.reflection_intensity_symmetry( anomalous_flag ).space_group()

    # Now we have to find a similarity transform that maps the two niggli cells onto each other
    c_inv_rs = nxs1.unit_cell().similarity_transformations(
      other=nxs2.unit_cell(),
      relative_length_tolerance=relative_length_tolerance,
      absolute_angle_tolerance=absolute_angle_tolerance)
    min_bases_msd = None
    self.similarity_cb_op = None

    for c_inv_r in c_inv_rs:
      # make the similarity transform into a cb_op
      c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r))
      cb_op = sgtbx.change_of_basis_op(c_inv).inverse()
      # compute the mean square difference for the bases
      bases_msd = nxs1.unit_cell() \
                  .bases_mean_square_difference(
        other=nxs2.unit_cell().change_basis(cb_op=cb_op))

      # and find cb_op correspondiong to the the minimum rmsd
      if (min_bases_msd is None
          or min_bases_msd > bases_msd):
        min_bases_msd = bases_msd
        self.similarity_cb_op = cb_op
    # if nothing is found, do not continue
    self.double_cosets = None
    if (self.similarity_cb_op is not None):
      # make the common lattice group please
      common_lattice_group = lat_sym_1
      for s in lat_sym_2.build_derived_acentric_group().change_basis(
        self.similarity_cb_op ):
        try: common_lattice_group.expand_smx(s)
        except RunTimeError:
          common_lattice_group=None
      if  common_lattice_group is not None:
        common_lattice_group.make_tidy()

        h1 = int_sym_1.build_derived_acentric_group().make_tidy()
        h2 = int_sym_2.build_derived_acentric_group().change_basis( self.similarity_cb_op ).make_tidy()

        # do the double coset decomposition
        self.double_cosets = cosets.double_cosets( common_lattice_group,
                                                   h1,
                                                   h2,
                                                   )
  def cb_ops_in_niggli_setting(self):
    result=[]
    if self.double_cosets is not None:
      for coset in self.double_cosets.double_cosets:
        result.append( sgtbx.change_of_basis_op(coset[0])*self.similarity_cb_op )
    return result

  def combined_cb_ops(self):
    results = []
    if self.double_cosets is not None:
    # now we need to make a combined operator wrt the incomming cell
      for coset in self.double_cosets.double_cosets:
        r =cosets.construct_nice_cb_op(coset,
                                       self.similarity_cb_op,
                                       self.cb_op_to_n_1 ,
                                     self.cb_op_to_n_2 )
        results.append( r )
    #and return it
    return( results )
