from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx import data_plots

import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import data_statistics

import libtbx.phil.command_line
from scitbx.python_utils import easy_pickle
import sys, os


class basic_analyses(object):
  def __init__(self,
               miller_array,
               phil_object,
               out=None,
               out_plot=None,
               verbose=0):
    if out is None:
      out=sys.stdout
    if verbose>0:
      print >> out
      print >> out
      print >> out, "Matthews coefficient and Solvent content statistics"

    n_copies_solc = 1.0
    matthews_results =matthews.matthews_rupp(
      miller_array = miller_array,
      n_residues = phil_object.scaling.input.parameters.asu_contents.n_residues,
      n_bases = phil_object.scaling.input.parameters.asu_contents.n_bases,
      out=out,verbose=1)
    phil_object.scaling.input.parameters.asu_contents.n_residues = matthews_results[0]
    phil_object.scaling.input.parameters.asu_contents.n_bases = matthews_results[1]
    n_copies_solc = matthews_results[2]

    if phil_object.scaling.input.parameters.asu_contents.n_copies_per_asu is not None:
      n_copies_solc = phil_object.scaling.input.parameters.asu_contents.n_copies_per_asu
      if verbose>0:
        print >> out,"Number of copies per asyymetric unit provided"
        print >> out," Will use user specified value of ", n_copies_solc
    else:
      phil_object.scaling.input.parameters.asu_contents.n_copies_per_asu = n_copies_solc


    # first report on I over sigma
    miller_array_new = miller_array
    if miller_array.sigmas() is not None:
      data_strength=data_statistics.i_sigi_completeness_stats(
        miller_array,
        isigi_cut = phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.isigi_cut,
        completeness_cut = phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.completeness_cut)
      data_strength.show(out)
      if phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution is None:
        phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution = data_strength.resolution_cut


    ## Isotropic wilson scaling
    if verbose>0:
      print >> out
      print >> out
      print >> out, "Maximum likelihood isotropic Wilson scaling "

    n_residues =  phil_object.scaling.input.parameters.asu_contents.n_residues
    n_bases = phil_object.scaling.input.parameters.asu_contents.n_bases
    if n_residues is None:
      n_residues = 0
    if n_bases is None:
      n_bases = 0
    if n_bases+n_residues==0:
      raise Sorry("No scatterers available")
    iso_scale_and_b = absolute_scaling.ml_iso_absolute_scaling(
      miller_array = miller_array_new,
      n_residues = n_residues*
      miller_array.space_group().order_z()*n_copies_solc,
      n_bases=n_bases*
      miller_array.space_group().order_z()*n_copies_solc)
    iso_scale_and_b.show(out=out,verbose=verbose)
    ## Store the b and scale values from isotropic ML scaling
    self.iso_p_scale = iso_scale_and_b.p_scale
    self.iso_b_wilson =  iso_scale_and_b.b_wilson


    ## Anisotropic ml wilson scaling
    if verbose>0:
      print >> out
      print >> out
      print >> out, "Maximum likelihood anisotropic Wilson scaling "
    aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array = miller_array_new,
      n_residues = n_residues*miller_array.space_group().order_z()*n_copies_solc,
      n_bases = n_bases*miller_array.space_group().order_z()*n_copies_solc)
    aniso_scale_and_b.show(out=out,verbose=1)

    self.aniso_p_scale = aniso_scale_and_b.p_scale
    self.aniso_u_star  = aniso_scale_and_b.u_star
    self.aniso_b_cart  = aniso_scale_and_b.b_cart

    ## Correcting for anisotropicity
    if verbose>0:
      print >> out,"Correcting for anisotropy in the data"
      print >> out

    b_cart_observed = aniso_scale_and_b.b_cart

    b_trace_average = (b_cart_observed[0]+
                       b_cart_observed[1]+
                       b_cart_observed[2])/3.0
    b_trace_min = b_cart_observed[0]
    if  b_cart_observed[1] <b_trace_min: b_trace_min=b_cart_observed[1]
    if  b_cart_observed[2] <b_trace_min: b_trace_min=b_cart_observed[2]

    if phil_object.scaling.input.optional.aniso.final_b == "eigen_min":
       b_use=aniso_scale_and_b.eigen_values[2]
    elif phil_object.scaling.input.optional.aniso.final_b == "eigen_mean" :
       b_use=flex.mean(aniso_scale_and_b.eigen_values)
    elif phil_object.scaling.input.optional.aniso.final_b == "use_b_iso":
       assert phil_object.scaling.input.optional.b_iso is not None
       b_use=phil_object.scaling.input.optional.b_iso
    else:
       b_use = 30

    b_cart_aniso_removed = [ -b_use,
                             -b_use,
                             -b_use,
                             0,
                             0,
                             0]
    u_star_aniso_removed = adptbx.u_cart_as_u_star(
      miller_array.unit_cell(),
      adptbx.b_as_u( b_cart_aniso_removed  ) )
    ## I do things in two steps, but can easely be done in 1 step
    ## just for clarity, thats all.
    self.no_aniso_array = absolute_scaling.anisotropic_correction(
      miller_array_new,0.0,aniso_scale_and_b.u_star )
    self.no_aniso_array = absolute_scaling.anisotropic_correction(
      self.no_aniso_array,0.0,u_star_aniso_removed)
    self.no_aniso_array = self.no_aniso_array.set_observation_type(
      miller_array )

    """
    ## write out this miller array as sca if directed to do so:
    output_file=phil_object.scaling.input.optional.hklout
    if output_file is not None and \
     phil_object.scaling.input.optional.aniso.action=="remove_aniso":
      if phil_object.scaling.input.optional.hklout_type == "sca":
        import iotbx.scalepack.merge
        iotbx.scalepack.merge.write(
            file_name=output_file,miller_array=self.no_aniso_array)
      elif phil_object.scaling.input.optional.hklout_type == "mtz":
        mtz_dataset = self.no_aniso_array.as_mtz_dataset(
            column_root_label='F_aniso')
        mtz_dataset.mtz_object().write(output_file)
      print>>out, "Wrote aniso-corrected data to ",output_file
    """

    ## Make normalised structure factors please
    normalistion = absolute_scaling.kernel_normalisation(
      self.no_aniso_array,auto_kernel=True)
    self.normalised_miller = normalistion.normalised_miller.deep_copy()


    self.phil_object=phil_object

    ## Some basic statistics and sanity checks follow
    if verbose>0:
      print >> out,"Some basic intensity statistics follow."
      print >> out

    basic_data_stats = data_statistics.basic_intensity_statistics(
      miller_array,
      aniso_scale_and_b.p_scale,
      aniso_scale_and_b.u_star,
      iso_scale_and_b.scat_info,
      out=out,
      out_plot=out_plot)

    self.miller_array = basic_data_stats.new_miller


    if verbose>0:
      print >> out, "Basic analyses completed"
