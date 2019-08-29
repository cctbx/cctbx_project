from __future__ import absolute_import, division, print_function
from cctbx import adptbx
from cctbx.array_family import flex
from libtbx.utils import Sorry, show_exception_info_if_full_testing
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, relative_wilson
from mmtbx.scaling import matthews
from mmtbx.scaling import data_statistics
import iotbx.bioinformatics
import sys


class basic_analyses(object):  # XXX is this ever used?
  def __init__(self,
               miller_array,
               phil_object,
               out=None,
               out_plot=None, miller_calc=None,
               original_intensities=None,
               completeness_as_non_anomalous=None,
               verbose=0):
    if out is None:
      out=sys.stdout
    if verbose>0:
      print(file=out)
      print(file=out)
      print("Matthews coefficient and Solvent content statistics", file=out)
    n_copies_solc = 1.0
    self.nres_known = False
    if (phil_object.scaling.input.asu_contents.n_residues is not None or
        phil_object.scaling.input.asu_contents.n_bases is not None):
      self.nres_known = True
      if (phil_object.scaling.input.asu_contents.sequence_file is not None):
        print("  warning: ignoring sequence file", file=out)
    elif (phil_object.scaling.input.asu_contents.sequence_file is not None):
      print("  determining composition from sequence file %s" % \
        phil_object.scaling.input.asu_contents.sequence_file, file=out)
      seq_comp = iotbx.bioinformatics.composition_from_sequence_file(
        file_name=phil_object.scaling.input.asu_contents.sequence_file,
        log=out)
      if (seq_comp is not None):
        phil_object.scaling.input.asu_contents.n_residues = seq_comp.n_residues
        phil_object.scaling.input.asu_contents.n_bases = seq_comp.n_bases
        self.nres_known = True
    matthews_results =matthews.matthews_rupp(
      crystal_symmetry = miller_array,
      n_residues = phil_object.scaling.input.asu_contents.n_residues,
      n_bases = phil_object.scaling.input.asu_contents.n_bases,
      out=out,verbose=1)
    phil_object.scaling.input.asu_contents.n_residues = matthews_results[0]
    phil_object.scaling.input.asu_contents.n_bases = matthews_results[1]
    n_copies_solc = matthews_results[2]
    self.matthews_results = matthews_results

    if phil_object.scaling.input.asu_contents.n_copies_per_asu is not None:
      n_copies_solc = phil_object.scaling.input.asu_contents.n_copies_per_asu
      self.defined_copies = n_copies_solc
      if verbose>0:
        print("Number of copies per asymmetric unit provided", file=out)
        print(" Will use user specified value of ", n_copies_solc, file=out)
    else:
      phil_object.scaling.input.asu_contents.n_copies_per_asu = n_copies_solc
      self.guessed_copies = n_copies_solc

    # first report on I over sigma
    miller_array_new = miller_array
    self.data_strength = None
    miller_array_intensities = miller_array
    if (original_intensities is not None):
      assert original_intensities.is_xray_intensity_array()
      miller_array_intensities = original_intensities
    if miller_array_intensities.sigmas() is not None:
      data_strength=data_statistics.i_sigi_completeness_stats(
        miller_array_intensities,
        isigi_cut = phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.isigi_cut,
        completeness_cut = phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.completeness_cut,
      completeness_as_non_anomalous=completeness_as_non_anomalous)
      data_strength.show(out)
      self.data_strength = data_strength
      if phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution is None:
        if data_strength.resolution_cut > data_strength.resolution_at_least:
          phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution = data_strength.resolution_at_least
        else:
           phil_object.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution = data_strength.resolution_cut

    ## Isotropic wilson scaling
    if verbose>0:
      print(file=out)
      print(file=out)
      print("Maximum likelihood isotropic Wilson scaling ", file=out)

    n_residues =  phil_object.scaling.input.asu_contents.n_residues
    n_bases = phil_object.scaling.input.asu_contents.n_bases
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
    self.iso_scale_and_b = iso_scale_and_b
    ## Store the b and scale values from isotropic ML scaling
    self.iso_p_scale = iso_scale_and_b.p_scale
    self.iso_b_wilson =  iso_scale_and_b.b_wilson


    ## Anisotropic ml wilson scaling
    if verbose>0:
      print(file=out)
      print(file=out)
      print("Maximum likelihood anisotropic Wilson scaling ", file=out)
    aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array = miller_array_new,
      n_residues = n_residues*miller_array.space_group().order_z()*n_copies_solc,
      n_bases = n_bases*miller_array.space_group().order_z()*n_copies_solc)
    aniso_scale_and_b.show(out=out,verbose=1)

    self.aniso_scale_and_b = aniso_scale_and_b

    try: b_cart = aniso_scale_and_b.b_cart
    except AttributeError as e:
      print("*** ERROR ***", file=out)
      print(str(e), file=out)
      show_exception_info_if_full_testing()
      return

    self.aniso_p_scale = aniso_scale_and_b.p_scale
    self.aniso_u_star  = aniso_scale_and_b.u_star
    self.aniso_b_cart  = aniso_scale_and_b.b_cart
    # XXX: for GUI
    self.overall_b_cart = getattr(aniso_scale_and_b, "overall_b_cart", None)

    ## Correcting for anisotropy
    if verbose>0:
      print("Correcting for anisotropy in the data", file=out)
      print(file=out)

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
    elif phil_object.scaling.input.optional.aniso.final_b == "user_b_iso":
       assert phil_object.scaling.input.optional.aniso.b_iso is not None
       b_use=phil_object.scaling.input.optional.aniso.b_iso
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


    ## Make normalised structure factors please

    sel_big = self.no_aniso_array.data() > 1.e+50
    self.no_aniso_array = self.no_aniso_array.array(
      data = self.no_aniso_array.data().set_selected(sel_big, 0))
    self.no_aniso_array = self.no_aniso_array.set_observation_type(
      miller_array )

    normalistion = absolute_scaling.kernel_normalisation(
      self.no_aniso_array,auto_kernel=True)
    self.normalised_miller = normalistion.normalised_miller.deep_copy()


    self.phil_object=phil_object

    ## Some basic statistics and sanity checks follow
    if verbose>0:
      print("Some basic intensity statistics follow.", file=out)
      print(file=out)

    basic_data_stats = data_statistics.basic_intensity_statistics(
      miller_array,
      aniso_scale_and_b.p_scale,
      aniso_scale_and_b.u_star,
      iso_scale_and_b.scat_info,
      out=out,
      out_plot=out_plot)
    self.basic_data_stats = basic_data_stats
    self.miller_array = basic_data_stats.new_miller

    #relative wilson plot
    self.rel_wilson = None
    if (miller_calc is not None) and (miller_calc.d_min() < 4.0):
      try :
        self.rel_wilson = relative_wilson.relative_wilson(
          miller_obs=miller_array,
          miller_calc=miller_calc)
      except RuntimeError as e :
        print("*** Error calculating relative Wilson plot - skipping.", file=out)
        print("", file=out)

    if verbose>0:
      print("Basic analyses completed", file=out)
