
from __future__ import absolute_import, division, print_function
from mmtbx.scaling import outlier_rejection
from mmtbx.scaling import absolute_scaling
import mmtbx.scaling
import iotbx.phil
from cctbx.array_family import flex
from cctbx import miller
from cctbx import adptbx
from libtbx.utils import Sorry, null_out
from libtbx import Auto
import os.path
import sys

output_params_str = """
  hklout = None
    .type = path
  hklout_type=mtz sca *Auto
    .type = choice
  label_extension="massaged"
    .type = str
"""

master_params = iotbx.phil.parse("""
    aniso
      .help="Parameters dealing with anisotropy correction"
      .short_caption = Anisotropy correction
      .style = box auto_align
    {
      action=*remove_aniso None
        .type=choice
        .caption = Remove_anisotropy None
        .help="Remove anisotropy?"
        .style = bold
      final_b=*eigen_min eigen_mean user_b_iso
        .help="Final b value"
        .type=choice
        .caption = Minimum_eigenvalue Mean_eigenvalue User_specified
        .short_caption = Final B-factor source
      b_iso=None
        .type=float
        .help="User specified B value"
        .short_caption = User-specified B-factor
    }
    outlier
      .help="Outlier analyses"
      .short_caption = Outlier analyses
      .style = box auto_align
    {
      action=*extreme basic beamstop None
        .help="Outlier protocol"
        .type=choice
        .short_caption = Outlier rejection protocol
      parameters
        .help="Parameters for outlier detection"
      {
        basic_wilson{
          level = 1E-6
            .type=float
            .short_caption = Basic Wilson level
        }
        extreme_wilson{
          level = 0.01
            .type=float
            .short_caption = Extreme Wilson level
        }
        beamstop{
          level = 0.001
            .type=float
            .short_caption = Beamstop level
          d_min = 10.0
            .type=float
            .short_caption = Max. resolution
            .style = resolution
        }
      }
    }
    symmetry
      .short_caption = Detwinning
      .style = box auto_align
    {
      action = detwin twin *None
        .type = choice
        .short_caption = Action
        .style = bold
      twinning_parameters{
        twin_law = None
          .type = str
          .short_caption = Twin law
          .input_size = 120
          .style = bold
        fraction = None
          .type = float
          .short_caption = Detwinning fraction
          .style = bold
      }
    }
""")

class massage_data(object):
  def __init__(self,
               miller_array,
               parameters,
               out=None,
               n_residues=100,
               n_bases=0):

    self.params=parameters
    self.miller_array=miller_array.deep_copy().set_observation_type(
      miller_array).merge_equivalents().array()
    self.out = out
    if self.out is None:
      self.out = sys.stdout
    if self.out == "silent":
      self.out = null_out()


    self.no_aniso_array = self.miller_array
    if self.params.aniso.action == "remove_aniso":
      # first perfom aniso scaling
      aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(
        miller_array = self.miller_array,
        n_residues = n_residues,
        n_bases = n_bases)
      aniso_scale_and_b.p_scale = 0 # set the p_scale back to 0!
      aniso_scale_and_b.show(out=out)
      # now do aniso correction please
      self.aniso_p_scale = aniso_scale_and_b.p_scale
      self.aniso_u_star  = aniso_scale_and_b.u_star
      self.aniso_b_cart  = aniso_scale_and_b.b_cart
      if self.params.aniso.final_b == "eigen_min":
        b_use=aniso_scale_and_b.eigen_values[2]
      elif self.params.aniso.final_b == "eigen_mean" :
        b_use=flex.mean(aniso_scale_and_b.eigen_values)
      elif self.params.aniso.final_b == "user_b_iso":
        assert self.params.aniso.b_iso is not None
        b_use=self.params.aniso.b_iso
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
        self.miller_array,0.0,aniso_scale_and_b.u_star )
      self.no_aniso_array = absolute_scaling.anisotropic_correction(
        self.no_aniso_array,0.0,u_star_aniso_removed)
      self.no_aniso_array = self.no_aniso_array.set_observation_type(
        miller_array )

    # that is done now, now we can do outlier detection if desired
    outlier_manager = outlier_rejection.outlier_manager(
      self.no_aniso_array,
      None,
      out=self.out)


    self.new_miller_array = self.no_aniso_array
    if self.params.outlier.action == "basic":
      print("Non-outliers found by the basic wilson statistics", file=self.out)
      print("protocol will be written out.", file=self.out)
      basic_array = outlier_manager.basic_wilson_outliers(
        p_basic_wilson = self.params.outlier.parameters.basic_wilson.level,
        return_data = True)
      self.new_miller_array = basic_array

    if self.params.outlier.action == "extreme":
      print("Non-outliers found by the extreme value wilson statistics", file=self.out)
      print("protocol will be written out.", file=self.out)
      extreme_array = outlier_manager.extreme_wilson_outliers(
      p_extreme_wilson = self.params.outlier.parameters.extreme_wilson.level,
      return_data = True)
      self.new_miller_array = extreme_array

    if self.params.outlier.action == "beamstop":
      print("Outliers found for the beamstop shadow", file=self.out)
      print("problems detection protocol will be written out.", file=self.out)
      beamstop_array = outlier_manager.beamstop_shadow_outliers(
        level = self.params.outlier.parameters.beamstop.level,
        d_min = self.params.outlier.parameters.beamstop.d_min,
        return_data=True)
      self.new_miller_array = beamstop_array

    if self.params.outlier.action == "None":
      self.new_miller_array =  self.no_aniso_array



    # now we can twin or detwin the data if needed
    self.final_array = self.new_miller_array
    if self.params.symmetry.action == "twin":
      alpha = self.params.symmetry.twinning_parameters.fraction
      if (alpha is None):
        raise Sorry("Twin fraction not specified, not twinning data")
      elif not (0 <= alpha <= 0.5):
        raise Sorry("Twin fraction must be between 0 and 0.5.")
      print(file=self.out)
      print("Twinning given data", file=self.out)
      print("-------------------", file=self.out)
      print(file=self.out)
      print("Artifically twinning the data with fraction %3.2f" %\
        alpha, file=self.out)

      self.final_array = self.new_miller_array.twin_data(
        twin_law = self.params.symmetry.twinning_parameters.twin_law,
        alpha=alpha).as_intensity_array()

    elif (self.params.symmetry.action == "detwin"):
      twin_law = self.params.symmetry.twinning_parameters.twin_law
      alpha = self.params.symmetry.twinning_parameters.fraction
      if (alpha is None):
        raise Sorry("Twin fraction not specified, not detwinning data")
      elif not (0 <= alpha <= 0.5):
        raise Sorry("Twin fraction must be between 0 and 0.5.")
      print("""

Attempting to detwin data
-------------------------
Detwinning data with:
  - twin law:      %s
  - twin fraciton: %.2f

BE WARNED! DETWINNING OF DATA DOES NOT SOLVE YOUR TWINNING PROBLEM!
PREFERABLY, REFINEMENT SHOULD BE CARRIED OUT AGAINST ORIGINAL DATA
ONLY USING A TWIN SPECIFIC TARGET FUNCTION!

""" % (twin_law, alpha), file=self.out)
      self.final_array = self.new_miller_array.detwin_data(
        twin_law=twin_law,
        alpha=alpha).as_intensity_array()

    assert self.final_array is not None

  def return_data(self):
    return self.final_array

  def write_data(self,
      file_name,
      output_type=Auto,
      label_extension="massaged"):
    ## write out this miller array as sca if directed to do so:
    if (str(output_type) == "Auto"):
      base, ext = os.path.splitext(file_name)
      if ext in [".mtz",".sca"]:
        output_type = ext[1:]
      else:
        raise Sorry("Unknown or unsupported output type")
    assert (output_type in ["mtz", "sca"]), output_type
    if output_type == "sca":
      import iotbx.scalepack.merge
      iotbx.scalepack.merge.write(
        file_name=file_name,
        miller_array=self.final_array,
        scale_intensities_for_scalepack_merge=True) # scales only if necessary
    elif output_type == "mtz":
      base_label=None
      if self.final_array.is_xray_intensity_array():
        base_label = "I"
      if self.final_array.is_xray_amplitude_array():
        base_label = "F"
      mtz_dataset = self.final_array.as_mtz_dataset(
        column_root_label=base_label+label_extension)
      mtz_dataset.mtz_object().write(file_name)
