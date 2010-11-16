"""
The following phil scope is expected for the massage_data class:
     hklout = None
     hklout_type=mtz sca *mtz_or_sca
     label_extension="massaged"
     aniso{
       action=*remove_aniso None
       final_b=*eigen_min eigen_mean user_b_iso
       b_iso=None
     }
     outlier{
       action=*extreme basic beamstop None
       parameters{
         basic_wilson{
          level=1E-6
         }
         extreme_wilson{
           level=0.01
         }
         beamstop{
           level=0.001
           d_min=10.0
         }
       }
     }
     symmetry{
       action=detwin twin *None
       twinning_parameters{
         twin_law=None
         fraction=None
       }
     }
   }

"""




from cctbx import miller
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import outlier_rejection
from libtbx.utils import null_out
import sys


master_params = iotbx.phil.parse("""
    hklout = None
      .type=path
      .help="HKL out"
      .short_caption = Output reflections file
      .style = bold new_file
    hklout_type=mtz sca *mtz_or_sca
      .type=choice
      .help="Output format"
      .caption = MTZ Scalepack MTZ_or_Scalepack
      .short_caption = Output format
    label_extension="massaged"
      .type=str
      .help="Label extension"
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



class twin_data(object):
  def __init__(self,
               miller_array,
               twin_law,
               out=None):
    self.out=out
    if self.out is None:
      self.out=sys.stdout

    print >> self.out
    print >> self.out, "Twinning given data"
    print >> self.out, "-------------------"
    print >> self.out
    self.miller_array = miller_array.deep_copy().set_observation_type(
      miller_array ).map_to_asu()

    self.twin_law = twin_law
    assert (self.twin_law is not None)
    self.twin_law=sgtbx.rt_mx(self.twin_law, r_den=24,t_den=288 )

    if self.twin_law.r().determinant() != 1:
      raise Sorry("The determinant of the provided twin law is not equal to unity")

  def twin_it(self,alpha):
    print >> self.out, "Artifically twinning the data with fraction %3.2f"%(alpha)

    assert alpha is not None
    assert alpha<=0.5
    assert alpha>=0.0
    # make sure we have intensities
    if self.miller_array.is_real_array():
      if not self.miller_array.is_xray_intensity_array():
        self.miller_array = self.miller_array.f_as_f_sq()
    assert self.miller_array.is_xray_intensity_array()

    cb_op = sgtbx.change_of_basis_op( self.twin_law )
    print >> self.out, "using twin law (%s)"%( cb_op.as_hkl() )

    self.new_miller = self.miller_array.change_basis( cb_op ).map_to_asu()
    xa,xb = self.miller_array.common_sets( self.new_miller )
    new_data = (1.0-alpha)*xa.data() + alpha*xb.data()
    xa = xa.customized_copy(data=new_data,
                            sigmas=new_data/100.0).set_observation_type( self.miller_array )
    return xa


class detwin_data(object):
  def __init__(self,
               miller_array,
               twin_law,
               out=None
               ):
    self.out=out
    if self.out is None:
      self.out=sys.stdout
    print >> self.out
    print >> self.out
    print >> self.out, "Attempting to detwin data"
    print >> self.out, "-------------------------"
    print >> self.out, "Detwinning data with:"
    print >> self.out, "  - twin law:      %s"%(twin_law)
    print >> self.out
    print >> self.out, "BE WARNED! DETWINNING OF DATA DOES NOT SOLVE YOUR TWINNING PROBLEM!"
    print >> self.out, "PREFERABLY, REFINEMENT SHOULD BE CARRIED OUT AGAINST ORIGINAL DATA "
    print >> self.out, "ONLY USING A TWIN SPECIFIC TARGET FUNCTION!"
    print >> self.out

    self.miller_array = miller_array.deep_copy().set_observation_type(
      miller_array )

    self.twin_law = twin_law
    assert (self.twin_law is not None)
    self.twin_law=sgtbx.rt_mx(self.twin_law, r_den=24,t_den=288 )
    if self.twin_law.r().determinant() != 1:
      raise Sorry("The determinant of the provided twin law is not equal to unity")

    # make sure we have intensities
    if self.miller_array.is_real_array():
      if not self.miller_array.is_xray_intensity_array():
        self.miller_array = self.miller_array.f_as_f_sq()
    assert self.miller_array.is_xray_intensity_array()

  def detwin_it(self,alpha):
    print alpha
    print >> self.out, "Detwinning the data with fraction %3.2f"%(alpha)

    assert alpha is not None
    assert alpha<0.5
    assert alpha>=0.0

    detwin_object = mmtbx.scaling.detwin(self.miller_array.indices(),
                                         self.miller_array.data(),
                                         self.miller_array.sigmas(),
                                         self.miller_array.space_group(),
                                         self.miller_array.anomalous_flag(),
                                         self.twin_law.r().as_double() )
    detwin_object.detwin_with_alpha( alpha )
    new_intensities = detwin_object.detwinned_i()
    new_sigmas = detwin_object.detwinned_sigi()
    new_hkl = detwin_object.detwinned_hkl()
    new_miller_array =  self.miller_array.customized_copy(
      indices = new_hkl,
      data =  new_intensities,
      sigmas = new_sigmas ).set_observation_type( self.miller_array )
    return new_miller_array


class massage_data(object):
  def __init__(self,
               miller_array,
               parameters,
               out=None,
               n_residues=100,
               n_bases=0):

    self.params=parameters
    self.miller_array=miller_array.deep_copy().set_observation_type(miller_array).merge_equivalents().array()
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
      aniso_scale_and_b.show(out=out,verbose=1)
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
      print >> self.out, "Non-outliers found by the basic wilson statistics"
      print >> self.out, "protocol will be written out."
      basic_array = outlier_manager.basic_wilson_outliers(
        p_basic_wilson = self.params.outlier.parameters.basic_wilson.level,
        return_data = True)
      self.new_miller_array = basic_array

    if self.params.outlier.action == "extreme":
      print >> self.out, "Non-outliers found by the extreme value wilson statistics"
      print >> self.out, "protocol will be written out."
      extreme_array = outlier_manager.extreme_wilson_outliers(
      p_extreme_wilson = self.params.outlier.parameters.extreme_wilson.level,
      return_data = True)
      self.new_miller_array = extreme_array

    if self.params.outlier.action == "beamstop":
      print >> self.out, "Outliers found for the beamstop shadow"
      print >> self.out, "problems detection protocol will be written out."
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
      if self.params.symmetry.twinning_parameters.fraction is None:
        raise Sorry("Twin fraction not specified, not twinning data")
      twinner = twin_data(miller_array = self.new_miller_array,
                            twin_law = self.params.symmetry.twinning_parameters.twin_law,
                            out = self.out)
      self.final_array = twinner.twin_it(alpha=self.params.symmetry.twinning_parameters.fraction)


    if self.params.symmetry.action == "detwin":
      if self.params.symmetry.twinning_parameters.fraction is None:
        raise Sorry("Twin fraction not specified, not detwinning data")

      detwinner = detwin_data(miller_array = self.new_miller_array,
                            twin_law = self.params.symmetry.twinning_parameters.twin_law,
                            out = self.out)
      self.final_array = detwinner.detwin_it(alpha=self.params.symmetry.twinning_parameters.fraction)


    assert self.final_array is not None

  def return_data(self):
    return self.final_array

  def write_data(self):
    ## write out this miller array as sca if directed to do so:
    output_file=self.params.hklout
    n=len(output_file)
    auto_output_type=output_file[n-3:n]
    output_type = self.params.hklout_type
    if output_type == "mtz_or_sca":
      if auto_output_type in ["mtz","sca"]:
        output_type = auto_output_type
      else:
        raise Sorry("Unknown or unsupported output type")

    if output_type == "sca":
      import iotbx.scalepack.merge
      iotbx.scalepack.merge.write(
        file_name=output_file,miller_array=self.final_array)
    if output_type == "mtz":
      base_label=None
      if self.final_array.is_xray_intensity_array():
        base_label = "I"
      if self.final_array.is_xray_amplitude_array():
        base_label = "F"
      mtz_dataset = self.final_array.as_mtz_dataset(
        column_root_label=base_label+self.params.label_extension)
      mtz_dataset.mtz_object().write(output_file)
