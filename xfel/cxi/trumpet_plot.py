from __future__ import absolute_import, division, print_function
from six.moves import range
from matplotlib import pyplot as plt
from matplotlib import gridspec
import math
import os
from dials.array_family import flex

class trumpet_plot(object):

  def __init__(self, nrows=2, ncols=2):
    self.AD1TF7B_MAX2T = 45.
    self.AD1TF7B_MAXDP = 0.5
    self.AD1TF7B_MAXDP = 0.1
    self.color_encoding = "I/sigma"
    plt.figure(figsize=(12,10))  #use 9 for laptop
    self.gs = gridspec.GridSpec(nrows=nrows,ncols=ncols,width_ratios=[3,1])
    self.ncols = ncols

  def set_reduction(self,reduction): #see data_utils for definition of this container
    self.reduction = reduction
  def set_refined(self,info): # a repository for the refined parameters
    self.refined = info

  def plot_one_model(self,nrow,out):
    fig = plt.subplot(self.gs[nrow*self.ncols])
    two_thetas = self.reduction.get_two_theta_deg()
    degrees = self.reduction.get_delta_psi_deg()

    if self.color_encoding=="conventional":
          positive = (self.reduction.i_sigi>=0.)
          fig.plot(two_thetas.select(positive), degrees.select(positive), "bo")
          fig.plot(two_thetas.select(~positive), degrees.select(~positive), "r+")
    elif self.color_encoding=="I/sigma":
          positive = (self.reduction.i_sigi>=0.)
          tt_selected = two_thetas.select(positive)
          dp_selected = degrees.select(positive)
          i_sigi_select = self.reduction.i_sigi.select(positive)
          order = flex.sort_permutation(i_sigi_select)
          tt_selected = tt_selected.select(order)
          dp_selected = dp_selected.select(order)
          i_sigi_selected = i_sigi_select.select(order)
          from matplotlib.colors import Normalize
          dnorm = Normalize()
          dcolors = i_sigi_selected.as_numpy_array()
          dnorm.autoscale(dcolors)
          N = len(dcolors)
          CMAP = plt.get_cmap("rainbow")
          if self.refined.get("partiality_array",None) is None:
            for n in range(N):
              fig.plot([tt_selected[n]],[dp_selected[n]],
              color=CMAP(dnorm(dcolors[n])),marker=".", markersize=10)
          else:
            partials = self.refined.get("partiality_array")
            partials_select = partials.select(positive)
            partials_selected = partials_select.select(order)
            assert len(partials)==len(positive)
            for n in range(N):
              fig.plot([tt_selected[n]],[dp_selected[n]],
              color=CMAP(dnorm(dcolors[n])),marker=".", markersize=20*partials_selected[n])
              # change the markersize to indicate partiality.
          negative = (self.reduction.i_sigi<0.)
          fig.plot(two_thetas.select(negative), degrees.select(negative), "r+", linewidth=1)
    else:
          strong = (self.reduction.i_sigi>=10.)
          positive = ((~strong) & (self.reduction.i_sigi>=0.))
          negative = (self.reduction.i_sigi<0.)
          assert (strong.count(True)+positive.count(True)+negative.count(True) ==
                  len(self.reduction.i_sigi))
          fig.plot(two_thetas.select(positive), degrees.select(positive), "bo")
          fig.plot(two_thetas.select(strong), degrees.select(strong), marker='.',linestyle='None',
           markerfacecolor='#00ee00', markersize=10)
          fig.plot(two_thetas.select(negative), degrees.select(negative), "r+")

    # indicate the imposed resolution filter
    wavelength = self.reduction.experiment.beam.get_wavelength()
    imposed_res_filter = self.reduction.get_imposed_res_filter(out)
    resolution_markers = [
      a for a in [imposed_res_filter,self.reduction.measurements.d_min()] if a is not None]
    for RM in resolution_markers:
          two_th = (180./math.pi)*2.*math.asin(wavelength/(2.*RM))
          plt.plot([two_th, two_th],[self.AD1TF7B_MAXDP*-0.8,self.AD1TF7B_MAXDP*0.8],'k-')
          plt.text(two_th,self.AD1TF7B_MAXDP*-0.9,"%4.2f"%RM)

    #indicate the linefit
    mean = flex.mean(degrees)
    minplot = flex.min(two_thetas)
    plt.plot([0,minplot],[mean,mean],"k-")
    LR = flex.linear_regression(two_thetas, degrees)
    model_y = LR.slope()*two_thetas + LR.y_intercept()
    plt.plot(two_thetas, model_y, "k-")

    #Now let's take care of the red and green lines.
    half_mosaic_rotation_deg = self.refined["half_mosaic_rotation_deg"]
    mosaic_domain_size_ang = self.refined["mosaic_domain_size_ang"]
    red_curve_domain_size_ang = self.refined.get("red_curve_domain_size_ang",mosaic_domain_size_ang)
    a_step = self.AD1TF7B_MAX2T / 50.
    a_range = flex.double([a_step*x for x in range(1,50)]) # domain two-theta array
    #Bragg law [d=L/2sinTH]
    d_spacing = (wavelength/(2.*flex.sin(math.pi*a_range/360.)))
    # convert two_theta to a delta-psi.  Formula for Deffective [Dpsi=d/2Deff]
    inner_phi_deg = flex.asin((d_spacing / (2.*red_curve_domain_size_ang)) )*(180./math.pi)
    outer_phi_deg = flex.asin((d_spacing / (2.*mosaic_domain_size_ang)) + \
      half_mosaic_rotation_deg*math.pi/180. )*(180./math.pi)
    plt.title("ML: mosaicity FW=%4.2f deg, Dsize=%5.0fA on %d spots\n%s"%(
          2.*half_mosaic_rotation_deg, mosaic_domain_size_ang, len(two_thetas),
          os.path.basename(self.reduction.filename)))
    plt.plot(a_range, inner_phi_deg, "r-")
    plt.plot(a_range,-inner_phi_deg, "r-")
    plt.plot(a_range, outer_phi_deg, "g-")
    plt.plot(a_range, -outer_phi_deg, "g-")
    plt.xlim([0,self.AD1TF7B_MAX2T])
    plt.ylim([-self.AD1TF7B_MAXDP,self.AD1TF7B_MAXDP])

    #second plot shows histogram
    fig = plt.subplot(self.gs[1+nrow*self.ncols])
    plt.xlim([-self.AD1TF7B_MAXDP,self.AD1TF7B_MAXDP])
    nbins = 50
    n,bins,patches = plt.hist(dp_selected, nbins,
           range=(-self.AD1TF7B_MAXDP,self.AD1TF7B_MAXDP),
           weights=self.reduction.i_sigi.select(positive),
           normed=0, facecolor="orange", alpha=0.75)
    #ersatz determine the median i_sigi point:
    isi_positive = self.reduction.i_sigi.select(positive)
    isi_order = flex.sort_permutation(isi_positive)
    reordered = isi_positive.select(isi_order)
    isi_median = reordered[int(len(isi_positive)*0.9)]
    isi_top_half_selection = (isi_positive>isi_median)
    n,bins,patches = plt.hist(dp_selected.select(isi_top_half_selection), nbins,
           range=(-self.AD1TF7B_MAXDP,self.AD1TF7B_MAXDP),
           weights=isi_positive.select(isi_top_half_selection),
           normed=0, facecolor="#ff0000", alpha=0.75)
    plt.xlabel("(degrees)")
    plt.title("Weighted histogram of Delta-psi")

  def set_color_coding(self,choice):
    self.color_encoding = choice

  def show(self):
    plt.show()
        #plt.tight_layout()
        #plt.savefig('grid_figure.pdf')

def trumpet_wrapper(result, postx, file_name, params, out):
      from scitbx import matrix
      from dxtbx.model import BeamFactory
      beam = BeamFactory.make_beam(s0=(0,0,-1./result["wavelength"]))
      from dxtbx.model import Experiment
      from dxtbx.model import crystal

      obs_to_plot = postx.observations_original_index_pair1_selected # XXX uses a private interface

      HKL=obs_to_plot.indices()
      i_sigi=obs_to_plot.data()/obs_to_plot.sigmas()
      direct_matrix = result["current_orientation"][0].direct_matrix()
      real_a = direct_matrix[0:3]
      real_b = direct_matrix[3:6]
      real_c = direct_matrix[6:9]
      SG = obs_to_plot.space_group()
      crystal = crystal.crystal_model(real_a, real_b, real_c, space_group=SG)
      q = Experiment(beam=beam, crystal=crystal)

      TPL = trumpet_plot()
      from xfel.cxi.data_utils import reduction
      RED = reduction(filename = file_name, experiment = q, HKL = HKL, i_sigi = i_sigi,
                    measurements = obs_to_plot, params = params)
      TPL.set_reduction(RED)

      # first plot: before postrefinement
      TPL.reduction.experiment.crystal.set_A(matrix.sqr(result["current_orientation"][0].reciprocal_matrix()))
      TPL.set_refined(dict(half_mosaic_rotation_deg=result["ML_half_mosaicity_deg"][0],
                             mosaic_domain_size_ang=result["ML_domain_size_ang"][0]))
      TPL.plot_one_model(nrow=0,out=out)

      # second plot after postrefinement
      values = postx.get_parameter_values()
      TPL.reduction.experiment.crystal.set_A(postx.refined_mini.refinery.get_eff_Astar(values))
      TPL.refined["red_curve_domain_size_ang"]=1/values.RS
      TPL.refined["partiality_array"]=postx.refined_mini.refinery.get_partiality_array(values)
      TPL.plot_one_model(nrow=1,out=out)

      TPL.show()
