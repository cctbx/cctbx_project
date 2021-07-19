#ifndef SIMTBX_PIXEL_STATS_H
#define SIMTBX_PIXEL_STATS_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/math/mean_and_variance.h>

namespace simtbx {
namespace pixel {

namespace af = scitbx::af;

struct pixel_stats {
  inline
  pixel_stats(){
  }

  void set_whitelist(
    af::shared<double> a,af::shared<double> b,af::shared<double> c){
    lunus_filtered_data = a;
    exp_data = b;
    sim_mock = c;
    SCITBX_ASSERT(lunus_filtered_data.size()==exp_data.size());
    SCITBX_ASSERT(exp_data.size()==sim_mock.size());
  }

  void set_shoebox_iterator(
    af::shared<int> a,af::shared<int> b, af::shared<std::size_t> c){
    shoebox_offset = a;
    shoebox_size = b;
    spots_pixels = c;
    SCITBX_ASSERT(shoebox_offset.size()==shoebox_size.size());
  }

  void analyze3(
    af::shared<double> whitelist_kernel_model, af::shared<double> reference_shoebox_sums,
    const int& slow_size, const int& panel_size, const double& keV_per_photon){

    //whitelist_kernel_model is the 1d list of pixel value of interest, in the order of interest
    std::size_t* spptr = spots_pixels.begin();

    af::shared<double> proposal_shoebox_mean_Z;
    af::shared<double> proposal_shoebox_sigma_Z;
    proposal_ctr_of_mass = af::shared<scitbx::vec3<double> >();

    af::shared<double> all_Z_values;
    LLG = 0;
    int firstpass_idx = -1;
    int whiteidx = -1;
    double* white_ptr = whitelist_kernel_model.begin();

    for (int sidx=0; sidx<shoebox_offset.size(); ++sidx){ //loop through the shoeboxes

      //first pass through shoebox, add all shoebox pixels to get proposal sum
      double proposal_shoebox_sum = 0.;
      for (int pidx=shoebox_offset[sidx]; pidx<shoebox_offset[sidx]+shoebox_size[sidx]; ++pidx){
        firstpass_idx += 1;
        double proposal_value = std::max(0.1,white_ptr[firstpass_idx]);
        white_ptr[firstpass_idx] = proposal_value;
        proposal_shoebox_sum += proposal_value;
      }
      SCITBX_ASSERT(firstpass_idx <= whitelist_kernel_model.size());
      double refrence_shoebox_sum = reference_shoebox_sums[sidx];
      double scale_factor = refrence_shoebox_sum / proposal_shoebox_sum;
      //result verified up to here

      //second pass through shoebox, determine center of mass, Z-statistics
      af::shared<double> shoebox_Z_values;
      scitbx::vec2<double> SUM_VEC(0,0);
      double SUM_wt = 0.;
      for (int pidx=shoebox_offset[sidx]; pidx<shoebox_offset[sidx]+shoebox_size[sidx]; ++pidx){
        int idxpx = spptr[pidx];
        int panelpx = idxpx%panel_size;
        int islow = panelpx/slow_size;
        int ifast = panelpx%slow_size;

        whiteidx += 1;
        double renormalized_proposal = white_ptr[whiteidx] * scale_factor;
        SUM_VEC = SUM_VEC + (renormalized_proposal * scitbx::vec2<double>(islow,ifast));
        SUM_wt += renormalized_proposal;
        double renormalize_bragg_plus_background =
          lunus_filtered_data[whiteidx] + renormalized_proposal;

        double experimental_pixel = exp_data[whiteidx];
        double abs_exp_pixel = std::abs(experimental_pixel);
        double abs_exp_pixel_photons = abs_exp_pixel/keV_per_photon;
        double poisson_noise_sigma = std::sqrt(abs_exp_pixel_photons);
        double std_dev_denominator_photons = poisson_noise_sigma;
        double renormalize_bragg_plus_background_photons =
          renormalize_bragg_plus_background/keV_per_photon;
        double mock_model_photons = sim_mock[whiteidx]/keV_per_photon;
        double diff_pixel_photons =
          renormalize_bragg_plus_background_photons - mock_model_photons;
        double Z = (diff_pixel_photons/std_dev_denominator_photons);
        shoebox_Z_values.push_back(Z);
        all_Z_values.push_back(Z);
if (renormalize_bragg_plus_background_photons <=0.0) {continue;} //cannot take log of negative value
        double pixel_ll = renormalize_bragg_plus_background_photons -
          mock_model_photons * std::log(renormalize_bragg_plus_background_photons); //sauter_eq_15_likelihood
        LLG+=pixel_ll;
      }
      SCITBX_ASSERT(whiteidx <= whitelist_kernel_model.size());
      // properties calculated per shoebox
      scitbx::math::mean_and_variance<double> stats(shoebox_Z_values.const_ref());
      proposal_shoebox_mean_Z.push_back( stats.mean() );
      proposal_shoebox_sigma_Z.push_back( stats.unweighted_sample_standard_deviation() );
      scitbx::vec2<double> c_o_m = SUM_VEC/SUM_wt;
      // there is a half pixel offset in our understanding of position
      proposal_ctr_of_mass.push_back(scitbx::vec3<double> (c_o_m[1]+0.5,c_o_m[0]+0.5,0.0));
    }
    // properties calculated over the entire set of shoeboxes
    // in python: self.refl_table["temp_values"] = proposal_ctr_of_mass
    // in python: rmsd = self.simple_rmsd(calc_data="temp_values",plot=False)
    scitbx::math::mean_and_variance<double> stats(all_Z_values.const_ref());
    mnz = stats.mean();
    sgz = stats.unweighted_sample_standard_deviation();
  }


  double get_LLG(){return LLG;}
  double get_mnz(){return mnz;}
  double get_sgz(){return sgz;}
  af::shared<scitbx::vec3<double> > get_proposal_center_of_mass(){return proposal_ctr_of_mass;}

  af::shared<double> lunus_filtered_data, exp_data, sim_mock;
  af::shared<int> shoebox_offset, shoebox_size;
  af::shared<std::size_t> spots_pixels;
  double LLG, rmsd, mnz/* mean Z over shoeboxes */, sgz /* std dev */;
  af::shared<scitbx::vec3<double> > proposal_ctr_of_mass;

};
} // pixel
} // simtbx
#endif // SIMTBX_PIXEL_STATS_H
