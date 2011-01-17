//! Peter Zwart April 05, 2005
#ifndef MMTBX_SCALING_SCALING_H
#define MMTBX_SCALING_SCALING_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cctbx/miller/sym_equiv.h>
#include <scitbx/math/chebyshev.h>

namespace mmtbx { namespace scaling {

  static const double d_star_sq_low_limit = 0.008;
  static const double d_star_sq_high_limit = 0.690;

  //! Returns low resolution limit of gamma lookup array
  inline
  double
  get_d_star_sq_low_limit()
  {
    return(d_star_sq_low_limit);
  }

  //! Returns low resolution limit of gamma lookup array
  inline
  double
  get_d_star_sq_high_limit()
  {
    return(d_star_sq_high_limit);
  }

  //! The gamma_prot look up table
  /*! describing the variation in a Wilson plot
   *  resolution ranges: bin 0=0.008 A^-2, increment=0.003478A^-2
   *
   *  This table only makes sense when you have only protein.
   *  Protein DNA/RNA complexes do not follow this table
   *
   *  See Zwart & Lamzin, Acta Cryst D60,220-226, 2004
   *  for details. This curve has been obtained from experimental data
   *  rather then from atomic coordinates.
   */

  static const double gamma_prot_lookup [] = {
    -0.340167,-0.349085,-0.432714,-0.537415,-0.597142,-0.612421,-0.585563,
    -0.5405  ,-0.471652,-0.345664,-0.180805,-0.016259, 0.118749, 0.200193,
     0.273765, 0.308676, 0.296059, 0.293262, 0.30348 , 0.340488, 0.360010,
     0.332146, 0.265772, 0.232256, 0.228609, 0.182422, 0.114466, 0.063855,
     0.0282269,0.005662,-0.023313,-0.054328,-0.086619,-0.141222,-0.165737,
    -0.164658,-0.163909,-0.186208,-0.213532,-0.215588,-0.213031,-0.219341,
    -0.219647,-0.20684 ,-0.200495,-0.195989,-0.188317,-0.157044,-0.124482,
    -0.121555,-0.110698,-0.093429,-0.086539,-0.083881,-0.071929,-0.067351,
    -0.056303,-0.031620,-0.006373,-0.001948,-0.005484,-0.013537,-0.012858,
    -0.006412,-0.000812,-0.008278,-0.037985,-0.063102,-0.048759,-0.037321,
    -0.038132,-0.035638,-0.037332,-0.052133,-0.073953,-0.096395,-0.079632,
    -0.060959,-0.082468,-0.113025,-0.131687,-0.138621,-0.185577,-0.224368,
    -0.224954,-0.203820,-0.192993,-0.205107,-0.250665,-0.248667,-0.233933,
    -0.245783,-0.248550,-0.239248,-0.234695,-0.244409,-0.239380,-0.236651,
    -0.264919,-0.266553,-0.252436,-0.251786,-0.241578,-0.236372,-0.230764,
    -0.225601,-0.216291,-0.198012,-0.215888,-0.230430,-0.222284,-0.217729,
    -0.214744,-0.225040,-0.245056,-0.239166,-0.227135,-0.232121,-0.225811,
    -0.224273,-0.240172,-0.258692,-0.264998,-0.257350,-0.247191,-0.231820,
    -0.220285,-0.228137,-0.211337,-0.185410,-0.173713,-0.187623,-0.203179,
    -0.217712,-0.246310,-0.246268,-0.219625,-0.194464,-0.187950,-0.193789,
    -0.197559,-0.183625,-0.169733,-0.172527,-0.192799,-0.181501,-0.137027,
    -0.117329,-0.133849,-0.161875,-0.164322,-0.166884,-0.170102,-0.165672,
    -0.152898,-0.139653,-0.134125,-0.142775,-0.134868,-0.108417,-0.096221,
    -0.104541,-0.107085,-0.071845,-0.062470,-0.078524,-0.094327,-0.064350,
    -0.040639,-0.043437,-0.015628, 0.000506,-0.001390, 0.022633, 0.045009,
     0.059946, 0.071254, 0.066774, 0.055433, 0.054914, 0.087948, 0.112070,
     0.100363, 0.069766, 0.067690, 0.088543, 0.087586, 0.126542, 0.171640,
     0.150073, 0.124215, 0.126672, 0.165951, 0.191294, 0.170923, 0.119596,
     0.130650, 0.139998, 0.149757, 0.214022};






  //! Computes the gamma_prot value for a given resolution
  /*! Uses a simple minded linear interpolation scheme.
   *  Clearly not optimal, but sufficient for practical purposes.
   */
  template <typename FloatType>
  FloatType
  gamma_prot(FloatType const& d_star_sq)
  {

    SCITBX_ASSERT (d_star_sq > d_star_sq_low_limit);
    SCITBX_ASSERT (d_star_sq < d_star_sq_high_limit);

    typedef FloatType f_t;
    f_t start = 0.008;
    f_t increment = 0.003478;

    int bin_low = static_cast<int>(
      std::floor( (d_star_sq - start - increment/2.0)/increment +0.5));
    int bin_high = bin_low+1;

    f_t d_star_sq_low = bin_low*increment+start;
    f_t d_star_sq_high = bin_high*increment+start;
    f_t gamma_prot_low = gamma_prot_lookup[bin_low];
    f_t gamma_prot_high = gamma_prot_lookup[bin_high];

    f_t dx = d_star_sq_high-d_star_sq_low;
    f_t dy = gamma_prot_high-gamma_prot_low;


    f_t result = (dy/dx)*(d_star_sq - d_star_sq_low )+gamma_prot_low;

    return (result);
  }



  //! Returns an array of gamma_prot values
  /*!
   * !Be carefull!, if outside of reso limit of lookup table,
   * zeros are returned.
   * It makes more sense to handle this at the stage of
   * preparing data then at this point.
   */
  template <typename FloatType>
  scitbx::af::shared<FloatType>
  get_gamma_prot(scitbx::af::const_ref<FloatType> const& d_star_sq)
  {

    scitbx::af::shared<FloatType> gamma_prot_array(d_star_sq.size(),0);

    for (unsigned i=0; i<d_star_sq.size();i++){

      if (d_star_sq[i]>d_star_sq_low_limit){
        if (d_star_sq[i]<d_star_sq_high_limit){
          gamma_prot_array[i] =  gamma_prot(d_star_sq[i]);
        }
      }

    }
    return(gamma_prot_array);

  }



  static const double a_h [] = {
    -0.117103661459, 0.00934859438867, 0.270068598623,
    0.284341399492, 0.552871721265, 0};
  static const double b_h [] = {
    3.05984661463, 0.746557743911, 3.29178616951,
    32.6456509125, 11.5463566492, 0};

  static const double a_c [] = { 2.657506, 1.078079, 1.490909,
                                -4.241070, 0.713791, 4.297983};
  static const double b_c [] = {14.780758, 0.776775, 42.086842,
                                -0.000294, 0.239535, 0};

  static const double a_n [] = {11.893780, 3.277479, 1.858092,
                                0.858927, 0.912985, -11.804902};
  static const double b_n [] = {0.000158, 10.232723, 30.344690,
                                0.656065, 0.217287, 0};

  static const double a_o [] = {2.960427, 2.508818, 0.637853,
                                0.722838, 1.142756, 0.027014};
  static const double b_o [] = {14.182259, 5.936858, 0.112726,
                                34.958481, 0.390240, 0};



  //!Returns the estimated summed square of form factors given the number of residues
  /*!
   *  For simplicity, I copied numbers from
   *      cctbx/eltbx/xray_scattering/wk1995.cpp
   *  rather then calling some functionality that i do not understand.
   *  This might change to something more efficient though.
   */

  template <typename FloatType>
  FloatType
  sigma_prot_sq(FloatType const& d_star_sq, FloatType const& n_residues)
  {
    typedef FloatType f_t;
    f_t result;
    f_t f_h_sq=0, f_c_sq=0, f_n_sq=0,f_o_sq=0;
    for (unsigned i=0; i<6;i++){
      f_h_sq +=a_h[i]*std::exp(-b_h[i]*d_star_sq/4.0);
      f_c_sq +=a_c[i]*std::exp(-b_c[i]*d_star_sq/4.0);
      f_n_sq +=a_n[i]*std::exp(-b_n[i]*d_star_sq/4.0);
      f_o_sq +=a_o[i]*std::exp(-b_o[i]*d_star_sq/4.0);
    }

    f_h_sq *= f_h_sq*8.0;
    f_c_sq *= f_c_sq*5.0;
    f_n_sq *= f_n_sq*1.5;
    f_o_sq *= f_o_sq*1.2;

    result = n_residues*(f_h_sq+f_c_sq+f_n_sq+f_o_sq);

    return(result);
  }

  //! Returns an array with sum of squared form factors
  template <typename FloatType>
  scitbx::af::shared<FloatType>
  get_sigma_prot_sq(scitbx::af::const_ref<FloatType> const& d_star_sq,
                    FloatType const& n_residues)
  {
    scitbx::af::shared<FloatType> sigma_prot_sq_array (d_star_sq.size(),0);
    for (unsigned i=0; i<d_star_sq.size(); i++){
      sigma_prot_sq_array[i] = sigma_prot_sq(d_star_sq[i],n_residues);
    }
    return(sigma_prot_sq_array);

  }




}}  // namespace mmtbx::scaling
#endif // MMTBX_SCALING_SCALING_H
