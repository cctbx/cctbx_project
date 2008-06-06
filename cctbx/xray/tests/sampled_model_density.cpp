#include <cctbx/xray/sampled_model_density.h>
#include <scitbx/array_family/simple_io.h>
#include <iostream>

using namespace cctbx;

int main() {
  uctbx::unit_cell uc(af::double6(2., 3., 4., 80., 91., 105.));
  typedef xray::sampled_model_density<double> smd_t;
  smd_t::grid_point_type fft_n_real(20, 25, 36), fft_m_real(20, 25, 38);
  xray::scattering_type_registry scatt_t_registry;
  af::shared<xray::scatterer<> > scatterers;
  {
    typedef fractional<double> frac_t;
    typedef scitbx::sym_mat3<double> u_star_t;
    xray::scatterer<> h1("H1", frac_t(0.4206, 0.2589, 0.5113), 0.2849, 0.82,
                         "H", 0, 0);
    h1.u_star = adptbx::u_cart_as_u_star(
      uc, u_star_t(0.282, 0.185, 0.281, -0.010, 0.007, -0.006));
    scatterers.push_back(h1);

    xray::scatterer<> h2("H2", frac_t(0.4049, 0.7838, 0.3033), 0.1739, 0.54,
                         "H", 0, 0);
    h2.u_star =  adptbx::u_cart_as_u_star(
      uc, u_star_t(0.191, 0.117, 0.032, 0.050, 0.021, -0.045));
    scatterers.push_back(h2);

    xray::scatterer<> h3("H3", frac_t(0.4766, 0.5834, 0.9081), 0.1352, 0.97,
                         "H", 0, 0);
    h3.u_star =  adptbx::u_cart_as_u_star(
      uc, u_star_t(0.079, 0.119, 0.233, 0.016, 0.023, 0.057));
    scatterers.push_back(h3);

    xray::scatterer<> h4("H4", frac_t(0.5047, 0.2818, 0.7558), 0.1981, 0.60,
                         "H", 0, 0);
    h4.u_star = adptbx::u_cart_as_u_star(
      uc, u_star_t(0.086, 0.073, 0.095, 0.050, 0.003, -0.022));
    scatterers.push_back(h4);

    xray::scatterer<> h5("H5", frac_t(0.6184, 0.2505, 0.9097), 0.2989, 0.50,
                         "H", 0, 0);
    h5.u_star = adptbx::u_cart_as_u_star(
      uc, u_star_t(0.192, 0.208, 0.171, 0.035, -0.016, -0.018));
    scatterers.push_back(h5);

    xray::scatterer<> h6("H6", frac_t(0.9828, 0.8102, 0.9022), 0.2751, 0.69,
                         "H", 0, 0);
    h6.u_star = adptbx::u_cart_as_u_star(
      uc, u_star_t(0.149, 0.095, 0.132, -0.011, -0.031, -0.008));
    scatterers.push_back(h6);

    xray::scatterer<> h7("H7", frac_t(0.3101, 0.7298, 0.8988), 0.2380, 0.93,
                         "H", 0, 0);
    h7.u_star = adptbx::u_cart_as_u_star(
      uc, u_star_t(0.244, 0.125, 0.148, 0.006, -0.030, 0.106));
    scatterers.push_back(h7);

    xray::scatterer<> h8("H8", frac_t(0.6840, 0.4721, 0.1007), 0.0247, 0.63,
                         "H", 0, 0);
    h8.u_star = adptbx::u_cart_as_u_star(
      uc, u_star_t(0.185, 0.098, 0.129, 0.005, 0.047, -0.017));
    scatterers.push_back(h8);
  }
  scatt_t_registry.process(scatterers.const_ref());
  xray::sampled_model_density<double> smd(
    uc, scatterers.const_ref(),
    scatt_t_registry,
    fft_n_real, fft_m_real,
    0.033773727880779258, //u_base
    1e-8, // wing_cutoff
    -100, // exp_table_one_over_step_size
    false, // force_complex
    false, // sampled_density_must_be_positive
    1e-5, // tolerance_positive_definite
    false // use_u_base_as_u_extra
    // store_grid_indices_for_each_scatterer
  );
  std::cout << smd.real_map().as_1d().const_ref() << std::endl;
}
