#include <boost/python.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec2.h>

namespace prime { namespace index_ambiguity {
  scitbx::af::shared<double>
  calc_BD_alg_2(scitbx::af::versa<double, scitbx::af::flex_grid<> > r_grid, scitbx::af::shared<scitbx::vec2<double> > x_vec_set) {
    //calculate residuals of the Brehm & Diederichs Algorithm 2
    scitbx::af::shared<double> r_out;
    for (int i = 0; i < x_vec_set.size() - 1; i++) {
      double sum_r = 0;
      for (int j = i+1; j < x_vec_set.size(); j++) {
        sum_r += std::abs(r_grid(i,j) - (x_vec_set[i] * x_vec_set[j]));
      }
      r_out.push_back(sum_r);
    }
    return r_out;
  }
  double
  calc_BD_alg_2_sum_sqr(scitbx::af::versa<double, scitbx::af::flex_grid<> > r_grid, scitbx::af::shared<scitbx::vec2<double> > x_vec_set) {
    //calculate squared residuals of the Brehm & Diederichs Algorithm 2
    double sum_r_sqr = 0;
    for (int i = 0; i < x_vec_set.size() - 1; i++) {
      for (int j = i+1; j < x_vec_set.size(); j++) {
        double tmp = r_grid(i,j) - (x_vec_set[i] * x_vec_set[j]);
        sum_r_sqr += (tmp * tmp);
      }
    }
    return sum_r_sqr;
  }
  scitbx::af::shared<double>
  calc_gradient_BD_alg_2(scitbx::af::versa<double, scitbx::af::flex_grid<> > r_grid, scitbx::af::shared<scitbx::vec2<double> > x_vec_set) {
    //calculate finite differences of each pair of parameters
    scitbx::af::shared<double> g_out;
    double f = calc_BD_alg_2_sum_sqr(r_grid, x_vec_set);
    double delta = 1.0e-7;
    for (int i = 0; i < x_vec_set.size(); i++) {
      //each component of the vector has its own gradient
      x_vec_set[i][0] += delta;
      double df_x = calc_BD_alg_2_sum_sqr(r_grid, x_vec_set);
      g_out.push_back((df_x - f)/delta);
      x_vec_set[i][0] -= delta;
      x_vec_set[i][1] += delta;
      double df_y = calc_BD_alg_2_sum_sqr(r_grid, x_vec_set);
      g_out.push_back((df_y - f)/delta);
      x_vec_set[i][1] -= delta;
    }
    return g_out;
  }
}}

namespace prime { namespace index_ambiguity {
namespace boost_python {
  void
  wrap_all() {
    boost::python::def("calc_BD_alg_2", calc_BD_alg_2);
    boost::python::def("calc_BD_alg_2_sum_sqr", calc_BD_alg_2_sum_sqr);
    boost::python::def("calc_gradient_BD_alg_2", calc_gradient_BD_alg_2);
  }
}}}

BOOST_PYTHON_MODULE(prime_index_ambiguity_ext)
{
  prime::index_ambiguity::boost_python::wrap_all();
}
