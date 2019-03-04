#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/rigu.h>
#include <boost/multi_array.hpp>

namespace cctbx { namespace adp_restraints {
namespace {

  struct rigu_proxy_wrappers
  {
    typedef rigu_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("rigu_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg("i_seqs"),
           arg("weight"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_rigu_proxy")
        ;
      }
    }
  };

  // original procedure for tests only
  //! Calculate all the partial derivatives of the adp
  //! in the local cartesian coordinate system.
  af::shared<rigu::sym_mat3d> CalcReferenceGradients(const rigu::mat3d& RM) {
    /** See Parois, et. al. (2018). J. Appl. Cryst. 51, 1059-1068.
    *  for the detail on how it works
    *
    *  The U tensors are vectorized as 9 elements vectors
    *  column by column to simplify the calculations.
    *  Operations like V = A U B translate then to vec(V) =
    *  (B^t @ A) U with @ the kronecker product
    */

    /** dUcart contains the 6 partial derivatives dU_cart/dUij_cart
    *  writen as vectors.
    *  The 6 columns vectors are then stacked column wise to form dUcart
    */
    const int dUcart[9][6] = {
      //dU11 dU22 dU33 dU12 dU13 dU23
      { 1,   0,   0,   0,   0,   0 }, // U11
      { 0,   0,   0,   1,   0,   0 }, // U21
      { 0,   0,   0,   0,   1,   0 }, // U23
      { 0,   0,   0,   1,   0,   0 }, // U21
      { 0,   1,   0,   0,   0,   0 }, // U22
      { 0,   0,   0,   0,   0,   1 }, // U23
      { 0,   0,   0,   0,   1,   0 }, // U31
      { 0,   0,   0,   0,   0,   1 }, // U32
      { 0,   0,   1,   0,   0,   0 }, // U33
    };

    boost::multi_array<double, 2> kron(boost::extents[9][9]);

    //! calulating the kronecker product
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        int startRow = i * 3;
        int startCol = j * 3;
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            kron[startRow + k][startCol + l] = RM(i, j)*RM(k, l);
          }
        }
      }
    }

    /** vec(dRUcart) = kron vec(dUcart).
    *  Derivatives of (RUcart = RM Ucart RM^t)
    */
    boost::multi_array<double, 2> dRUcart_(boost::extents[9][6]);
    std::fill(dRUcart_.data(), dRUcart_.data() + dRUcart_.num_elements(), 0);
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 6; j++) {
        for (int k = 0; k < 9; k++) {
          dRUcart_[i][j] += kron[i][k] * dUcart[k][j];
        }
      }
    }
    af::shared<rigu::sym_mat3d> dRUcart(3);
    for (int i = 0; i < 3; i++) {
      int r_idx[3] = { 8, 6, 7 };
      for (int j = 0; j < 6; j++) {
        dRUcart[i][j] = dRUcart_[r_idx[i]][j];
      }
    }
    return dRUcart;
  }
}

struct rigu_wrappers {
  typedef rigu w_t;

  static void wrap() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    class_<w_t>("rigu", no_init)
      .def(init<
        af::tiny<scitbx::vec3<double>, 2> const&,
        af::tiny<scitbx::sym_mat3<double>, 2> const&,
        double>(
        (arg("sites"),
          arg("u_cart"),
          arg("weight"))))
      .def(init<
        adp_restraint_params<double> const &,
        rigu_proxy const&>(
        (arg("params"),
          arg("proxy"))))
      .add_property("weight", make_getter(&w_t::weight, rbv()))
      .def("delta_33", &w_t::delta_33)
      .def("delta_13", &w_t::delta_13)
      .def("delta_23", &w_t::delta_23)
      .def("residual33", &w_t::residual33)
      .def("residual13", &w_t::residual13)
      .def("residual23", &w_t::residual23)
      .def("gradients33", &w_t::gradients33)
      .def("gradients13", &w_t::gradients13)
      .def("gradients23", &w_t::gradients23)
      .def("RM", &w_t::getRM)
      .def("raw_gradients", &w_t::raw_gradients)
      .def("reference_gradients", &CalcReferenceGradients)
      .staticmethod("reference_gradients")
      ;
  }
};

void wrap_all() {
  using namespace boost::python;
  rigu_wrappers::wrap();

  rigu_proxy_wrappers::wrap();
  def("rigu_residual_sum",
    adp_restraint_residual_sum_aniso<rigu_proxy, rigu>::impl,
    (arg("params"),
      arg("proxies"),
      arg("gradients_aniso_cart")));
  def("rigu_residuals",
    adp_restraint_residuals<rigu_proxy, rigu>::impl,
    (arg("params"),
      arg("proxies")));
}

namespace boost_python {

  void wrap_rigu() { wrap_all(); }

}}}
