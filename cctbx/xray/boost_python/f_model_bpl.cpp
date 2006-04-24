#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/f_model.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace xray { namespace f_model { namespace boost_python {
  namespace {


    struct f_model_derivative_holder_wrappers
    {
      typedef f_model_derivative_holder<double> w_t;

      static void
      wrap()
      {
        using namespace boost::python;
        class_<w_t>("f_model_derivative_holder")//, no_init)
          .def("ksol",(double(w_t::*)())  &w_t::ksol )
          .def("ksol",(void (w_t::*)(double))  &w_t::ksol )

          .def("usol",(double(w_t::*)())  &w_t::usol )
          .def("usol",(void(w_t::*)(double))  &w_t::usol )

          .def("kpart",(double(w_t::*)()) &w_t::kpart )
          .def("kpart",(void(w_t::*)(double)) &w_t::kpart )

          .def("upart",(double(w_t::*)()) &w_t::upart )
          .def("upart",(void(w_t::*)(double)) &w_t::upart )

          .def("koverall", (double(w_t::*)()) &w_t::koverall )
          .def("koverall", (void(w_t::*)(double)) &w_t::koverall )

          .def("ustar", (scitbx::sym_mat3<double>(w_t::*)()) &w_t::ustar )
          .def("ustar", (void(w_t::*)(scitbx::sym_mat3<double>)) &w_t::ustar )

          .def("accumulate", &w_t::accumulate)
          ;
      }


    };




    struct f_model_wrappers
    {
      typedef f_model<double> w_t;

      static void
      wrap()
      {
        using namespace boost::python;

        class_<w_t>("f_model", no_init)
          .def(init<
               scitbx::af::const_ref< cctbx::miller::index<> > const&,  // 1 indices
               scitbx::af::const_ref< std::complex< double > > const&,  // 2 f_atmos
               scitbx::af::const_ref< std::complex< double > > const&,  // 3 f_mask
               cctbx::uctbx::unit_cell const&,                          // 4 unit cell
               double const&,                                           // 7 overall scale
               scitbx::sym_mat3< double > const&,                       // 8 u star
               double const&,                                           // 9 bs scale
               double const&,                                           // 0 bs u
               scitbx::af::const_ref< std::complex< double > > const&,  // 1 f_mask
               double const&,                                           // 2 bs scale
               double const&                                            // 3 bs u
               >
               ((arg_("hkl"),
                 arg_("f_atoms"),
                 arg_("f_mask"),
                 arg_("unit_cell"),
                 arg_("k_overall"),
                 arg_("u_star"),
                 arg_("k_sol"),
                 arg_("u_sol"),
                 arg_("f_part"),
                 arg_("k_part"),
                 arg_("u_part") )))

          .def(init<
               scitbx::af::const_ref< cctbx::miller::index<> > const&,  // 1 indices
               scitbx::af::const_ref< std::complex< double > > const&,  // 2 f_atmos
               scitbx::af::const_ref< std::complex< double > > const&,  // 3 f_mask
               cctbx::uctbx::unit_cell const&,                          // 4 unit cell
               double const&,                                           // 7 overall scale
               scitbx::sym_mat3< double > const&,                       // 8 u star
               double const&,                                           // 9 bs scale
               double const&,                                           // 0 bs b
               scitbx::af::const_ref< std::complex< double > > const&,  // 1 f_mask
               double const&,                                           // 2 bs scale
               double const&                                            // 3 bs b
               >
               ((arg_("hkl"),
                 arg_("f_atoms"),
                 arg_("f_mask"),
                 arg_("d_star_sq"),
                 arg_("k_overall"),
                 arg_("u_star"),
                 arg_("k_sol"),
                 arg_("u_sol"),
                 arg_("f_part"),
                 arg_("k_part"),
                 arg_("u_part") )))

          // derivatives using gradientflags
          .def("d_target_d_all",
               ( f_model_derivative_holder<double>(w_t::*)
                 (double const&,
                  double const&,
                  std::size_t const&,
                  scitbx::af::const_ref<bool> const& )
               ) &w_t::d_target_d_all)
          .def("d_target_d_all",
               ( f_model_derivative_holder<double>(w_t::*)
                 (scitbx::af::const_ref<double> const&,
                  scitbx::af::const_ref<double> const&,
                  scitbx::af::const_ref<bool> const&)
               ) &w_t::d_target_d_all)
          /*
          .def("d_target_d_all",
               ( f_model_derivative_holder<double>(w_t::*)
                 (scitbx::af::const_ref<std::complex<double> > const&,
                  scitbx::af::const_ref<bool> const&)
               ) &w_t::d_target_d_all)
          */
          .def("d_target_d_all",
               ( f_model_derivative_holder<double>(w_t::*)
                 (scitbx::af::const_ref<double> const&,
                  scitbx::af::const_ref<bool> const&)
                 ) &w_t::d_target_d_all)


          // updating variables
          .def("refresh", &w_t::refresh)
          .def("renew_fatoms", &w_t::renew_fatoms)
          .def("renew_fmask", &w_t::renew_fatoms)
          .def("renew_fpart", &w_t::renew_fpart)
          .def("renew_overall_scale_parameters",
               &w_t::renew_overall_scale_parameters)
          .def("renew_bulk_solvent_scale_parameters",
               &w_t::renew_bulk_solvent_scale_parameters)
          .def("renew_partial_structure_scale_parameters",
               &w_t::renew_partial_structure_scale_parameters)
          // get individual scales and sf's
          .def("f_atoms", &w_t::f_atoms)
          .def("f_bulk", &w_t::f_bulk)
          .def("f_part", &w_t::f_part)
          .def("bulk_scale", &w_t::bulk_scale)
          .def("part_scale", &w_t::part_scale)
          .def("f_model", &w_t::get_f_model)
          .def("overall_scale", &w_t::overall_scale)
          .def("fb_cart", &w_t::fb_cart)
          .def("d_f_model_d_f_atoms", &w_t::d_f_model_d_f_atoms)

          .def("ustar", (scitbx::sym_mat3<double>(w_t::*)()) &w_t::ustar)
          .def("koverall",(double(w_t::*)()) &w_t::koverall)
          .def("ksol", (double(w_t::*)()) &w_t::ksol)
          .def("kpart",(double(w_t::*)()) &w_t::kpart)
          .def("usol", (double(w_t::*)()) &w_t::usol)
          .def("upart",(double(w_t::*)()) &w_t::upart)

          .def("ustar", (void(w_t::*)(scitbx::sym_mat3<double>)) &w_t::ustar)
          .def("koverall",(void(w_t::*)(double)) &w_t::koverall)
          .def("ksol", (void(w_t::*)(double)) &w_t::ksol)
          .def("kpart",(void(w_t::*)(double)) &w_t::kpart)
          .def("usol", (void(w_t::*)(double)) &w_t::usol)
          .def("upart",(void(w_t::*)(double)) &w_t::upart)

          .def("d_star_sq", &w_t::d_star_sq)
          // selection method
          .def("select", &w_t::select)
          ;
      }
    };

  }

}} // namespace cctbx::xray

  namespace boost_python {

    void wrap_f_model()
    {
      f_model::boost_python::f_model_wrappers::wrap();
      f_model::boost_python::f_model_derivative_holder_wrappers::wrap();
    }

  }

}}
