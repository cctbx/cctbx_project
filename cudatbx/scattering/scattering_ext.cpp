#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cudatbx/scattering/direct_summation.h>

namespace cudatbx { namespace scattering {

  struct direct_summation_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<cudatbx::scattering::direct_summation>("direct_summation",init<>() )
        .def("add",&cudatbx::scattering::direct_summation::add)
        .def("get_sum",&cudatbx::scattering::direct_summation::get_sum)
        .def("prepare_saxs",&cudatbx::scattering::direct_summation::prepare_saxs)
        .def("set_xyz",&cudatbx::scattering::direct_summation::set_xyz)
        .def("clear_xyz",&cudatbx::scattering::direct_summation::clear_xyz)
        .def("set_solvent_weights",
             &cudatbx::scattering::direct_summation::set_solvent_weights)
        .def("clear_solvent_weights",
             &cudatbx::scattering::direct_summation::clear_solvent_weights)
        .def("set_hkl",&cudatbx::scattering::direct_summation::set_hkl)
        .def("clear_hkl",&cudatbx::scattering::direct_summation::clear_hkl)
        .def("set_q",&cudatbx::scattering::direct_summation::set_q)
        .def("clear_q",&cudatbx::scattering::direct_summation::clear_q)
        .def("set_lattice",&cudatbx::scattering::direct_summation::set_lattice)
        .def("clear_weights",&cudatbx::scattering::direct_summation::clear_weights)
        .def("clear_lattice",&cudatbx::scattering::direct_summation::clear_lattice)
        .def("set_rotations_translations",
             &cudatbx::scattering::direct_summation::set_rotations_translations)
        .def("clear_rotations_translations",
             &cudatbx::scattering::direct_summation::clear_rotations_translations)
        .def("set_scattering_types",
             &cudatbx::scattering::direct_summation::set_scattering_types)
        .def("clear_scattering_types",
             &cudatbx::scattering::direct_summation::clear_scattering_types)
        .def("set_scattering_type_registry",
             &cudatbx::scattering::direct_summation::set_scattering_type_registry)
        .def("clear_scattering_type_registry",
             &cudatbx::scattering::direct_summation::clear_scattering_type_registry)
        .def("allocate_amplitudes",
             &cudatbx::scattering::direct_summation::allocate_amplitudes)
        .def("reset_amplitudes",&cudatbx::scattering::direct_summation::reset_amplitudes)
        .def("clear_amplitudes",&cudatbx::scattering::direct_summation::clear_amplitudes)
        .def("allocate_workspace",
             &cudatbx::scattering::direct_summation::allocate_workspace)
        .def("clear_workspace",&cudatbx::scattering::direct_summation::clear_workspace)
        .def("run_kernel",&cudatbx::scattering::direct_summation::run_kernel)
        .def("run_saxs_kernel",&cudatbx::scattering::direct_summation::run_saxs_kernel)
        .def("run_solvent_saxs_kernel",
             &cudatbx::scattering::direct_summation::run_solvent_saxs_kernel)
        .def("run_collect_solvent_saxs_kernel",
             &cudatbx::scattering::direct_summation::run_collect_solvent_saxs_kernel)
        .def("sum_over_lattice",
             &cudatbx::scattering::direct_summation::sum_over_lattice)
        ;
    }
  };

  }

  BOOST_PYTHON_MODULE(cudatbx_scattering_ext)
  {
    cudatbx::scattering::direct_summation_wrapper::wrap();
  }
}
