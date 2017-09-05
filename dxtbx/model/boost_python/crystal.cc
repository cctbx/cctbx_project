/*
 * crystal.cc
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <string>
#include <sstream>
#include <dxtbx/model/crystal.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  static
  Crystal* make_crystal_default(
      const vec3<double> &real_space_a,
      const vec3<double> &real_space_b,
      const vec3<double> &real_space_c,
      const cctbx::sgtbx::space_group &space_group) {

    Crystal *crystal = new Crystal(
      real_space_a,
      real_space_b,
      real_space_c,
      space_group);

    return crystal;
  }

  static
  Crystal* make_crystal_with_symbol(
      const vec3<double> &real_space_a,
      const vec3<double> &real_space_b,
      const vec3<double> &real_space_c,
      const std::string &space_group_symbol) {

    Crystal *crystal = new Crystal(
      real_space_a,
      real_space_b,
      real_space_c,
      cctbx::sgtbx::space_group(
        cctbx::sgtbx::space_group_symbols(
          space_group_symbol)));

    return crystal;
  }

  static
  MosaicCrystalKabsch2010* make_kabsch2010_mosaic_crystal_default(
      const vec3<double> &real_space_a,
      const vec3<double> &real_space_b,
      const vec3<double> &real_space_c,
      const cctbx::sgtbx::space_group &space_group) {

    MosaicCrystalKabsch2010 *crystal = new MosaicCrystalKabsch2010(
      real_space_a,
      real_space_b,
      real_space_c,
      space_group);

    return crystal;
  }

  static
  MosaicCrystalKabsch2010* make_kabsch2010_mosaic_crystal_with_symbol(
      const vec3<double> &real_space_a,
      const vec3<double> &real_space_b,
      const vec3<double> &real_space_c,
      const std::string &space_group_symbol) {

    MosaicCrystalKabsch2010 *crystal = new MosaicCrystalKabsch2010(
      real_space_a,
      real_space_b,
      real_space_c,
      cctbx::sgtbx::space_group(
        cctbx::sgtbx::space_group_symbols(
          space_group_symbol)));

    return crystal;
  }

  static
  MosaicCrystalSauter2014* make_sauter2014_mosaic_crystal_default(
      const vec3<double> &real_space_a,
      const vec3<double> &real_space_b,
      const vec3<double> &real_space_c,
      const cctbx::sgtbx::space_group &space_group) {

    MosaicCrystalSauter2014 *crystal = new MosaicCrystalSauter2014(
      real_space_a,
      real_space_b,
      real_space_c,
      space_group);

    return crystal;
  }

  static
  MosaicCrystalSauter2014* make_sauter2014_mosaic_crystal_with_symbol(
      const vec3<double> &real_space_a,
      const vec3<double> &real_space_b,
      const vec3<double> &real_space_c,
      const std::string &space_group_symbol) {

    MosaicCrystalSauter2014 *crystal = new MosaicCrystalSauter2014(
      real_space_a,
      real_space_b,
      real_space_c,
      cctbx::sgtbx::space_group(
        cctbx::sgtbx::space_group_symbols(
          space_group_symbol)));

    return crystal;
  }

  static
  void Crystal_set_A_at_scan_points_from_tuple(CrystalBase &self, boost::python::tuple l) {
    scitbx::af::shared< mat3<double> > A_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      mat3<double> A = boost::python::extract< mat3<double> >(l[i]);
      A_list.push_back(A);
    }
    self.set_A_at_scan_points(A_list.const_ref());
  }

  static
  void Crystal_set_A_at_scan_points_from_list(CrystalBase &self, boost::python::list l) {
    scitbx::af::shared< mat3<double> > A_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      mat3<double> A = boost::python::extract< mat3<double> >(l[i]);
      A_list.push_back(A);
    }
    self.set_A_at_scan_points(A_list.const_ref());
  }

  static
  void Crystal_set_B_covariance_from_tuple(CrystalBase &self, boost::python::object obj) {
    scitbx::af::versa< double, scitbx::af::c_grid<2> > B_cov(scitbx::af::c_grid<2>(9,9));
    DXTBX_ASSERT(boost::python::len(obj) == 9*9);
    for (std::size_t i = 0; i < boost::python::len(obj); ++i) {
      B_cov[i] = boost::python::extract<double>(obj[i]);
    }
    self.set_B_covariance(B_cov.const_ref());
  }

  struct CrystalPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const Crystal &obj) {
      scitbx::af::shared< vec3<double> > real_space_v = obj.get_real_space_vectors();
      return boost::python::make_tuple(
          real_space_v[0],
          real_space_v[1],
          real_space_v[2],
          obj.get_space_group());
    }

    static
    boost::python::tuple getstate(boost::python::object obj)
    {
      const Crystal &crystal = boost::python::extract<const Crystal &>(obj)();
      return boost::python::make_tuple(
          obj.attr("__dict__"),
          crystal.get_A_at_scan_points(),
          crystal.get_B_covariance());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state)
    {
      Crystal &crystal = boost::python::extract<Crystal&>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 3);

      // restore the object's __dict__
      boost::python::dict d = boost::python::extract<boost::python::dict>(
          obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref< mat3<double> > A_list = boost::python::extract<
        scitbx::af::const_ref< mat3<double> > >(state[1]);
      scitbx::af::const_ref< double, scitbx::af::c_grid<2> > cov_B = boost::python::extract<
        scitbx::af::const_ref< double, scitbx::af::c_grid<2> > >(state[2]);
      crystal.set_A_at_scan_points(A_list);
      crystal.set_B_covariance(cov_B);
    }

    static bool getstate_manages_dict() { return true; }
  };

  struct MosaicCrystalKabsch2010PickleSuite: CrystalPickleSuite {
    static
    boost::python::tuple getinitargs(const MosaicCrystalKabsch2010 &obj) {
      scitbx::af::shared< vec3<double> > real_space_v = obj.get_real_space_vectors();
      return boost::python::make_tuple(
          real_space_v[0],
          real_space_v[1],
          real_space_v[2],
          obj.get_space_group());
    }

    static
    boost::python::tuple getstate(boost::python::object obj)
    {
      const MosaicCrystalKabsch2010 &crystal = boost::python::extract<const MosaicCrystalKabsch2010 &>(obj)();
      return boost::python::make_tuple(
          obj.attr("__dict__"),
          crystal.get_A_at_scan_points(),
          crystal.get_B_covariance(),
          crystal.get_mosaicity());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state)
    {
      MosaicCrystalKabsch2010 &crystal = boost::python::extract<MosaicCrystalKabsch2010&>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 4);

        // restore the object's __dict__
      boost::python::dict d = boost::python::extract<boost::python::dict>(
          obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref< mat3<double> > A_list = boost::python::extract<
        scitbx::af::const_ref< mat3<double> > >(state[1]);
      scitbx::af::const_ref< double, scitbx::af::c_grid<2> > cov_B = boost::python::extract<
        scitbx::af::const_ref< double, scitbx::af::c_grid<2> > >(state[2]);
      crystal.set_A_at_scan_points(A_list);
      crystal.set_B_covariance(cov_B);
      double mosaicity = boost::python::extract<double>(state[3]);
      crystal.set_mosaicity(mosaicity);
    }
  };

  struct MosaicCrystalSauter2014PickleSuite: CrystalPickleSuite {
    static
    boost::python::tuple getinitargs(const MosaicCrystalSauter2014 &obj) {
      scitbx::af::shared< vec3<double> > real_space_v = obj.get_real_space_vectors();
      return boost::python::make_tuple(
          real_space_v[0],
          real_space_v[1],
          real_space_v[2],
          obj.get_space_group());
    }

    static
    boost::python::tuple getstate(boost::python::object obj)
    {
      const MosaicCrystalSauter2014 &crystal = boost::python::extract<const MosaicCrystalSauter2014 &>(obj)();
      return boost::python::make_tuple(
          obj.attr("__dict__"),
          crystal.get_A_at_scan_points(),
          crystal.get_B_covariance(),
          crystal.get_half_mosaicity_deg(),
          crystal.get_domain_size_ang());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state)
    {
      MosaicCrystalSauter2014 &crystal = boost::python::extract<MosaicCrystalSauter2014&>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 5);

        // restore the object's __dict__
      boost::python::dict d = boost::python::extract<boost::python::dict>(
          obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref< mat3<double> > A_list = boost::python::extract<
        scitbx::af::const_ref< mat3<double> > >(state[1]);
      scitbx::af::const_ref< double, scitbx::af::c_grid<2> > cov_B = boost::python::extract<
        scitbx::af::const_ref< double, scitbx::af::c_grid<2> > >(state[2]);
      crystal.set_A_at_scan_points(A_list);
      crystal.set_B_covariance(cov_B);
      double half_mosaicity_deg = boost::python::extract<double>(state[3]);
      crystal.set_half_mosaicity_deg(half_mosaicity_deg);
      double domain_size_ang = boost::python::extract<double>(state[4]);
      crystal.set_domain_size_ang(domain_size_ang);
    }
  };

  void export_crystal()
  {
    class_ <CrystalBase, boost::noncopyable> ("CrystalBase", no_init)
      .def("set_unit_cell", &CrystalBase::set_unit_cell)
      .def("update_B", &CrystalBase::update_B)
      .def("set_U", &CrystalBase::set_U)
      .def("get_U", &CrystalBase::get_U)
      .def("set_B", &CrystalBase::set_B)
      .def("get_B", &CrystalBase::get_B)
      .def("set_A", &CrystalBase::set_A)
      .def("get_A", &CrystalBase::get_A)
      .def("get_unit_cell", &CrystalBase::get_unit_cell)
      .def("get_real_space_vectors", &CrystalBase::get_real_space_vectors)
      .def("set_space_group", &CrystalBase::set_space_group)
      .def("get_space_group", &CrystalBase::get_space_group)
      .add_property("num_scan_points", &CrystalBase::get_num_scan_points)
      .def("get_num_scan_points", &CrystalBase::get_num_scan_points)
      .def("set_A_at_scan_points", &CrystalBase::set_A_at_scan_points)
      .def("set_A_at_scan_points", &Crystal_set_A_at_scan_points_from_tuple)
      .def("set_A_at_scan_points", &Crystal_set_A_at_scan_points_from_list)
      .def("get_A_at_scan_point", &CrystalBase::get_A_at_scan_point)
      .def("get_B_at_scan_point", &CrystalBase::get_B_at_scan_point)
      .def("get_U_at_scan_point", &CrystalBase::get_U_at_scan_point)
      .def("get_unit_cell_at_scan_point", &CrystalBase::get_unit_cell_at_scan_point)
      .def("reset_scan_points", &CrystalBase::reset_scan_points)
      .def("change_basis", &CrystalBase::change_basis)
      .def("update", &CrystalBase::update)
      .def("rotate_around_origin", &CrystalBase::rotate_around_origin, (
            arg("axis"),
            arg("angle"),
            arg("deg")=true))
      .def("is_similar_to", &CrystalBase::is_similar_to, (
            arg("other"),
            arg("angle_tolerance")=0.01,
            arg("uc_rel_length_tolerance")=0.01,
            arg("uc_abs_angle_tolerance")=1.0))
      .def("get_B_covariance", &CrystalBase::get_B_covariance)
      .def("set_B_covariance", &CrystalBase::set_B_covariance)
      .def("set_B_covariance", &Crystal_set_B_covariance_from_tuple)
      .def("get_cell_parameter_sd", &CrystalBase::get_cell_parameter_sd)
      .def("get_cell_volume_sd", &CrystalBase::get_cell_volume_sd)
      .def("reset_unit_cell_errors", &CrystalBase::reset_unit_cell_errors)
      .def("__eq__", &CrystalBase::operator==)
      .def("__ne__", &CrystalBase::operator!=);

    class_ <Crystal, bases <CrystalBase> > ("Crystal", no_init)
      .def(init<const Crystal&>())
      .def("__init__",
          make_constructor(
          &make_crystal_default,
          default_call_policies(), (
            arg("real_space_a"),
            arg("real_space_b"),
            arg("real_space_c"),
            arg("space_group"))))
      .def("__init__",
          make_constructor(
          &make_crystal_with_symbol,
          default_call_policies(), (
            arg("real_space_a"),
            arg("real_space_b"),
            arg("real_space_c"),
            arg("space_group_symbol"))))
      .def_pickle(CrystalPickleSuite());

    class_ <MosaicCrystalKabsch2010, bases <CrystalBase> > ("MosaicCrystalKabsch2010", no_init)
      .def(init<const MosaicCrystalKabsch2010&>())
      .def(init<const Crystal&>())
      .def("__init__",
          make_constructor(
          &make_kabsch2010_mosaic_crystal_default,
          default_call_policies(), (
            arg("real_space_a"),
            arg("real_space_b"),
            arg("real_space_c"),
            arg("space_group"))))
      .def("__init__",
          make_constructor(
          &make_kabsch2010_mosaic_crystal_with_symbol,
          default_call_policies(), (
            arg("real_space_a"),
            arg("real_space_b"),
            arg("real_space_c"),
            arg("space_group_symbol"))))
      .def("is_similar_to", &MosaicCrystalKabsch2010::is_similar_to, (
            arg("other"),
            arg("angle_tolerance")=0.01,
            arg("uc_rel_length_tolerance")=0.01,
            arg("uc_abs_angle_tolerance")=1.0,
            arg("mosaicity_tolerance")=0.8))
      .def("get_mosaicity", &MosaicCrystalKabsch2010::get_mosaicity, (
            arg("deg")=true))
      .def("set_mosaicity", &MosaicCrystalKabsch2010::set_mosaicity, (
            arg("mosaicity"),
            arg("deg")=true))
      .def_pickle(MosaicCrystalKabsch2010PickleSuite());

    class_ <MosaicCrystalSauter2014, bases <CrystalBase> > ("MosaicCrystalSauter2014", no_init)
      .def(init<const MosaicCrystalSauter2014&>())
      .def(init<const Crystal&>())
      .def("__init__",
          make_constructor(
          &make_sauter2014_mosaic_crystal_default,
          default_call_policies(), (
            arg("real_space_a"),
            arg("real_space_b"),
            arg("real_space_c"),
            arg("space_group"))))
      .def("__init__",
          make_constructor(
          &make_sauter2014_mosaic_crystal_with_symbol,
          default_call_policies(), (
            arg("real_space_a"),
            arg("real_space_b"),
            arg("real_space_c"),
            arg("space_group_symbol"))))
      .def("is_similar_to", &MosaicCrystalSauter2014::is_similar_to, (
            arg("other"),
            arg("angle_tolerance")=0.01,
            arg("uc_rel_length_tolerance")=0.01,
            arg("uc_abs_angle_tolerance")=1.0,
            arg("half_mosaicity_tolerance")=0.4,
            arg("domain_size_tolerance")=1.0))
      .def("get_half_mosaicity_deg", &MosaicCrystalSauter2014::get_half_mosaicity_deg)
      .def("set_half_mosaicity_deg", &MosaicCrystalSauter2014::set_half_mosaicity_deg, (
            arg("half_mosaicity_deg")))
      .def("get_domain_size_ang", &MosaicCrystalSauter2014::get_domain_size_ang)
      .def("set_domain_size_ang", &MosaicCrystalSauter2014::set_domain_size_ang, (
            arg("domain_size_ang")))
      .def_pickle(MosaicCrystalSauter2014PickleSuite());

    register_ptr_to_python<boost::shared_ptr<CrystalBase> >();
    register_ptr_to_python<boost::shared_ptr<Crystal> >();
    register_ptr_to_python<boost::shared_ptr<MosaicCrystalKabsch2010> >();
    register_ptr_to_python<boost::shared_ptr<MosaicCrystalSauter2014> >();
  }

}}} // namespace = dxtbx::model::boost_python
