#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <boost/python/args.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <spotfinder/core_toolbox/distl.h>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct to_string_ir : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    to_string_ir()
    {
      int version = 1;
      *this << version;
    }

    to_string_ir& operator<<(Distl::icering const& val)
    {
      *this << val.lowerr2
            << val.upperr2
            << val.lowerresol
            << val.upperresol
            << val.strength;
      return *this;
    }
  };

  struct from_string_ir : pickle_double_buffered::from_string
  {
    from_string_ir(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr)
    {
          *this >> version;
    }

    using pickle_double_buffered::from_string::operator>>;

    from_string_ir& operator>>(Distl::icering& val)
    {
      Distl::point lpeak;

      *this >> val.lowerr2
            >> val.upperr2
            >> val.lowerresol
            >> val.upperresol
            >> val.strength;
      return *this;
    }

    int version;
  };

  struct to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    to_string()
    {
      int version = 2; //Current version for flex_distl_spot
      //printf("Pickle flex_distl_spot; version: %d\n", version);
      *this << version;
    }

    to_string& operator<<(scitbx::vec2<double> const& vec)
    {
      //special time saver for 2D profiles:  ignore z-direction
      *this << (float)vec[0] << (float)vec[1];
      return *this;
    }

    to_string& operator<<(spotfinder::distltbx::w_spot const& val)
    {
      *this << val.max_pxl_x()
            << val.max_pxl_y()
            << val.peakintensity
            << val.bodypixels.size();

      for (int i=0; i<val.bodypixels.size(); ++i){
      *this << val.bodypixels[i].x
            << val.bodypixels[i].y;
      }

      *this << val.total_mass;
      to_string::operator<<(val.model_m->center_of_mass());
      *this << (float)val.model_m->eigenvalue(0)
            << (float)val.model_m->eigenvalue(1);

      to_string::operator<<(val.model_m->eigenvector(0));
      to_string::operator<<(val.model_m->eigenvector(1));

      return *this;
    }
  };

  struct from_string : pickle_double_buffered::from_string
  {
    from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr)
    {
      *this >> version;
      //printf("Unpickle flex_distl_spot; version: %d\n", version);
    }

    from_string& operator>>(scitbx::vec2<double>& vec)
    {
      //special time saver for 2D profiles:  ignore z-direction
      float a,b;
      *this >> a >> b;
      vec[0] = (double)a; vec[1] = (double)b;
      return *this;
    }

    using pickle_double_buffered::from_string::operator>>;

    from_string& operator>>(spotfinder::distltbx::w_spot& val)
    {
      int iptx,ipty,bodylen;
      Distl::point lpeak;

      *this >> iptx
            >> ipty;
      lpeak = Distl::point(iptx,ipty);
      *this >> val.peakintensity
            >> bodylen;

      val.bodypixels = scitbx::af::shared<Distl::point>();
      for (int i=0; i<bodylen; ++i){
      *this >> iptx
            >> ipty;
        val.bodypixels.push_back(Distl::point(iptx,ipty));
      }

      val.setstate(lpeak);
      if (version == 2) {
        double mass;
        *this >> mass;
        scitbx::vec2<double> com,value,axis0,axis1;
        float a,b;
        from_string::operator>>(com);
        *this >> a >> b;
        value[0] = (double)a; value[1] = (double)b;
        from_string::operator>>(axis0);
        from_string::operator>>(axis1);
        val.ss_setstate(com,value,axis0,axis1,mass);
      }
      return *this;
    }

    int version;
  };

}}}} // namespace scitbx::af::boost_python::<anonymous>


namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_icering()
  {
    flex_wrapper<Distl::icering,
                 boost::python::return_internal_reference<>
                 >::plain("distl_icering")
    .def_pickle(flex_pickle_double_buffered<
                 Distl::icering, to_string_ir, from_string_ir>())
    ;
  }

  af::shared<double>
  ctr_mass_distances_from_direct_beam(
    af::const_ref<spotfinder::distltbx::w_spot> const& spots,
    scitbx::vec2<double> detector_size,
    scitbx::vec2<int> detector_pixels,
    scitbx::vec2<double> xy_beam)
  {
    af::shared<double> result(
      spots.size(),
      af::init_functor_null<double>());
    scitbx::vec2<double> sop;
    for(std::size_t i=0;i<2;i++) {
      sop[i] = detector_size[i] / detector_pixels[i];
    }
    double const* b = xy_beam.begin();
    for(std::size_t i=0;i<spots.size();i++) {
      spotfinder::distltbx::w_spot const& spot = spots[i];
      double dx = spot.ctr_mass_x() * sop[0] - b[0];
      double dy = spot.ctr_mass_y() * sop[1] - b[1];
      result.push_back(std::sqrt(dx*dx + dy*dy));
    }
    return result;
  }

  void wrap_flex_w_spot()
  {
    using boost::python::arg;
    flex_wrapper<spotfinder::distltbx::w_spot,
                 boost::python::return_internal_reference<>
                 >::plain("distl_spot")
    .def_pickle(flex_pickle_double_buffered<
                 spotfinder::distltbx::w_spot, to_string, from_string>())
    .def("ctr_mass_distances_from_direct_beam",
      ctr_mass_distances_from_direct_beam, (
        arg("detector_size"),
        arg("detector_pixels"),
        arg("xy_beam")))
    ;
  }

  void wrap_flex_point()
  {
    flex_wrapper<Distl::point>::plain("distl_point");
  }

}}} // namespace scitbx::af::boost_python
