#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <iotbx/mtz/batch.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace iotbx { namespace mtz {
namespace {

#define IOTBX_MTZ_BATCH_BPL_GET_SET(name) \
        .def(#name, &w_t::name) \
        .def("set_" #name, &w_t::set_##name, (arg_("value")), return_self<>())

#define IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(name) \
        .def(#name, &w_t::name) \
        .def("set_" #name, &w_t::set_##name, (arg_("values")), return_self<>())

  struct batch_wrappers
  {
    typedef batch w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("batch", no_init)
        .def(init<object const&, int>((arg_("mtz_object"), arg_("i_batch"))))
        .def("mtz_object", &w_t::mtz_object)
        .def("i_batch", &w_t::i_batch)
        IOTBX_MTZ_BATCH_BPL_GET_SET(num)
        IOTBX_MTZ_BATCH_BPL_GET_SET(title)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(gonlab)
        IOTBX_MTZ_BATCH_BPL_GET_SET(iortyp)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(lbcell)
        IOTBX_MTZ_BATCH_BPL_GET_SET(misflg)
        IOTBX_MTZ_BATCH_BPL_GET_SET(jumpax)
        IOTBX_MTZ_BATCH_BPL_GET_SET(ncryst)
        IOTBX_MTZ_BATCH_BPL_GET_SET(lcrflg)
        IOTBX_MTZ_BATCH_BPL_GET_SET(ldtype)
        IOTBX_MTZ_BATCH_BPL_GET_SET(jsaxs)
        IOTBX_MTZ_BATCH_BPL_GET_SET(nbscal)
        IOTBX_MTZ_BATCH_BPL_GET_SET(ngonax)
        IOTBX_MTZ_BATCH_BPL_GET_SET(lbmflg)
        IOTBX_MTZ_BATCH_BPL_GET_SET(ndet)
        IOTBX_MTZ_BATCH_BPL_GET_SET(nbsetid)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(cell)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(umat)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(phixyz)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(crydat)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(datum)
        IOTBX_MTZ_BATCH_BPL_GET_SET(phistt)
        IOTBX_MTZ_BATCH_BPL_GET_SET(phiend)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(scanax)
        IOTBX_MTZ_BATCH_BPL_GET_SET(time1)
        IOTBX_MTZ_BATCH_BPL_GET_SET(time2)
        IOTBX_MTZ_BATCH_BPL_GET_SET(bscale)
        IOTBX_MTZ_BATCH_BPL_GET_SET(bbfac)
        IOTBX_MTZ_BATCH_BPL_GET_SET(sdbscale)
        IOTBX_MTZ_BATCH_BPL_GET_SET(sdbfac)
        IOTBX_MTZ_BATCH_BPL_GET_SET(phirange)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(e1)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(e2)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(e3)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(source)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(so)
        IOTBX_MTZ_BATCH_BPL_GET_SET(alambd)
        IOTBX_MTZ_BATCH_BPL_GET_SET(delamb)
        IOTBX_MTZ_BATCH_BPL_GET_SET(delcor)
        IOTBX_MTZ_BATCH_BPL_GET_SET(divhd)
        IOTBX_MTZ_BATCH_BPL_GET_SET(divvd)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(dx)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(theta)
        IOTBX_MTZ_BATCH_BPL_GET_SET_ARRAY(detlm)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_batch");
      }
    }
  };

  void
  wrap_all()
  {
    batch_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_batch() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
