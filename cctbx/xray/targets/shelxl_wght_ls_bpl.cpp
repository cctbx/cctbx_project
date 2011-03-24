/*
http://tapenade.inria.fr:8080/tapenade/paste.jsp
  Paste shelxl_wght_ls.f
    top: kwt
    dependent: t
    independent: ic
    Tangent Multidirectional Mode
  Download differentiated file
  unzip TapenadeResults.zip default_dv-all.f
  scitbx.apply_tapenade_hints default_dv-all.f --no-comments
  Paste output
    top: kwt_dv
    dependent: t td
    independent: ic
    Tangent Multidirectional Mode
  scitbx.apply_tapenade_hints default_dv-all.f --no-comments > kwt_dv_dv.f
  fable.cout kwt_dv_dv.f --inline-all --no-fem-do-safe --namespace=cctbx::xray::targets > shelxl_wght_ls.hpp
*/

#include <boost/python.hpp>
#include <cctbx/xray/targets/shelxl_wght_ls.hpp>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/mat_grid.h>

namespace cctbx { namespace xray { namespace boost_python {

boost::python::tuple
kwt_dv_dv_wrapper(
  af::const_ref<double> const& f_obs,
  af::const_ref<double> const& i_obs,
  af::const_ref<double> const& i_sig,
  af::const_ref<double> const& ic,
  double wa,
  double wb)
{
  TBXX_ASSERT(i_obs.size() == f_obs.size());
  TBXX_ASSERT(i_sig.size() == f_obs.size());
  TBXX_ASSERT(ic.size() == f_obs.size());
  std::size_t nh = static_cast<int>(f_obs.size());
  double t;
  af::shared<double> td0(nh);
  af::shared<double> td(nh);
  af::versa<double, af::flex_grid<> > tdd(af::flex_grid<>(nh, nh));
  af::versa<double, af::flex_grid<> > icd(af::flex_grid<>(nh, nh));
  for(unsigned ih=0;ih<nh;ih++) icd(ih,ih) = 1;
  targets::kwt_dv_dv(
    t,
    td0.front(),
    td.front(),
    tdd.front(),
    nh,
    f_obs.front(),
    i_obs.front(),
    i_sig.front(),
    ic.front(),
    icd.front(),
    icd.front(),
    wa,
    wb,
    nh,
    nh);
  return boost::python::make_tuple(t, td0, td, tdd);
}

void
wrap_targets_shelxl_wght_ls()
{
  using namespace boost::python;
  def("targets_shelxl_wght_ls_kwt_dv_dv", kwt_dv_dv_wrapper, (
    arg("f_obs"), arg("i_obs"), arg("i_sig"),
    arg("ic"),
    arg("wa"), arg("wb")));
}

}}}
