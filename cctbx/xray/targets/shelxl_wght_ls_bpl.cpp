/*
http://tapenade.inria.fr:8080/tapenade/paste.jsp
  Paste shelxl_wght_ls.f
    top: kwt
    dependent: t
    independent: ic
    Reverse Mode
  unzip TapenadeResults.zip default_b-all.f
  scitbx.apply_tapenade_hints default_b-all.f > kwt_b.f
  diff -u kwt_b.f kwt_b_edited.f > kwt_b_patch
  cat shelxl_wght_ls.f kwt_b_edited.f | grep -v '^C'
  Paste output
    top: kwt_b
    dependent: t icb
    independent: ic
    Tangent Multidirectional Mode
  unzip TapenadeResults.zip default_dv-all.f
  scitbx.apply_tapenade_hints default_dv-all.f > kwt_b_dv.f
  fable.cout shelxl_wght_ls.f kwt_b_dv.f --inline-all --no-fem-do-safe --namespace=cctbx::xray::targets > shelxl_wght_ls.hpp
*/

#include <boost/python.hpp>
#include <cctbx/xray/targets/shelxl_wght_ls.hpp>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/mat_grid.h>

namespace cctbx { namespace xray { namespace boost_python {

boost::python::tuple
kwt_b_dv_wrapper(
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
  double tb = 1;
  af::shared<double> icb(nh);
  af::versa<double, af::flex_grid<> > icd(af::flex_grid<>(nh, nh));
  for(unsigned ih=0;ih<nh;ih++) icd(ih,ih) = 1;
  af::versa<double, af::flex_grid<> > icbd(af::flex_grid<>(nh, nh));
  targets::kwt_b_dv(
    0,
    tb,
    nh,
    f_obs.front(),
    i_obs.front(),
    i_sig.front(),
    ic.front(),
    icd.front(),
    icb.front(),
    icbd.front(),
    wa,
    wb,
    nh);
  return boost::python::make_tuple(icb, icbd);
}

void
wrap_targets_shelxl_wght_ls()
{
  using namespace boost::python;
  def("targets_shelxl_wght_ls_kwt_b_dv", kwt_b_dv_wrapper, (
    arg("f_obs"), arg("i_obs"), arg("i_sig"),
    arg("ic"),
    arg("wa"), arg("wb")));
}

}}}
