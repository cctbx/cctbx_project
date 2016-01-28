#include <mmtbx/bulk_solvent/bulk_solvent.h>
#include <cctbx/xray/targets.h>

namespace mmtbx { namespace bulk_solvent {


double fu_star(sym_mat3<double> const& u_star,
               cctbx::miller::index<> const& mi)
{
    double arg = -0.25 * (u_star[0]*mi[0]*mi[0] +
                          u_star[1]*mi[1]*mi[1] +
                          u_star[2]*mi[2]*mi[2] +
                       2.*u_star[3]*mi[0]*mi[1] +
                       2.*u_star[4]*mi[0]*mi[2] +
                       2.*u_star[5]*mi[1]*mi[2]);
    if(arg > 40.0) arg=40.0; // to avoid overflow problem
    return std::exp(arg);
}

af::shared<double> fb_cart(sym_mat3<double> const& b_cart,
                            af::const_ref<cctbx::miller::index<> > const& hkl,
                            cctbx::uctbx::unit_cell const& uc)
{
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (b_cart).tensor_transform(a);
    af::shared<double> fu_mem(hkl.size(), af::init_functor_null<double>());
    double* fu = fu_mem.begin();
    for(std::size_t i=0; i < hkl.size(); i++) {
      fu[i] = fu_star(u_star, hkl[i]);
    }
    return fu_mem;
}

}} // namespace mmtbx::bulk_solvent
