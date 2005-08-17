#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/reciprocal_space_reference_asu.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace sgtbx {

  phase_info::phase_info(
    sgtbx::space_group const& space_group,
    miller::index<> const& miller_index,
    bool no_test_sys_absent)
  :
    ht_(-1),
    t_den_(space_group.t_den()),
    sys_abs_was_tested_(!no_test_sys_absent)
  {
    miller::index<> const& h = miller_index;
    // systematically absent reflection: if hr == h and ht != 0 mod 1
    // restricted phase: if hr == -h: phi(h) = pi*ht + n*pi
    if (no_test_sys_absent) {
      // Fast determination of phase restriction without considering
      // conditions for systematically absent reflections.
      if (space_group.is_centric()) {
        ht_ = ht_mod_1(h, space_group.inv_t());
      }
      else {
        for(std::size_t i=0;i<space_group.n_smx();i++) {
          if (h * space_group.smx(i).r() == -h) {
            ht_ = ht_mod_1(h, space_group.smx(i).t());
            break;
          }
        }
      }
      return;
    }
    // Simulatenous determination of phase restriction and evaluation
    // conditions for systematically absent reflections.
    for(std::size_t i_smx=0;i_smx<space_group.n_smx();i_smx++) {
      rot_mx const& r = space_group.smx(i_smx).r();
      tr_vec const& t = space_group.smx(i_smx).t();
      tr_vec ts(0);
      tr_vec tr(0);
      miller::index<> hr = h * r;
      if      (h == hr) {
        ts = t;
        if (space_group.is_centric()) tr = space_group.inv_t() - t;
      }
      else if (h == -hr) {
        tr = t;
        if (space_group.is_centric()) ts = space_group.inv_t() - t;
      }
      if (ts.is_valid()) {
        for(std::size_t i_ltr=0;i_ltr<space_group.n_ltr();i_ltr++) {
          if ((h * (ts + space_group.ltr(i_ltr))) % ts.den() != 0) {
            ht_ = -2;
            return;
          }
        }
      }
      if (tr.is_valid()) {
        for(std::size_t i_ltr=0;i_ltr<space_group.n_ltr();i_ltr++) {
          int ht = ht_mod_1(h, tr + space_group.ltr(i_ltr));
          if      (ht_ < 0) ht_ = ht;
          else if (ht_ != ht) {
            ht_ = -2;
            return;
          }
        }
      }
    }
  }

  bool
  phase_info::is_valid_phase(double phi, bool deg, double tolerance) const
  {
    if (!is_centric()) return true;
    double pi_u = pi_unit(deg);
    double delta = std::fmod(phi - ht_angle(deg), pi_u);
    if (delta >  tolerance) delta -= pi_u;
    if (delta < -tolerance) delta += pi_u;
    if (delta <= tolerance) return true;
    return false;
  }

  double
  phase_info::nearest_valid_phase(double phi, bool deg) const
  {
    if (!is_centric()) return phi;
    double pi_u = pi_unit(deg);
    double phi_restr = ht_angle(deg);
    double delta = math::fmod_short(phi - phi_restr, 2*pi_u);
    if (delta <= -pi_u/2 || delta > pi_u/2) return phi_restr + pi_u;
    return phi_restr;
  }

  af::shared<bool>
  space_group::is_sys_absent(
    af::const_ref<miller::index<> > const& miller_indices) const
  {
    af::shared<bool> result(
      miller_indices.size(), af::init_functor_null<bool>());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      result[i] = is_sys_absent(miller_indices[i]);
    }
    return result;
  }

  bool space_group::is_centric(miller::index<> const& miller_index) const
  {
    if (is_centric()) return true;
    if (n_smx() > 1) {
      miller::index<> minus_miller_index = -miller_index;
      if (miller_index * smx_[1].r() == minus_miller_index) return true;
      for(std::size_t i=2;i<n_smx();i++) {
        if (miller_index * smx_[i].r() == minus_miller_index) return true;
      }
    }
    return false;
  }

  af::shared<bool>
  space_group::is_centric(
    af::const_ref<miller::index<> > const& miller_indices) const
  {
    af::shared<bool> result(
      miller_indices.size(), af::init_functor_null<bool>());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      result[i] = is_centric(miller_indices[i]);
    }
    return result;
  }

  int space_group::epsilon(miller::index<> const& miller_index) const
  {
    int result = 1;
    for(std::size_t i=1;i<n_smx();i++) {
      miller::index<> hr = miller_index * smx_[i].r();
      if (hr == miller_index || (is_centric() && hr == -miller_index))
        result++;
    }
    CCTBX_ASSERT(n_smx() % result == 0);
    return result;
  }

  af::shared<int>
  space_group::epsilon(
    af::const_ref<miller::index<> > const& miller_indices) const
  {
    af::shared<int> result(
      miller_indices.size(), af::init_functor_null<int>());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      result[i] = epsilon(miller_indices[i]);
    }
    return result;
  }

  int
  space_group
  ::multiplicity(
    miller::index<> const& miller_index,
    bool anomalous_flag) const
  {
    if (miller_index.is_zero()) return 1;
    int centro = (is_centric() || !anomalous_flag);
    int m = 1;
    int r = 0;
    for(std::size_t i=1;i<n_smx();i++) {
      miller::index<> hr = miller_index * smx_[i].r();
      if      (hr == miller_index) m++;
      else if (hr == -miller_index) r++;
    }
    CCTBX_ASSERT(n_smx() % m == 0 && (r == 0 || r == m));
    m = n_smx() / m;
    if (centro && r == 0) m *= 2;
    return m;
  }

  af::shared<int>
  space_group
  ::multiplicity(
    af::const_ref<miller::index<> > const& miller_indices,
    bool anomalous_flag) const
  {
    af::shared<int> result(
      miller_indices.size(), af::init_functor_null<int>());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      result[i] = multiplicity(miller_indices[i], anomalous_flag);
    }
    return result;
  }

}} // namespace cctbx::sgtbx
