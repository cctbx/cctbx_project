/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 Oct: Redesign: AsymIndex (rwgk)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/reciprocal_space_reference_asu.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace sgtbx {

  phase_info::phase_info(
    space_group const& sg,
    miller::index<> const& h,
    bool no_test_sys_absent)
    : ht_(-1), t_den_(sg.t_den()), sys_abs_was_tested_(!no_test_sys_absent)
  {
    // systematically absent reflection: if hr == h and ht != 0 mod 1
    // restricted phase: if hr == -h: phi(h) = pi*ht + n*pi
    if (no_test_sys_absent) {
      // Fast determination of phase restriction without considering
      // conditions for systematically absent reflections.
      if (sg.is_centric()) {
        ht_ = ht_mod_1(h, sg.inv_t());
      }
      else {
        for(std::size_t i=0;i<sg.n_smx();i++) {
          if (h * sg.smx(i).r() == -h) {
            ht_ = ht_mod_1(h, sg.smx(i).t());
            break;
          }
        }
      }
      return;
    }
    // Simulatenous determination of phase restriction and evaluation
    // conditions for systematically absent reflections.
    for(std::size_t i_smx=0;i_smx<sg.n_smx();i_smx++) {
      rot_mx const& r = sg.smx(i_smx).r();
      tr_vec const& t = sg.smx(i_smx).t();
      tr_vec ts(0);
      tr_vec tr(0);
      miller::index<> hr = h * r;
      if      (h == hr) {
        ts = t;
        if (sg.is_centric()) tr = sg.inv_t() - t;
      }
      else if (h == -hr) {
        tr = t;
        if (sg.is_centric()) ts = sg.inv_t() - t;
      }
      if (ts.is_valid()) {
        for(std::size_t i_ltr=0;i_ltr<sg.n_ltr();i_ltr++) {
          if ((h * (ts + sg.ltr(i_ltr))) % ts.den() != 0) {
            ht_ = -2;
            return;
          }
        }
      }
      if (tr.is_valid()) {
        for(std::size_t i_ltr=0;i_ltr<sg.n_ltr();i_ltr++) {
          int ht = ht_mod_1(h, tr + sg.ltr(i_ltr));
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
  space_group::is_sys_absent(af::const_ref<miller::index<> > const& h) const
  {
    af::shared<bool> result(h.size(), af::init_functor_null<bool>());
    for(std::size_t i=0;i<h.size();i++) result[i] = is_sys_absent(h[i]);
    return result;
  }

  bool space_group::is_centric(miller::index<> const& h) const
  {
    if (is_centric()) return true;
    for(std::size_t i=1;i<n_smx();i++) {
      if (h * smx_[i].r() == -h) return true;
    }
    return false;
  }

  af::shared<bool>
  space_group::is_centric(af::const_ref<miller::index<> > const& h) const
  {
    af::shared<bool> result(h.size(), af::init_functor_null<bool>());
    for(std::size_t i=0;i<h.size();i++) result[i] = is_centric(h[i]);
    return result;
  }

  int space_group::epsilon(miller::index<> const& h) const
  {
    int result = 1;
    for(std::size_t i=1;i<n_smx();i++) {
      miller::index<> hr = h * smx_[i].r();
      if (hr == h || (is_centric() && hr == -h))
        result++;
    }
    CCTBX_ASSERT(n_smx() % result == 0);
    return result;
  }

  af::shared<int>
  space_group::epsilon(af::const_ref<miller::index<> > const& h) const
  {
    af::shared<int> result(h.size(), af::init_functor_null<int>());
    for(std::size_t i=0;i<h.size();i++) result[i] = epsilon(h[i]);
    return result;
  }

  int
  space_group
  ::multiplicity(miller::index<> const& h, bool anomalous_flag) const
  {
    if (h.is_zero()) return 1;
    int centro = (is_centric() || !anomalous_flag);
    int m = 1;
    int r = 0;
    for(std::size_t i=1;i<n_smx();i++) {
      miller::index<> hr = h * smx_[i].r();
      if      (hr == h) m++;
      else if (hr == -h) r++;
    }
    CCTBX_ASSERT(n_smx() % m == 0 && (r == 0 || r == m));
    m = n_smx() / m;
    if (centro && r == 0) m *= 2;
    return m;
  }

  af::shared<int>
  space_group::multiplicity(af::const_ref<miller::index<> > const& h,
                            bool anomalous_flag) const
  {
    af::shared<int> result(h.size(), af::init_functor_null<int>());
    for(std::size_t i=0;i<h.size();i++) {
      result[i] = multiplicity(h[i], anomalous_flag);
    }
    return result;
  }

}} // namespace cctbx::sgtbx
