#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <vector>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace boost_python {

namespace {

  struct add_pair
  {
#if defined(__MACH__) && defined(__APPLE_CC__) && __APPLE_CC__ <= 1640
      bool dummy_;
#endif
    double tolerance_;
    af::shared<scitbx::vec3<double> > ref_model2_sites_;
    af::versa<bool, af::c_grid<2> > pair_flags_;
    std::vector<sgtbx::sym_equiv_sites<> > equiv1_;
    sgtbx::rt_mx eucl_symop_;
    int new_pair_1_;
    int new_pair_2_;

    add_pair(double tolerance,
             uctbx::unit_cell const& unit_cell,
             sgtbx::space_group const& space_group,
             double min_distance_sym_equiv,
             af::const_ref<scitbx::vec3<double> > const& ref_model1_sites,
             af::shared<scitbx::vec3<double> > const& ref_model2_sites)
    :
      tolerance_(tolerance),
      ref_model2_sites_(ref_model2_sites),
      pair_flags_(af::versa<bool, af::c_grid<2> >(
        af::c_grid<2>(ref_model1_sites.size(),
                      ref_model2_sites.size())))
    {
      for(std::size_t i=0;i<ref_model1_sites.size();i++) {
        sgtbx::site_symmetry site_symmetry(
          unit_cell, space_group,
          ref_model1_sites[i], min_distance_sym_equiv);
        equiv1_.push_back(sgtbx::sym_equiv_sites<>(site_symmetry));
      }
    }

    sgtbx::sym_equiv_sites<> const&
    equiv1(std::size_t i) const { return equiv1_[i]; }

    void
    next_pivot(af::tiny<bool, 3> const& continuous_shift_flags,
               sgtbx::rt_mx const& eucl_symop,
               scitbx::vec3<double> const& adjusted_shift,
               af::const_ref<int> const& singles1,
               af::const_ref<int> const& singles2)
    {
      CCTBX_ASSERT(ref_model2_sites_.size() == pair_flags_.accessor()[1]);
      pair_flags_.fill(false);
      eucl_symop_ = eucl_symop;
      CCTBX_ASSERT(equiv1_.size() > 0);
      bool no_continuous_shifts = continuous_shift_flags.all_eq(false);
      uctbx::unit_cell const& unit_cell = equiv1_[0].unit_cell();
      af::ref<bool, af::c_grid<2> > pair_flags_ref = pair_flags_.ref();
      for (const int* is2 = singles2.begin(); is2 != singles2.end(); is2++) {
        scitbx::vec3<double>
          c2 = eucl_symop_ * ref_model2_sites_[*is2] + adjusted_shift;
        for (const int* is1 = singles1.begin(); is1 != singles1.end(); is1++) {
          sgtbx::min_sym_equiv_distance_info<> dist_info(equiv1_[*is1], c2);
          double dist = dist_info.dist();
          if (no_continuous_shifts) {
            if (dist < tolerance_) {
              pair_flags_ref(*is1,*is2) = true;
            }
          }
          else if (dist < 4 * tolerance_) {
            // ensure that this pair can be matched within 1 * tolerance.
            // (not entirely sure that this is safe under all circumstances.)
            fractional<> const& diff = dist_info.diff();
            fractional<> diff_allowed;
            for (std::size_t j=0;j<3;j++) {
              diff_allowed[j] = (continuous_shift_flags[j] ? diff[j] : 0);
            }
            double dist_allowed = unit_cell.length(diff_allowed);
            if (dist_allowed < tolerance_) {
              pair_flags_ref(*is1,*is2) = true;
            }
          }
        }
      }
    }

    bool
    next_pair(scitbx::vec3<double> const& adjusted_shift,
              af::const_ref<int> const& singles1,
              af::const_ref<int> const& singles2)
    {
      double shortest_dist = 2 * tolerance_;
      bool have_new_pair = false;
      new_pair_1_ = 0;
      new_pair_2_ = 0;
      af::const_ref<bool, af::c_grid<2> >
        pair_flags_ref = pair_flags_.const_ref();
      for (const int* is2 = singles2.begin(); is2 != singles2.end(); is2++) {
        scitbx::vec3<double>
          c2 = eucl_symop_ * ref_model2_sites_[*is2] + adjusted_shift;
        for (const int* is1 = singles1.begin(); is1 != singles1.end(); is1++) {
          if (pair_flags_ref(*is1,*is2)) {
            sgtbx::min_sym_equiv_distance_info<> dist_info(equiv1_[*is1], c2);
            double dist = dist_info.dist();
            if (dist < shortest_dist) {
              shortest_dist = dist;
              have_new_pair = true;
              new_pair_1_ = *is1;
              new_pair_2_ = *is2;
            }
          }
        }
      }
      return have_new_pair;
    }

    int
    new_pair_1() const { return new_pair_1_; }

    int
    new_pair_2() const { return new_pair_2_; }
  };

  void init_module()
  {
    using namespace boost::python;

    typedef return_internal_reference<> rir;

    class_<add_pair>("add_pair", no_init)
      .def(init<double,
                uctbx::unit_cell const&,
                sgtbx::space_group const&,
                double,
                af::const_ref<scitbx::vec3<double> > const&,
                af::shared<scitbx::vec3<double> > const&>())
      .def("equiv1", &add_pair::equiv1, rir())
      .def("next_pivot", &add_pair::next_pivot)
      .def("next_pair", &add_pair::next_pair)
      .def("new_pair_1", &add_pair::new_pair_1)
      .def("new_pair_2", &add_pair::new_pair_2)
    ;
  }

} // namespace <anonymous>
}} // namespace cctbx::boost_python

BOOST_PYTHON_MODULE(cctbx_emma_ext)
{
  cctbx::boost_python::init_module();
}
