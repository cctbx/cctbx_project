#include <cctbx/miller/asu.h>
#include <set>

namespace cctbx { namespace miller {

  asym_index::asym_index(
    sgtbx::space_group const& sg,
    sgtbx::reciprocal_space::asu const& asu,
    index<> const& h)
  {
    t_den_ = sg.t_den();
    friedel_flag_ = false;
    for(std::size_t i_inv=0;i_inv<sg.f_inv();i_inv++) {
      for(std::size_t i_smx=0;i_smx<sg.n_smx();i_smx++) {
        sgtbx::rt_mx s = sg(0, i_inv, i_smx);
        hr_ = h * s.r();
        if (asu.is_inside(hr_)) {
          ht_ = sgtbx::ht_mod_1(h, s.t());
          return;
        }
      }
    }
    CCTBX_ASSERT(!sg.is_centric());
    for(std::size_t i_smx=0;i_smx<sg.n_smx();i_smx++) {
      sgtbx::rt_mx s = sg(0, 0, i_smx);
      hr_ = h * s.r();
      if (asu.is_inside(-hr_)) {
        ht_ = sgtbx::ht_mod_1(h, s.t());
        friedel_flag_ = true;
        return;
      }
    }
    throw CCTBX_INTERNAL_ERROR();
  }

  asym_index::asym_index(sym_equiv_indices const& h_eq)
  {
    t_den_ = h_eq.indices()[0].t_den();
    std::size_t i_selected = 0;
    index<> selected_h = h_eq.indices()[0].hr();
    friedel_flag_ = false;
    for(std::size_t i_indices=0;i_indices<h_eq.indices().size();i_indices++) {
      index<> trial_h = h_eq.indices()[i_indices].hr();
      for(std::size_t i_mate=0;i_mate<h_eq.f_mates(false);i_mate++) {
        if (i_mate) trial_h = -trial_h;
        if (trial_h < selected_h) {
          i_selected = i_indices;
          selected_h = trial_h;
          friedel_flag_ = (i_mate != 0);
        }
      }
    }
    hr_ = h_eq.indices()[i_selected].hr();
    ht_ = h_eq.indices()[i_selected].ht();
  }

  asym_index::asym_index(
    sgtbx::space_group const& sg,
    index<> const& h)
  {
    *this = asym_index(sym_equiv_indices(sg, h));
  }

  void
  map_to_asu(
    sgtbx::space_group_type const& sg_type,
    bool anomalous_flag,
    af::ref<index<> > const& miller_indices)
  {
    sgtbx::reciprocal_space::asu asu(sg_type);
    sgtbx::space_group const& sg = sg_type.group();
    for(std::size_t i=0;i<miller_indices.size();i++) {
      asym_index ai(sg, asu, miller_indices[i]);
      index_table_layout_adaptor ila = ai.one_column(anomalous_flag);
      miller_indices[i] = ila.h();
    }
  }

  bool
  is_unique_set_under_symmetry(
    sgtbx::space_group_type const& space_group_type,
    bool anomalous_flag,
    af::const_ref<index<> > const& miller_indices)
  {
    std::set<index<>, fast_less_than<> > unique_set;
    sgtbx::reciprocal_space::asu asu(space_group_type);
    sgtbx::space_group const& sg = space_group_type.group();
    for(std::size_t i=0;i<miller_indices.size();i++) {
      asym_index ai(sg, asu, miller_indices[i]);
      index_table_layout_adaptor ila = ai.one_column(anomalous_flag);
      if (unique_set.find(ila.h()) != unique_set.end()) return false;
      unique_set.insert(ila.h());
    }
    return true;
  }

  af::shared<std::size_t>
  unique_under_symmetry_selection(
    sgtbx::space_group_type const& space_group_type,
    bool anomalous_flag,
    af::const_ref<index<> > const& miller_indices)
  {
    af::shared<std::size_t> result((af::reserve(miller_indices.size())));
    std::set<index<>, fast_less_than<> > unique_set;
    sgtbx::reciprocal_space::asu asu(space_group_type);
    sgtbx::space_group const& sg = space_group_type.group();
    for(std::size_t i=0;i<miller_indices.size();i++) {
      asym_index ai(sg, asu, miller_indices[i]);
      index_table_layout_adaptor ila = ai.one_column(anomalous_flag);
      if (unique_set.find(ila.h()) == unique_set.end()) {
        unique_set.insert(ila.h());
        result.push_back(i);
      }
    }
    return result;
  }

}} // namespace cctbx::miller
