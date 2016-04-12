#ifndef CCTBX_CRYSTAL_INCREMENTAL_PAIRS_H
#define CCTBX_CRYSTAL_INCREMENTAL_PAIRS_H

#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace crystal {

  template <typename FloatType=double, typename IntShiftType=int>
  class incremental_pairs
  {
    protected:
      FloatType distance_cutoff_;
      FloatType distance_cutoff_sq_;
      FloatType cubicle_epsilon_;
      typedef
        direct_space_asu::asu_mappings<FloatType, IntShiftType>
          asu_mappings_t;
      boost::shared_ptr<asu_mappings_t> asu_mappings_owner_;
      asu_mappings_t* asu_mappings_;
      typedef std::vector<direct_space_asu::asu_mapping_index>
        cubicle_content_t;
      scitbx::cubicles<cubicle_content_t, FloatType> cubicles_;
      typedef crystal::pair_asu_table<FloatType, IntShiftType>
        pair_asu_table_t;
      boost::shared_ptr<pair_asu_table_t> pair_asu_table_owner_;
      pair_asu_table_t* pair_asu_table_;

    public:
      //! Default constructor. Some data members are not initialized!
      incremental_pairs() {}

      incremental_pairs(
        sgtbx::space_group const& space_group,
        direct_space_asu::float_asu<FloatType> const& asu,
        FloatType const& distance_cutoff,
        FloatType const& asu_mappings_buffer_thickness=-1,
        FloatType const& cubicle_epsilon=-1)
      :
        distance_cutoff_(distance_cutoff),
        distance_cutoff_sq_(distance_cutoff*distance_cutoff),
        cubicle_epsilon_(cubicle_epsilon >= 0
          ? cubicle_epsilon
          : asu.is_inside_epsilon()),
        asu_mappings_owner_(new asu_mappings_t(
          space_group,
          asu,
          asu_mappings_buffer_thickness >= 0
            ? asu_mappings_buffer_thickness
            : distance_cutoff)),
        asu_mappings_(asu_mappings_owner_.get()),
        cubicles_(
          asu_mappings_->asu_buffer().box_min(/*cartesian*/ true),
          asu_mappings_->asu_buffer().box_span(/*cartesian*/ true),
          distance_cutoff_,
          cubicle_epsilon_),
        pair_asu_table_owner_(new pair_asu_table_t(asu_mappings_owner_)),
        pair_asu_table_(pair_asu_table_owner_.get()),
        min_distance_sym_equiv(0.5),
        assert_min_distance_sym_equiv(true)
      {
        CCTBX_ASSERT(distance_cutoff_ > 0);
        CCTBX_ASSERT(asu_mappings_->buffer_thickness() >= distance_cutoff_);
        init_space_group = space_group;
        init_asu = asu;
        init_distance_cutoff = distance_cutoff;
        init_asu_mappings_buffer_thickness = asu_mappings_buffer_thickness;
        init_cubicle_epsilon = cubicle_epsilon;
      }
      sgtbx::space_group init_space_group;
      direct_space_asu::float_asu<FloatType> init_asu;
      FloatType init_distance_cutoff;
      FloatType init_asu_mappings_buffer_thickness;
      FloatType init_cubicle_epsilon;

      boost::shared_ptr<
        direct_space_asu::asu_mappings<FloatType, IntShiftType> >
      asu_mappings() const { return asu_mappings_owner_; }

      boost::shared_ptr<
        crystal::pair_asu_table<FloatType, IntShiftType> >
      pair_asu_table() const { return pair_asu_table_owner_; }

      FloatType min_distance_sym_equiv;
      bool assert_min_distance_sym_equiv;

      void
      process_site_frac(
        fractional<FloatType> const& original_site,
        sgtbx::site_symmetry_ops const& site_symmetry_ops)
      {
        direct_space_asu::asu_mapping_index mi;
        mi.i_seq = asu_mappings_->mappings().size();
        asu_mappings_->process(original_site, site_symmetry_ops);
        pair_asu_table_->grow(1);
        sgtbx::rt_mx
          rt_mx_i_inverse = asu_mappings_->get_rt_mx(mi.i_seq, 0).inverse();
        af::const_ref<
          typename asu_mappings_t::array_of_mappings_for_one_site> const&
            mappings = asu_mappings_->mappings_const_ref();
        typename asu_mappings_t::array_of_mappings_for_one_site const&
          mappings_i = mappings.back();
        scitbx::vec3<unsigned> n_cub = cubicles_.ref.accessor();
        for(mi.i_sym=0; mi.i_sym<mappings_i.size(); mi.i_sym++) {
          scitbx::vec3<unsigned>
            i_cub = cubicles_.i_cubicle(mappings_i[mi.i_sym].mapped_site());
          scitbx::vec3<unsigned> j_cub_min;
          scitbx::vec3<unsigned> j_cub_max;
          j_cub_min[0] = (i_cub[0] == 0 ? 0 : i_cub[0]-1);
          j_cub_max[0] = (i_cub[0] == n_cub[0]-1 ? i_cub[0] : i_cub[0]+1);
          j_cub_min[1] = (i_cub[1] == 0 ? 0 : i_cub[1]-1);
          j_cub_max[1] = (i_cub[1] == n_cub[1]-1 ? i_cub[1] : i_cub[1]+1);
          j_cub_min[2] = (i_cub[2] == 0 ? 0 : i_cub[2]-1);
          j_cub_max[2] = (i_cub[2] == n_cub[2]-1 ? i_cub[2] : i_cub[2]+1);
          scitbx::vec3<unsigned> j_cub;
          for(j_cub[0]=j_cub_min[0];j_cub[0]<=j_cub_max[0];j_cub[0]++) {
          for(j_cub[1]=j_cub_min[1];j_cub[1]<=j_cub_max[1];j_cub[1]++) {
          for(j_cub[2]=j_cub_min[2];j_cub[2]<=j_cub_max[2];j_cub[2]++) {
            const cubicle_content_t* cub_j = &cubicles_.ref(j_cub);
            typename cubicle_content_t::const_iterator cub_ji;
            if (mi.i_sym == 0) {
              for(cub_ji=cub_j->begin();
                  cub_ji!=cub_j->end();
                  cub_ji++) {
                scitbx::vec3<FloatType> diff_vec =
                    mappings[cub_ji->i_seq][cub_ji->i_sym].mapped_site()
                  - mappings_i[0].mapped_site();
                if (diff_vec.length_sq() > distance_cutoff_sq_) continue;
                pair_asu_table_->add_pair(
                  mi.i_seq,
                  cub_ji->i_seq,
                  rt_mx_i_inverse.multiply(
                    asu_mappings_->get_rt_mx(cub_ji->i_seq, cub_ji->i_sym)));
              }
            }
            else {
              for(cub_ji=cub_j->begin();
                  cub_ji!=cub_j->end();
                  cub_ji++) {
                if (cub_ji->i_sym != 0) continue;
                scitbx::vec3<FloatType> diff_vec =
                    mappings_i[mi.i_sym].mapped_site()
                  - mappings[cub_ji->i_seq][0].mapped_site();
                if (diff_vec.length_sq() > distance_cutoff_sq_) continue;
                pair_asu_table_->add_pair(
                  cub_ji->i_seq,
                  mi.i_seq,
                  asu_mappings_->get_rt_mx(cub_ji->i_seq, 0)
                    .inverse().multiply(
                      asu_mappings_->get_rt_mx(mi.i_seq, mi.i_sym)));
              }
            }
          }}}
          std::size_t i1d_cub = cubicles_.ref.accessor()(i_cub);
          cubicles_.ref[i1d_cub].push_back(mi);
        }
      }

      void
      process_site_frac(
        fractional<FloatType> const& original_site)
      {
        process_site_frac(
          original_site,
          sgtbx::site_symmetry(
            asu_mappings_->asu().unit_cell(),
            asu_mappings_->space_group(),
            original_site,
            min_distance_sym_equiv,
            assert_min_distance_sym_equiv));
      }

      //! Not available in Python.
      void
      reserve_additional(std::size_t n)
      {
        n += asu_mappings_->mappings_const_ref().size();
        asu_mappings_->reserve(n);
        pair_asu_table_->reserve(n);
      }

      void
      process_sites_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        sgtbx::site_symmetry_table const& site_symmetry_table)
      {
        CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
                  == original_sites.size());
        reserve_additional(original_sites.size());
        for(std::size_t i=0;i<original_sites.size();i++) {
          process_site_frac(original_sites[i], site_symmetry_table.get(i));
        }
      }

      void
      process_sites_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites)
      {
        reserve_additional(original_sites.size());
        for(std::size_t i=0;i<original_sites.size();i++) {
          process_site_frac(original_sites[i]);
        }
      }

      void
      process_sites_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        sgtbx::site_symmetry_table const& site_symmetry_table)
      {
        CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
                  == original_sites.size());
        reserve_additional(original_sites.size());
        uctbx::unit_cell const& uc = asu_mappings_->asu().unit_cell();
        for(std::size_t i=0;i<original_sites.size();i++) {
          process_site_frac(
            uc.fractionalize(original_sites[i]), site_symmetry_table.get(i));
        }
      }

      void
      process_sites_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites)
      {
        reserve_additional(original_sites.size());
        uctbx::unit_cell const& uc = asu_mappings_->asu().unit_cell();
        for(std::size_t i=0;i<original_sites.size();i++) {
          process_site_frac(uc.fractionalize(original_sites[i]));
        }
      }

      std::map<long, long>
      cubicle_size_counts() const { return cubicles_.cubicle_size_counts(); }
  };

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_INCREMENTAL_PAIRS_H
