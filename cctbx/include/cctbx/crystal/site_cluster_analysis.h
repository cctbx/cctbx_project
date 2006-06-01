#ifndef CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_H
#define CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_H

#include <cctbx/crystal/direct_space_asu.h>
#include <cctbx/crystal/cubicles.h>

namespace cctbx { namespace crystal {

  template <typename FloatType=double, typename IntShiftType=int>
  class site_cluster_analysis
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
      neighbors::cubicles<cubicle_content_t, FloatType> cubicles_;

    public:
      //! Default constructor. Some data members are not initialized!
      site_cluster_analysis() {}

      site_cluster_analysis(
        sgtbx::space_group const& space_group,
        direct_space_asu::float_asu<FloatType> const& asu,
        FloatType const& distance_cutoff,
        bool general_positions_only_=false,
        unsigned estimated_reduction_factor_=4,
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
        general_positions_only(general_positions_only_),
        estimated_reduction_factor(estimated_reduction_factor_),
        min_distance_sym_equiv(0.5),
        assert_min_distance_sym_equiv(true)
      {
        CCTBX_ASSERT(distance_cutoff_ > 0);
        CCTBX_ASSERT(asu_mappings_->buffer_thickness() >= distance_cutoff_);
      }

      boost::shared_ptr<
        direct_space_asu::asu_mappings<FloatType, IntShiftType> >
      asu_mappings() const { return asu_mappings_owner_; }

      bool general_positions_only;
      unsigned estimated_reduction_factor;
      FloatType min_distance_sym_equiv;
      bool assert_min_distance_sym_equiv;

    protected:
      void
      discard_new(std::vector<std::size_t> const& registry_new)
      {
        asu_mappings_->discard_last();
        for(std::size_t i=registry_new.size();i!=0;) {
          cubicles_.ref[registry_new[--i]].pop_back();
        }
      }

    public:
      bool
      process_site_frac(
        fractional<FloatType> const& original_site,
        sgtbx::site_symmetry_ops const& site_symmetry_ops)
      {
        if (general_positions_only
            && !site_symmetry_ops.is_point_group_1()) return false;
        direct_space_asu::asu_mapping_index mi;
        mi.i_seq = asu_mappings_->mappings().size();
        asu_mappings_->process(original_site, site_symmetry_ops);
        sgtbx::rt_mx
          rt_mx_i_inverse = asu_mappings_->get_rt_mx(mi.i_seq, 0).inverse();
        af::const_ref<
          typename asu_mappings_t::array_of_mappings_for_one_site> const&
            mappings = asu_mappings_->mappings_const_ref();
        typename asu_mappings_t::array_of_mappings_for_one_site const&
          mappings_i = mappings.back();
        scitbx::vec3<unsigned> n_cub = cubicles_.ref.accessor();
        std::vector<std::size_t> registry_new;
        registry_new.reserve(mappings_i.size());
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
                discard_new(registry_new);
                return false;
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
                discard_new(registry_new);
                return false;
              }
            }
          }}}
          std::size_t i1d_cub = cubicles_.ref.accessor()(i_cub);
          cubicles_.ref[i1d_cub].push_back(mi);
          registry_new.push_back(i1d_cub);
        }
        return true;
      }

      bool
      process_site_frac(
        fractional<FloatType> const& original_site)
      {
        return process_site_frac(
          original_site,
          sgtbx::site_symmetry(
            asu_mappings_->asu().unit_cell(),
            asu_mappings_->space_group(),
            original_site,
            min_distance_sym_equiv,
            assert_min_distance_sym_equiv));
      }

      //! Not available in Python.
      std::size_t
      reserve_additional(std::size_t n, bool apply_estimated_reduction_factor)
      {
        if (apply_estimated_reduction_factor
            && estimated_reduction_factor > 1) {
          n = (n+estimated_reduction_factor-1)/estimated_reduction_factor;
        }
        asu_mappings_->reserve(n + asu_mappings_->mappings_const_ref().size());
        return n;
      }

      af::shared<std::size_t>
      process_sites_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        std::size_t max_clusters=0)
      {
        CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
                  == original_sites.size());
#define CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_CONSTRUCT_SELECTION \
        af::shared<std::size_t> result; \
        if (max_clusters == 0) { \
          result.reserve(reserve_additional(original_sites.size(), true)); \
        } \
        else { \
          result.reserve(reserve_additional(max_clusters, false)); \
        }
        CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_CONSTRUCT_SELECTION
        for(std::size_t i=0;i<original_sites.size();i++) {
          if (process_site_frac(
                original_sites[i], site_symmetry_table.get(i))) {
            result.push_back(i);
            if (result.size() == max_clusters) break;
          }
        }
        return result;
      }

      af::shared<std::size_t>
      process_sites_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        std::size_t max_clusters=0)
      {
        CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_CONSTRUCT_SELECTION
        for(std::size_t i=0;i<original_sites.size();i++) {
          if (process_site_frac(original_sites[i])) {
            result.push_back(i);
            if (result.size() == max_clusters) break;
          }
        }
        return result;
      }

      af::shared<std::size_t>
      process_sites_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        std::size_t max_clusters=0)
      {
        CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
                  == original_sites.size());
        CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_CONSTRUCT_SELECTION
        scitbx::mat3<FloatType>
          frac = asu_mappings_->asu().unit_cell().fractionalization_matrix();
        for(std::size_t i=0;i<original_sites.size();i++) {
          if (process_site_frac(
                frac * original_sites[i], site_symmetry_table.get(i))) {
            result.push_back(i);
            if (result.size() == max_clusters) break;
          }
        }
        return result;
      }

      af::shared<std::size_t>
      process_sites_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        std::size_t max_clusters=0)
      {
        CCTBX_CRYSTAL_SITE_CLUSTER_ANALYSIS_CONSTRUCT_SELECTION
        scitbx::mat3<FloatType>
          frac = asu_mappings_->asu().unit_cell().fractionalization_matrix();
        for(std::size_t i=0;i<original_sites.size();i++) {
          if (process_site_frac(frac * original_sites[i])) {
            result.push_back(i);
            if (result.size() == max_clusters) break;
          }
        }
        return result;
      }
  };

}} // namespace cctbx::crystal

#endif // CCTBX_SITE_CLUSTER_ANALYSIS_H
