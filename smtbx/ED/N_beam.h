#pragma once
#include <cctbx/uctbx.h>
#include <cctbx/xray/thickness.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/miller/index_generator.h>

#include <smtbx/ED/ed_data.h>
#include <smtbx/ED/kinematic.h>
#include <smtbx/ED/dyn_calculator.h>
#include <smtbx/refinement/least_squares.h>

#define N_NEAM_SHARED_DATA_TYPEDES                                              \
  typedef f_calc_function_base<FloatType> f_calc_function_base_t;               \
  typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;\
  typedef builder_base<FloatType> data_t;                                       \
  typedef af::shared<const BeamInfo<FloatType>*> beam_at;                       \
  typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;      \
  typedef boost::shared_ptr<lookup_t> lookup_ptr_t;                             \
  typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;                \
  typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;              \
  typedef typename utils<FloatType>::ReflectionSlice ref_slice_t;

namespace smtbx { namespace ED
{
  using namespace cctbx;
  using namespace smtbx::refinement::least_squares;

  template <typename FloatType>
  struct N_beam_shared_data_base {
    ED_UTIL_TYPEDEFS;
    N_NEAM_SHARED_DATA_TYPEDES;

  protected:
    void init(af::shared<BeamGroup<FloatType> > beam_groups) {
      SMTBX_ASSERT(params.getBeamN() > 2);
      this->beam_groups = beam_groups;
      const typename utils<FloatType>::a_geometry& geom = *beam_groups[0].geometry;
      K = geom.Kl_as_K(Kl);
      af::shared<miller::index<> > groups_indices;
      // update min_d as the 'complete' list must have all of the measured reflections
      FloatType max_d_start_sq = 1./scitbx::fn::pow2(params.getTopUpD());
      for (size_t i = 0; i < beam_groups.size(); i++) {
        BeamGroup<FloatType>& beam_group = beam_groups[i];
        for (size_t hi = 0; hi < beam_group.strong_beams.size(); hi++) {
          const miller::index<>& h = beam_group.indices[beam_group.strong_beams[hi]];
          groups_indices.push_back(h);
          FloatType d_star_sq = unit_cell.d_star_sq(h);
          if (d_star_sq > max_d_start_sq) {
            max_d_start_sq = d_star_sq;
          }
        }
      }
      //// a single complete list of reflections up to d
      //{
      //  miller::index<> h;
      //  miller::index_generator h_generator(unit_cell, sgtbx::space_group_type("P1"),
      //    true,
      //    std::min(params.getTopUpD(), 1./std::sqrt(max_d_start_sq)) - 0.001,
      //    false);
      //  typedef typename utils<FloatType>::Reflection ref_t;
      //  while (!(h = h_generator.next()).is_zero()) {
      //    if (space_group.is_sys_absent(h)) {
      //      continue;
      //    }
      //    complete_sorted_set.push_back(
      //      new ref_t(h, std::abs(geom.get_diffraction_angle(h, K)))
      //    );
      //  }
      //  std::sort(
      //    complete_sorted_set.begin(), complete_sorted_set.end(),
      //    &ref_t::cmp_ptrs);
      //}
      // a single complete list of reflections up to d/2
      {
        miller::index<> h;
        miller::index_generator h_generator(unit_cell, sgtbx::space_group_type("P1"),
          true,
          (std::min(params.getTopUpD(), 1. / std::sqrt(max_d_start_sq)))/2 - 0.001,
          false);
        FloatType set_max_d_star_sq = std::max(1. / params.getTopUpD(), max_d_start_sq);
        typedef typename utils<FloatType>::Reflection ref_t;
        while (!(h = h_generator.next()).is_zero()) {
          indices.push_back(h);
          if (space_group.is_sys_absent(h)) {
            continue;
          }
          FloatType d_star_sq = unit_cell.d_star_sq(h);
          if (d_star_sq <= set_max_d_star_sq) {
            complete_sorted_set.push_back(
              new ref_t(h, std::abs(geom.get_diffraction_angle(h, K)))
            );
          }
        }
        std::sort(
          complete_sorted_set.begin(), complete_sorted_set.end(),
          &ref_t::cmp_ptrs);
      }
      sgtbx::space_group P1("P1");
      mi_lookup = lookup_t(
        indices.const_ref(),
        P1,
        true);
    }

    size_t find_sorted(const miller::index<>& h, FloatType da0,
      bool return_insertion_point) const
    {
      size_t from = 0, to = complete_sorted_set.size() - 1;
      if (complete_sorted_set[0]->da == da0) {
        return 0;
      }
      if (complete_sorted_set[to]->da == da0) {
        return to;
      }
      while ((to - from) != 1) {
        size_t index = from + (to - from) / 2;
        FloatType da = complete_sorted_set[index]->da;
        if (da < da0) {
          from = index;
        }
        else if (da > da0) {
          to = index;
        }
        else {
          return index;
        }
      }
      return return_insertion_point ? to : ~0;
    }
  public:
    N_beam_shared_data_base(
      f_calc_function_base_t& f_calc_function,
      sgtbx::space_group const& space_group,
      uctbx::unit_cell const& unit_cell,
      cctbx::xray::thickness<FloatType> const& thickness,
      const RefinementParams<FloatType>& params)
      : f_calc_function(f_calc_function),
      space_group(space_group),
      unit_cell(unit_cell),
      params(params),
      Kl(params.getKl()),
      Fc2Ug(params.getFc2Ug()),
      thickness(thickness),
      thread_n(params.getThreadN())
    {
      mtx = boost::shared_ptr<boost::mutex>(new boost::mutex());
    }

    ~N_beam_shared_data_base() {
      for (size_t i = 0; i < complete_sorted_set.size(); i++) {
        delete complete_sorted_set[i];
      }
    }

    void do_build_kin_mt() {
      if (thread_n < 0) {
        thread_n = builder_base<FloatType>::get_available_threads();
      }
      build_kin_mt(thread_n, Fc2Ug, f_calc_function, indices, Fcs_kin);
    }

    void build() {
      if (Fcs_kin.size() != indices.size()) {
        Fcs_kin.resize(indices.size());
        do_build_kin_mt();
        for (size_t i = 0; i < complete_sorted_set.size(); i++) {
          int j = mi_lookup.find_hkl(complete_sorted_set[i]->h);
          SMTBX_ASSERT(j >= 0);
          complete_sorted_set[i]->Ug = Fcs_kin[j];
        }
      }
    }

    af::shared<size_t> fetch_selection(const miller::index<>& h, const mat3_t& RM) const {
      af::shared<size_t> selection;
      FloatType da0 = std::abs(beam_groups[0].geometry->get_diffraction_angle(h, K));
      size_t sz = complete_sorted_set.size();
      size_t pos = find_sorted(h, da0, false);
      SMTBX_ASSERT(pos != ~0)(h.as_string());
      FloatType max_sg = params.getTopUpMaxSg(),
        ang_span = scitbx::deg_as_rad(params.getGroupWidth());
      selection.push_back(pos);

      /* alternative bracketing
      FloatType m_ang1 = scitbx::rad_as_deg(beam_groups[0].Sg_to_angle(-max_sg * 2, h, K));
      FloatType m_ang2 = scitbx::rad_as_deg(beam_groups[0].Sg_to_angle(max_sg * 2, h, K));
      */

      // go up
      size_t j = pos;
      //while (++j < sz) {
      while (++j < sz && (complete_sorted_set[j]->da - da0) <= ang_span) {
        FloatType sg = std::abs(utils<FloatType>::calc_Sg(RM, complete_sorted_set[j]->h, K));
        if (sg < max_sg) {
          selection.push_back(j);
        }
      }
      // go down
      j = pos;
      //while (--j != ~0) {
      while (--j != ~0 && (da0 - complete_sorted_set[j]->da) <= ang_span) {
        FloatType sg = std::abs(utils<FloatType>::calc_Sg(RM, complete_sorted_set[j]->h, K));
        if (sg < max_sg) {
          selection.push_back(j);
        }
      }
      return selection;
    }

    ref_slice_t get_selection(const miller::index<> &h, const mat3_t &RM, bool cache) const {
      if (cache) {
        typedef std::map<miller::index<>, ref_slice_t, miller::fast_less_than<> >
          ref_sel_cache_t;
        typename ref_sel_cache_t::const_iterator fi =
          selection_cache.find(h);
        if (fi != selection_cache.end()) {
          return fi->second;
        }
      }
      af::shared<size_t> selection = fetch_selection(h, RM);
      ref_slice_t rv(complete_sorted_set, selection);
      if (cache) {
        boost::mutex::scoped_lock lock(*mtx);
        selection_cache.insert(std::make_pair(h, rv));
      }
      return rv;
    }

    void cache_selection(const miller::index<>& h, const af::shared<mat3_t>& RMfs,
      bool reset) const
    {
      typedef std::map<miller::index<>, ref_slice_t, miller::fast_less_than<> >
        ref_sel_cache_t;
      typename ref_sel_cache_t::const_iterator fi =
        selection_cache.find(h);
      if (!reset && fi != selection_cache.end()) {
        return;
      }
      af::shared<size_t> all;
      for (size_t i = 0; i < RMfs.size(); i++) {
        af::shared<size_t> s = fetch_selection(h, RMfs[i]);
        all.extend(s.begin(), s.end());
      }
      std::sort(all.begin(), all.end());
      af::shared<size_t> selection;
      for (size_t i = 0; i < all.size(); i++) {
        selection.push_back(all[i]);
        size_t j = i + 1;
        while (j < all.size() && all[i] == all[j]) {
          j++;
        }
        i = j - 1;
      }
      boost::mutex::scoped_lock lock(*mtx);
      if (fi != selection_cache.end()) {
        selection_cache.erase(fi->first);
      }
      ref_slice_t rv(complete_sorted_set, selection);
      selection_cache.insert(std::make_pair(h, rv));
    }

    f_calc_function_base_t& f_calc_function;
    sgtbx::space_group const& space_group;
    uctbx::unit_cell const& unit_cell;
    af::shared<miller::index<> > indices;
    RefinementParams<FloatType> params;
    FloatType Kl, Fc2Ug;
    cart_t K;
    mutable std::map<miller::index<>, ref_slice_t, miller::fast_less_than<> >
      selection_cache;
    af::shared<BeamGroup<FloatType> > beam_groups;
    cctbx::xray::thickness<FloatType> const& thickness;

    af::shared<typename utils<FloatType>::Reflection*> complete_sorted_set;

    af::shared<complex_t> Fcs_kin;
    // 
    lookup_t mi_lookup;
    int thread_n;
    mutable boost::shared_ptr<boost::mutex> mtx;
  };

  template <typename FloatType>
  class dyn_calculator_n_beam {
  public:
    ED_UTIL_TYPEDEFS;
    typedef typename utils<FloatType>::ReflectionSlice ref_slice_t;
    typedef N_beam_shared_data_base<FloatType> data_t;

    dyn_calculator_n_beam(
      const data_t& data,
      const BeamGroup<FloatType>& beam_group,
      int mat_type)
      : data(data),
      dc_f(mat_type),
      beam_group(beam_group), beam_n(data.params.getBeamN()),
      K(data.K),
      thickness(data.thickness.value),
      wght(data.params.getNBeamWght()),
      useSg(data.params.useNBeamSg())
    {}

    /* builds the potential matrix in dc; init must be called before this
    function!! */
    void build(const mat3_t &R) {
      dc->reset(A, R);
    }

    complex_t calc_amp(const mat3_t& R, size_t idx=1) {
      return dc->reset(A, R)
        .calc_amps_1(idx);
    }

    //D_dyn has one row as output
    complex_t calc_amp_ext(const mat3_t& R,
      const af::shared<cmat_t>& Ds_kin,
      bool grad_thickness,
      mat_t& D_dyn)
    {
      return dc->reset(A, R)
        .calc_amps_ext_1(Ds_kin, grad_thickness, D_dyn, 1);
    }

    // recomputes the Eigen matrix
    dyn_calculator_n_beam& init(const miller::index<> &h, FloatType angle,
      const af::shared<complex_t> & Fcs_kin, const lookup_t& mi_lookup)
    {
      return init(h, beam_group.get_R(angle), Fcs_kin, mi_lookup);
    }

    // just caches the reflection selection
    dyn_calculator_n_beam& init(const miller::index<>& h,
      const af::shared<FloatType>& angles)
    {
      af::shared<mat3_t> RMfs(af::reserve(angles.size()));
      for (size_t i = 0; i < angles.size(); i++) {
        RMfs.push_back(beam_group.get_R(angles[i]));
      }
      return init(h, RMfs);
    }

    // recomputes the Eigen matrix
    dyn_calculator_n_beam& init(const miller::index<>& h,
      const mat3_t& RMf,
      const af::shared<complex_t>& Fcs_kin, const lookup_t& mi_lookup)
    {
      ref_slice_t reflections = data.get_selection(h, RMf, true);
      indices = utils<FloatType>::build_Ug_matrix_N(A, Fcs_kin, mi_lookup,
        reflections, K, h, RMf, beam_n, useSg, wght);
      dc = dc_f.make(indices, K, beam_group.get_N(), thickness);
      return *this;
    }
    
    // just initailises the reflection selection cache
    dyn_calculator_n_beam& init(const miller::index<>& h,
      const af::shared<mat3_t>& RMfs)
    {
      data.cache_selection(h, RMfs, false);
      return *this;
    }

    // indices selected for the Ug matrix - intialised by 'build'
    af::shared<miller::index<> > indices;
    const cmat_t& get_matrix() const {
      return A;
    }
    a_dyn_calculator<FloatType>&  get_dc() const {
      return *dc;
    }
  protected:
    const data_t& data;
    dyn_calculator_factory<FloatType> dc_f;
    const BeamGroup<FloatType>& beam_group;
    boost::shared_ptr<a_dyn_calculator<FloatType> > dc;
    size_t beam_n;
    cmat_t A;
    cart_t K;
    FloatType thickness, wght;
    bool useSg;
  };

}}
