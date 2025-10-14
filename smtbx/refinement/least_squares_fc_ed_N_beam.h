#pragma once
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <scitbx/lstbx/normal_equations.h>
#include <smtbx/ED/n_beam.h>
#include <smtbx/refinement/weighting_schemes.h>

namespace smtbx {  namespace refinement  { namespace least_squares
{
  using namespace smtbx::ED;

  template <typename FloatType>
  struct beam_width_cache {
    typedef std::map<miller::index<>, FloatType> cache_t;
    cache_t cache;

    FloatType find_width(const miller::index<> &h) const {
      typename std::map<miller::index<>, FloatType>::const_iterator
        width_itr = cache.find(h);
      return width_itr != cache.end() ? width_itr->second : 0;
    }
  };

  template <typename FloatType>
  struct N_beam_shared_data : public N_beam_shared_data_base<FloatType> {
    ED_UTIL_TYPEDEFS;
    N_NEAM_SHARED_DATA_TYPEDES;
    typedef N_beam_shared_data_base<FloatType> base_t;

    N_beam_shared_data(const scitbx::sparse::matrix<FloatType>&
        Jt_matching_grad_fc,
      f_calc_function_base_t& f_calc_function,
      sgtbx::space_group const& space_group,
      uctbx::unit_cell const &unit_cell,
      af::shared<BeamGroup<FloatType> > beam_groups,
      cctbx::xray::thickness<FloatType> const& thickness,
      const RefinementParams<FloatType>& params,
      bool compute_grad,
      bool do_build = true)
      : base_t(f_calc_function, space_group, unit_cell, thickness, params),
      Jt_matching_grad_fc(Jt_matching_grad_fc),
      compute_grad(compute_grad)
    {
      base_t::init(beam_groups);
      for (size_t i = 0; i < beam_groups.size(); i++) {
        BeamGroup<FloatType>& beam_group = beam_groups[i];
        beam_groups_map.insert(std::make_pair(beam_group.id, &beam_group));
      }
      BeamGroup<FloatType>::link_groups(beam_groups);
      if (do_build) {
        build();
      }
    }

    ~N_beam_shared_data() {
    }

    void do_build_kin_mt() {
      if (this->thread_n < 0) {
        this->thread_n = builder_base<FloatType>::get_available_threads();
      }
      build_kin_mt_ext(this->thread_n, Jt_matching_grad_fc, this->Fc2Ug,
        this->f_calc_function,
        this->indices, this->Fcs_kin, design_matrix_kin, compute_grad);
    }

    void build() {
      if (this->Fcs_kin.size() != this->indices.size()) {
        if (compute_grad) {
          size_t cn = Jt_matching_grad_fc.n_rows() - (this->thickness.grad ? 1 : 0);
          if (cn > 0) {
            design_matrix_kin.resize(
              af::mat_grid(this->indices.size(), cn));
          }
        }
        this->Fcs_kin.resize(this->indices.size());
        do_build_kin_mt();
        for (size_t i = 0; i < this->complete_sorted_set.size(); i++) {
          int j = this->mi_lookup.find_hkl(this->complete_sorted_set[i]->h);
          SMTBX_ASSERT(j >= 0);
          this->complete_sorted_set[i]->Ug = this->Fcs_kin[j];
        }
      }
    }

    void build_width_cache(bool rebuild = false);

    /* Is is from compute_dynI - raw intensities + scales */
    FloatType compute_OSF(const af::shared<FloatType>& Is,
      const IWeightingScheme<FloatType>& ws,
      FloatType OSF, bool use_scales) const;
    
    /* Is is from compute_dynI - raw intensities + scales */
    FloatType compute_group_wR2(const af::shared<FloatType>& Is,
      const IWeightingScheme<FloatType>& ws,
      FloatType OSF,
      size_t group_idx) const;

    /* returns raw intensities and scales that would be applied for the L.S.
    */
    af::shared<FloatType> compute_dynI(
      const af::shared<miller::index<> > &indices) const;

    af::shared<FloatType> compute_linear_scales(
      const af::shared<miller::index<> >& indices) const;

    af::shared<FloatType> optimise_linear_scales(
      const IWeightingScheme<FloatType> &ws, FloatType OSF, bool compute_wR2);

    af::shared<FloatType> optimise_flat_scales(
      const IWeightingScheme<FloatType>& ws, FloatType OSF, bool compute_wR2);

    void set_width_cache(const beam_width_cache<FloatType>& width_cache) {
      this->width_cache = width_cache;
    }

    scitbx::sparse::matrix<FloatType> Jt_matching_grad_fc;

    typename std::map<int, BeamGroup<FloatType>*> beam_groups_map;
    bool compute_grad;
    // 
    cmat_t design_matrix_kin;
    beam_width_cache<FloatType> width_cache;
  };

  template <typename FloatType>
  class f_calc_function_ed_N_beam : public f_calc_function_base<FloatType> {
  public:
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef N_beam_shared_data<FloatType> data_t;
    typedef af::versa_plain<FloatType> one_dim_type;
    typedef typename one_dim_type::accessor_type one_dim_accessor_type;

    f_calc_function_ed_N_beam(data_t const& data)
      : data(data),
      index(-1),
      observable_updated(false),
      computed(false),
      raw_I(false)
    {
      beam_group = 0;
    }

    f_calc_function_ed_N_beam(f_calc_function_ed_N_beam const& other)
      : data(other.data),
      observable_updated(false),
      computed(false), raw_I(other.raw_I)
    {}

    virtual void compute(
      miller::index<> const& h,
      boost::optional<complex_t> const& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      SMTBX_ASSERT(fraction != 0);
      index = data.mi_lookup.find_hkl(h);
      if (index == -1) {
        if (!data.space_group.is_sys_absent(h)) {
          SMTBX_ASSERT(index >= 0)(h.as_string());
        }
        Fc = 0;
        Fsq = 0;
        observable_updated = true;
      }
      else {
        observable_updated = false;
        Fc = data.Fcs_kin[index];
        Fsq = std::norm(Fc);
      }
      typename std::map<int, BeamGroup<FloatType>*>::const_iterator fi =
        data.beam_groups_map.find(fraction->tag);
      SMTBX_ASSERT(fi != data.beam_groups_map.end());
      beam_group = fi->second;
      this->h = h;
      this->compute_grad = compute_grad;
      computed = true;
    }

    void setup_compute(
      const miller::index<> & h,
      const BeamGroup<FloatType>* beam_group)
    {
      SMTBX_ASSERT(beam_group != 0);
      this->beam_group = beam_group;
      this->h = h;
      this->compute_grad = false;
      computed = true;
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed_N_beam(*this));
    }

    f_calc_function_ed_N_beam<FloatType> * raw_fork() const {
      return new f_calc_function_ed_N_beam(*this);
    }

    FloatType get_observable_N() const {
      FloatType da = beam_group->geometry->get_diffraction_angle(h, data.K);
      af::shared<FloatType> angles;

      FloatType width_Sg = 0;
      typename std::map<miller::index<>, FloatType>::const_iterator
        width_itr = data.width_cache.cache.find(h);
      if (width_itr != data.width_cache.cache.end()) {
        width_Sg = width_itr->second;
      }
      //
      if (width_Sg > 0) {
        //angles = beam_group->get_angles_Sg(h, data.K,
        //  width_Sg,
        //  width_Sg * data.params.getIntStep() / data.params.getIntSpan());
        angles = beam_group->get_angles_Sg_N(h, data.K,
          width_Sg,
          data.params.getIntSpan() / data.params.getIntStep());
      }
      else {
        if (data.params.isAngleInt()) {
          angles = beam_group->get_angles(da,
            data.params.getIntSpan(),
            data.params.getIntStep());
        }
        else {
          angles = beam_group->get_angles_Sg(h, data.K,
            data.params.getIntSpan(),
            data.params.getIntStep());
        }
      }
      mat3_t da_rm = beam_group->geometry->get_RM(da);
      dyn_calculator_n_beam<FloatType> n_beam_dc(
        data, *beam_group, data.params.getMatrixType());

      if (!data.params.isNBeamFloating()) {
        n_beam_dc.init(h, beam_group->geometry->get_RMf(da_rm), data.Fcs_kin,
          data.mi_lookup);
      }
      else {
        n_beam_dc.init(h, angles);
      }

      mat_t D_dyn;
      af::shared<cmat_t> Ds_kin;
      if (compute_grad && !data.params.isNBeamFloating()) {
        utils<FloatType>::build_D_matrices(data.mi_lookup, n_beam_dc.indices,
          data.design_matrix_kin, Ds_kin);
      }
      size_t n_param = data.Jt_matching_grad_fc.n_rows();
      size_t coln = data.design_matrix_kin.accessor().n_columns() +
        (data.thickness.grad ? 1 : 0);
      af::shared<FloatType> grads_sum(coln),
        grads1;

      FloatType I1 = -1, g1 = -1, I_sum = 0, step_sum = 0,
        I = 0;
      for (size_t ai = 0; ai < angles.size(); ai++) {
        mat3_t R = beam_group->get_R(angles[ai]);
        cart_t g = R * cart_t(h[0], h[1], h[2]);
        cart_t K_g = g + data.K;
        FloatType K_g_l = K_g.length();
        if (compute_grad) {
          if (data.params.isNBeamFloating()) {
            n_beam_dc.init(h, R, data.Fcs_kin, data.mi_lookup);
            utils<FloatType>::build_D_matrices(data.mi_lookup, n_beam_dc.indices,
              data.design_matrix_kin, Ds_kin);
          }
          I = std::norm(
            n_beam_dc.calc_amp_ext(R, Ds_kin, data.thickness.grad, D_dyn)
          );
        }
        else {
          if (data.params.isNBeamFloating()) {
            n_beam_dc.init(h, R, data.Fcs_kin, data.mi_lookup);
          }
          I = std::norm(n_beam_dc.calc_amp(R));
        }
        if (g1 >= 0) {
          FloatType d = std::abs(K_g_l - g1) / 2;
          I_sum += (I + I1) * d;
          step_sum += d;
          if (compute_grad) {
            for (size_t i = 0; i < coln; i++) {
              grads_sum[i] += (grads1[i] + D_dyn(0, i)) * d;
            }
          }
        }
        I1 = I;
        g1 = K_g_l;
        if (compute_grad) {
          grads1 = af::shared<FloatType>(&D_dyn(0, 0), &D_dyn(0, n_param));
        }
      }
      if (compute_grad) {
        grads = grads_sum;
      }
      scale = compute_scale();
      if (!data.params.useFlatScales() && scale != 1 && !raw_I) {
        I_sum *= scale;
        if (compute_grad) {
          for (size_t i = 0; i < coln; i++) {
            grads_sum[i] *= scale;
          }
        }
      }
      return I_sum;
    }

    FloatType compute_scale() const {
      SMTBX_ASSERT(beam_group != 0);
      if (data.params.useFlatScales()) {
        return 1;
      }
      FloatType sc = 1;
      FloatType da = beam_group->geometry->get_diffraction_angle(h, data.K);
      if (beam_group->next != 0) {
        if (da > beam_group->angle && da < beam_group->next->angle) {
          FloatType Sg1 = std::abs(beam_group->calc_Sg(h, data.K));
          FloatType Sg2 = std::abs(beam_group->next->calc_Sg(h, data.K));
          FloatType Sgr = Sg1 + Sg2;
          sc = (Sg2 * beam_group->scale
            + Sg1 * beam_group->next->scale);
          sc /= Sgr * beam_group->scale;
        }
        return sc;
      }
      if (beam_group->prev != 0) {
        if (da < beam_group->angle && da > beam_group->prev->angle) {
          FloatType Sg1 = std::abs(beam_group->calc_Sg(h, data.K));
          FloatType Sg2 = std::abs(beam_group->prev->calc_Sg(h, data.K));
          FloatType Sgr = Sg1 + Sg2;
          sc = (Sg2 * beam_group->scale
            + Sg1 * beam_group->prev->scale);
          sc /= Sgr * beam_group->scale;
        }
      }
      return sc;
    }

    virtual FloatType get_observable() const {
      if (observable_updated || !computed) {
        return Fsq;
      }
      SMTBX_ASSERT(beam_group != 0);
      Fsq = get_observable_N();
      observable_updated = true;
      return Fsq;
    }

    virtual complex_t get_f_calc() const {
      if (!observable_updated) {
        get_observable();
      }
      return Fc;
    }

    virtual af::const_ref<complex_t> get_grad_f_calc() const {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
    }

    virtual af::const_ref<FloatType> get_grad_observable() const {
      if (!observable_updated) {
        get_observable();
      }
      return grads.const_ref();
    }

    virtual bool raw_gradients() const { return false; }

    FloatType compute_width(const miller::index<>& index,
      const BeamGroup<FloatType>& group,
      FloatType test_span, size_t test_points) const
    {
      FloatType da = group.geometry->get_diffraction_angle(index, data.K);
      af::shared<FloatType> angles = group.get_angles_Sg(index, data.K,
        test_span,
        test_span / test_points);

      mat3_t da_rm = group.geometry->get_RM(da);
      dyn_calculator_n_beam<FloatType> n_beam_dc(
        data, group, data.params.getMatrixType());

      if (!data.params.isNBeamFloating()) {
        n_beam_dc.init(index, group.geometry->get_RMf(da_rm), data.Fcs_kin,
          data.mi_lookup);
      }
      else {
        n_beam_dc.init(index, angles);
      }
      size_t start_idx = ~0, end_idx = ~0;
      af::shared<FloatType> amps(angles.size());
      FloatType maxI = 0, minI = 1e6;
      for (size_t ai = 0; ai < angles.size(); ai++) {
        mat3_t R = group.get_R(angles[ai]);
        if (data.params.isNBeamFloating()) {
          n_beam_dc.init(index, R, data.Fcs_kin, data.mi_lookup);
        }
        FloatType I = std::norm(n_beam_dc.calc_amp(R));
        if (I > maxI) {
          maxI = I;
        }
        if (I < minI) {
          minI = I;
        }
        amps[ai] = I;
      }
      FloatType I_th = maxI * data.params.getIntProfileStartTh();
      for (size_t i = 0; i < amps.size(); i++) {
        if ((amps[i] - minI) >= I_th) {
          if (start_idx == ~0) {
            start_idx = (i == 0 ? i : i - 1);
            if (end_idx != ~0) {
              break;
            }
          }
        }
        if ((amps[amps.size() - i - 1] - minI) >= I_th) {
          if (end_idx == ~0) {
            end_idx = (i > 0 ? amps.size() - i : amps.size() - 1);
            if (start_idx != ~0) {
              break;
            }
          }
        }
      }
      if (start_idx == ~0 || end_idx == ~0 || start_idx == end_idx) {
        return 0;
      }
      else {
        FloatType eps = maxI/1e6;
        // extend left & right
        while (start_idx != 0) {
          FloatType diff = amps[start_idx] - amps[start_idx-1];
          if (std::abs(diff) > eps && diff < 0) {
            break;
          }
          start_idx--;
        }
        while (end_idx + 1 < amps.size()) {
          FloatType diff = amps[end_idx] - amps[end_idx+1];
          if (std::abs(diff) > eps && diff < 0) {
            break;
          }
          end_idx++;
        }
        std::pair<FloatType, FloatType> k = group.Sg_to_angle_k(index, data.K);
        FloatType s = k.second * (angles[start_idx] - k.first);
        FloatType e = k.second * (angles[end_idx] - k.first);
        return std::max(std::abs(s), std::abs(e));
      }
    }
    const data_t& get_data() const {
      return data;
    }

    FloatType get_scale() const {
      return scale;
    }
  private:
    const data_t& data;
    long index;
    const BeamGroup<FloatType>* beam_group;
    bool compute_grad;
    mutable bool observable_updated, computed;
    mutable complex_t Fc;
    mutable FloatType Fsq, scale;
    mutable af::shared<FloatType> grads;
    miller::index<> h;
  public:
    bool raw_I;
  };


  template <typename FloatType>
  struct beam_width_thread {
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_ed_N_beam<FloatType> f_calc_function_t;
    typedef typename boost::shared_ptr<f_calc_function_t> f_calc_function_ptr_t;
    typedef std::map<miller::index<>, FloatType> cache_t;

    beam_width_thread(const BeamGroup<FloatType>& beam_group,
      const f_calc_function_t& f_calc_function)
    : beam_group(beam_group),
      f_calc_function(f_calc_function.raw_fork())
    {}
    
    void operator ()() {
      try {
        FloatType Sg_span = f_calc_function->get_data().params.getIntProfileSpan_Sg();
        size_t pts_N = f_calc_function->get_data().params.getIntProfilePoints();
        for (size_t bi = 0; bi < beam_group.beams.size(); bi++) {
          FloatType width = f_calc_function->compute_width(
            beam_group.beams[bi].h, beam_group, Sg_span, pts_N);
          cache.insert(std::make_pair(beam_group.beams[bi].h, width));
        }
      }
      catch (smtbx::error const& e) {
        exception_.reset(new smtbx::error(e));
      }
      catch (std::exception const& e) {
        exception_.reset(new smtbx::error(e.what()));
      }
    }

    boost::scoped_ptr<smtbx::error> exception_;
    const BeamGroup<FloatType>& beam_group;
    f_calc_function_ptr_t f_calc_function;
    cache_t cache;
  };

  template <typename FloatType>
  void N_beam_shared_data<FloatType>::build_width_cache(bool rebuild) {
    typedef beam_width_thread<FloatType> builder_t;
    typedef typename boost::shared_ptr<builder_t> build_bw_t;
    if (!rebuild && width_cache.cache.size() > 0) {
      return;
    }
    if (this->thread_n < 0) {
      this->thread_n = builder_base<FloatType>::get_available_threads();
    }
    boost::thread_group pool;
    f_calc_function_ed_N_beam<FloatType> calc_func(*this);
    size_t start = 0;
    while (start < this->beam_groups.size()) {
      std::vector<build_bw_t> accumulators;
      size_t end = std::min(this->beam_groups.size() - start, (size_t)this->thread_n);
      for (size_t th = 0; th < end; th++) {
        build_bw_t pf(
          new builder_t(this->beam_groups[start + th], calc_func)
        );
        accumulators.push_back(pf);
        pool.create_thread(boost::ref(*pf));
      }
      pool.join_all();
      for (size_t th = 0; th < end; th++) {
        if (accumulators[th]->exception_) {
          throw* accumulators[th]->exception_.get();
        }
        width_cache.cache.insert(
          accumulators[th]->cache.begin(),
          accumulators[th]->cache.end());
      }
      start += this->thread_n;
    }
  }

  template <typename FloatType>
  struct dynI_thread {
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_ed_N_beam<FloatType> f_calc_function_t;
    typedef typename boost::shared_ptr<f_calc_function_t> f_calc_function_ptr_t;
    typedef std::map<miller::index<>, FloatType> cache_t;

    dynI_thread(const af::shared<miller::index<> >& indices,
      size_t start, size_t end,
      const std::map<miller::index<>, const BeamGroup<FloatType>*>& lookup,
      const f_calc_function_t& f_calc_function)
      : indices(indices),
      start(start), end(end),
      lookup(lookup),
      f_calc_function(f_calc_function.raw_fork())
    {
      Is.reserve(end - start);
      scales.reserve(end - start);
    }

    void operator ()() {
      try {
        f_calc_function_t& f_calc = *f_calc_function;
        for (size_t i = start; i < end; i++) {
          const miller::index<>& h = indices[i];
          typename std::map<miller::index<>, const BeamGroup<FloatType>*>::const_iterator fi =
            lookup.find(h);
          SMTBX_ASSERT(fi != lookup.end());
          f_calc.setup_compute(h, fi->second);
          FloatType I = f_calc.get_observable_N();
          Is.push_back(I);
          scales.push_back(f_calc.get_scale());
        }
      }
      catch (smtbx::error const& e) {
        exception_.reset(new smtbx::error(e));
      }
      catch (std::exception const& e) {
        exception_.reset(new smtbx::error(e.what()));
      }
    }
    boost::scoped_ptr<smtbx::error> exception_;
    const af::shared<miller::index<> >& indices;
    size_t start, end;
    const std::map<miller::index<>, const BeamGroup<FloatType>*>& lookup;
    f_calc_function_ptr_t f_calc_function;
    std::vector<FloatType> Is, scales;
  };

  template <typename FloatType>
  af::shared<FloatType> N_beam_shared_data<FloatType>::compute_dynI(
    const af::shared<miller::index<> >& indices) const
  {
    typedef dynI_thread<FloatType> thread_t;
    typedef typename boost::shared_ptr<thread_t> build_bw_t;
    int thN = this->thread_n < 0 ? builder_base<FloatType>::get_available_threads()
      : this->thread_n;
    
    std::map<miller::index<>, const BeamGroup<FloatType>*> lookup;
    for (size_t i = 0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      for (size_t bi = 0; bi < group.beams.size(); bi++) {
        lookup.insert(std::make_pair(group.beams[bi].h, &group));
      }
    }
    boost::thread_group pool;
    f_calc_function_ed_N_beam<FloatType> calc_func(*this);
    calc_func.raw_I = true;
    size_t t_cnt = (indices.size() + thN - 1) / thN;
    std::vector<build_bw_t> accumulators;
    for (size_t i = 0; i < indices.size(); i += t_cnt) {
      build_bw_t pf(
        new thread_t(indices, i, std::min(i + t_cnt, indices.size()), lookup, calc_func)
      );
      accumulators.push_back(pf);
      pool.create_thread(boost::ref(*pf));
    }
    pool.join_all();
    af::shared<FloatType> rv(indices.size() * 2);
    size_t st = 0;
    for (size_t i = 0; i < accumulators.size(); i++) {
      if (accumulators[i]->exception_) {
        throw* accumulators[i]->exception_.get();
      }
      memcpy(&rv[st], &accumulators[i]->Is[0],
        accumulators[i]->Is.size() * sizeof(FloatType));
      st += accumulators[i]->Is.size();
    }
    for (size_t i = 0; i < accumulators.size(); i++) {
      memcpy(&rv[st], &accumulators[i]->scales[0],
        accumulators[i]->scales.size() * sizeof(FloatType));
      st += accumulators[i]->scales.size();
    }
    return rv;
  }

  template <typename FloatType>
  af::shared<FloatType> N_beam_shared_data<FloatType>::
    compute_linear_scales(
      const af::shared<miller::index<> >& indices) const
  {
    af::shared<FloatType> rv(indices.size());
    f_calc_function_ed_N_beam<FloatType> calc_func(*this);
    for (size_t i = 0, ii = 0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      for (size_t bi = 0; bi < group.beams.size(); bi++) {
        calc_func.setup_compute(group.beams[bi].h, &group);
        rv[ii] = calc_func.compute_scale();
      }
    }
    return rv;
  }

  template <typename FloatType>
  FloatType N_beam_shared_data<FloatType>::
    compute_OSF(const af::shared<FloatType> &Is,
      const IWeightingScheme<FloatType>& ws,
      FloatType OSF, bool use_scales) const
  {
    FloatType yo_dot_yc = 0, yc_sq = 0;
    size_t scales_off = Is.size() / 2;
    for (size_t i = 0, ii = 0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      for (size_t bi = 0; bi < group.beams.size(); bi++, ii++) {
        //SMTBX_ASSERT(ii < Is.size());
        FloatType Ic = Is[ii], Io = group.beams[bi].I;
        if (use_scales) {
          Ic *= Is[scales_off + ii];
        }
        FloatType w = ws.compute(Io, group.beams[bi].sig, Ic, group.beams[bi].h, OSF);
        yo_dot_yc += w * Io * Ic;
        yc_sq += w * Ic * Ic;
      }
    }
    return yo_dot_yc / yc_sq;
  }

  template <typename FloatType>
  FloatType N_beam_shared_data<FloatType>::
    compute_group_wR2(const af::shared<FloatType>& Is,
      const IWeightingScheme<FloatType>& ws,
      FloatType OSF,
      size_t group_idx) const
  {
    SMTBX_ASSERT(group_idx < this->beam_groups.size());
    size_t off = 0;
    for (size_t i = 0; i < this->beam_groups.size(); i++) {
      if (i == group_idx) {
        break;
      }
      off += this->beam_groups[i].beams.size();
    }
    size_t sz = Is.size() / 2;
    const BeamGroup<FloatType>& group = this->beam_groups[group_idx];
    FloatType wR2_up = 0, wR2_dn = 0;
    for (size_t bi = 0; bi < group.beams.size(); bi++) {
      FloatType Ic = Is[off+bi], Io = group.beams[bi].I;
      FloatType w = ws.compute(Io, group.beams[bi].sig, Ic, group.beams[bi].h, OSF);
      FloatType x = Ic * OSF * Is[sz + off + bi] * group.scale;
      wR2_up += w * scitbx::fn::pow2(Io - x);
      wR2_dn += w * x * x;
    }
    return std::sqrt(wR2_up / wR2_dn);
  }

  template <typename FloatType>
  af::shared<FloatType> N_beam_shared_data<FloatType>::
    optimise_linear_scales(const IWeightingScheme<FloatType>& ws, FloatType OSF,
      bool compute_wR2)
  {
    af::shared<miller::index<> > idx_obs;
    for (size_t i = 0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      for (size_t bi = 0; bi < group.beams.size(); bi++) {
        idx_obs.push_back(group.beams[bi].h);
      }
    }
    af::shared<FloatType> Is = this->compute_dynI(idx_obs);
    size_t sz = idx_obs.size();
    FloatType bare_OSF = compute_OSF(Is, ws, OSF, false);
    scitbx::lstbx::normal_equations::linear_ls<FloatType> ls(this->beam_groups.size());
    for (size_t i = 0, ii=0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      size_t beam_sz = group.beams.size();
      for (size_t bi = 0; bi < beam_sz; bi++, ii++) {
        const miller::index<>& h = idx_obs[ii];
        FloatType da = group.geometry->get_diffraction_angle(h, this->K);
        FloatType Sg1 = std::abs(group.calc_Sg(h, this->K)),
          I = bare_OSF * Is[ii];

        af::shared<FloatType> row(this->beam_groups.size());
        if (group.prev != 0 && da < group.angle && da > group.prev->angle) {
          FloatType Sg2 = std::abs(group.prev->calc_Sg(h, this->K));
          FloatType Sgr = Sg1 + Sg2;
          row[group.id-1] = I * Sg2 / Sgr;
          row[group.prev->id-1] = I * Sg1 / Sgr;
          I = row[group.id - 1] * group.scale + row[group.prev->id - 1] * group.prev->scale;
        }
        else if (group.next != 0 && da > group.angle && da < group.next->angle) {
          FloatType Sg2 = std::abs(group.next->calc_Sg(h, this->K));
          FloatType Sgr = Sg1 + Sg2;
          row[group.id-1] = I * Sg2 / Sgr;
          row[group.next->id-1] = I * Sg1 / Sgr;
          I = row[group.id - 1] * group.scale + row[group.next->id - 1] * group.next->scale;
        }
        else {
          row[group.id - 1] = I;
        }
        ls.add_equation(group.beams[bi].I, row.const_ref(),
          ws(group.beams[bi].I, group.beams[bi].sig, I / bare_OSF, group.beams[bi].h, bare_OSF));
      }
    }
    ls.solve();
    af::shared<FloatType> rv = ls.solution();
    for (size_t i = 0; i < this->beam_groups.size(); i++) {
      this->beam_groups[i].scale = rv[i];
    }
    FloatType new_OSF = compute_OSF(Is, ws, OSF, true);
    if (compute_wR2) {
      af::shared<FloatType> l_scales = this->compute_linear_scales(idx_obs);
      for (size_t i = 0; i < sz; i++) {
        Is[sz + i] = l_scales[i];
      }
      for (size_t i = 0; i < this->beam_groups.size(); i++) {
        rv.push_back(
          compute_group_wR2(Is, ws, new_OSF, i));
      }
    }
    rv.push_back(new_OSF);
    return rv;
  }

  template <typename FloatType>
  af::shared<FloatType> N_beam_shared_data<FloatType>::
    optimise_flat_scales(const IWeightingScheme<FloatType>& ws, FloatType OSF,
      bool compute_wR2)
  {
    af::shared<miller::index<> > idx_obs;
    for (size_t i = 0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      for (size_t bi = 0; bi < group.beams.size(); bi++) {
        idx_obs.push_back(group.beams[bi].h);
      }
    }
    af::shared<FloatType> Is = this->compute_dynI(idx_obs);
    size_t sz = idx_obs.size();
    FloatType bare_OSF = compute_OSF(Is, ws, OSF, false);
    af::shared<FloatType> rv(af::reserve(this->beam_groups.size()));
    for (size_t i = 0, ii = 0; i < this->beam_groups.size(); i++) {
      const BeamGroup<FloatType>& group = this->beam_groups[i];
      size_t beam_sz = group.beams.size();
      FloatType left = 0, right = 0;
      for (size_t bi = 0; bi < beam_sz; bi++, ii++) {
        const miller::index<>& h = idx_obs[ii];
        FloatType I = bare_OSF * Is[ii];
        FloatType w = ws.compute(group.beams[bi].I, group.beams[bi].sig, Is[ii], h, bare_OSF);
        right += w * I * group.beams[bi].I;
        left += w * I * I;
      }
      rv.push_back(right / left);
    }
    for (size_t i = 0; i < this->beam_groups.size(); i++) {
      this->beam_groups[i].scale = rv[i];
    }
    FloatType new_OSF = bare_OSF;// compute_OSF(Is, ws, OSF, true);
    if (compute_wR2) {
      for (size_t i = 0; i < this->beam_groups.size(); i++) {
        rv.push_back(
          compute_group_wR2(Is, ws, new_OSF, i));
      }
    }
    rv.push_back(new_OSF);
    return rv;
  }
}}}
