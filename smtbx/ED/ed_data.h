#pragma once
#include <smtbx/ED/utils.h>
#include <scitbx/constants.h>

namespace smtbx { namespace ED
{
  using namespace cctbx;

  template <typename FloatType> struct BeamInfo;

  template <typename FloatType>
  class FrameInfo {
  public:
    ED_UTIL_TYPEDEFS;

    FrameInfo() {}

    FrameInfo(int id, const cart_t &c_normal,
      FloatType alpha, FloatType beta, FloatType omega,
      FloatType angle, FloatType scale, mat3_t const& UB)
      : id(id), tag(-1), original_normal(c_normal),
      UB(UB),
      alpha(alpha), beta(beta), omega(omega),
      angle(angle), scale(scale),
      offset(~0)
    {
      update_alpha(alpha);
    }

    std::pair<mat3_t, cart_t> compute_RMf_N(FloatType alpha_) const {
      FloatType ca = std::cos(alpha_), sa = std::sin(alpha_),
        cb = std::cos(beta), sb = std::sin(beta),
        co = std::cos(omega), so = std::sin(omega);
      mat3_t rxa(1, 0, 0, 0, ca, -sa, 0, sa, ca),
        ryb(cb, 0, sb, 0, 1, 0, -sb, 0, cb),
        rzo(co, -so, 0, so, co, 0, 0, 0, 1);
      mat3_t rm = rzo * rxa * ryb;
      mat3_t rmf = rm * UB;
      cart_t N = rm * original_normal;
      return std::make_pair(rmf, N / N.length());
    }

    void update_alpha(FloatType alpha_) {
      alpha = alpha_;
      std::pair<mat3_t, cart_t> r = compute_RMf_N(alpha_);
      RMf = r.first;
      normal = r.second;
    }

    void update_angles(FloatType alpha_, FloatType beta_, FloatType omega_) {
      alpha = alpha_;
      beta = beta_;
      omega = omega_;
      std::pair<mat3_t, cart_t> r = compute_RMf_N(alpha_);
      RMf = r.first;
      normal = r.second;
    }

    // returns angle in rads at which the excitation error is 0
    FloatType get_diffraction_angle(const miller::index<>& h,
      const cart_t& K, FloatType sweep_angle=3) const
    {
      return Sg_to_angle(0, h, K, sweep_angle);
    }

    /* returns angle in rads at which the excitation error is Sg as
    angle = alpha + (Sg-rv.first)/rv.second
    */
    std::pair<FloatType, FloatType> Sg_to_angle_k(
      const miller::index<>& h,
      const cart_t& K, FloatType sweep_angle = 3) const
    {
#ifdef _DEBUG
      SMTBX_ASSERT(std::abs(K.length() - std::abs(K[2])) < 1e-6);
#endif
      FloatType Sg1 = utils<FloatType>::calc_Sg(RMf * h, K);
      FloatType ang_diff = scitbx::deg_as_rad(sweep_angle);
      std::pair<mat3_t, cart_t> r = compute_RMf_N(alpha + ang_diff);
      FloatType Sg2 = utils<FloatType>::calc_Sg(r.first * h, K);
      FloatType k = (Sg2 - Sg1) / ang_diff;
      //FloatType a = Sg1 - k*alpha;
      //return (Sg - a) / k = (Sg + k*alpha - sg1)/k = alpha + (Sg - Sg1)/k;
      return std::make_pair(Sg1, k);
    }

    /* returns angle in rads at which the excitation error is Sg */
    FloatType Sg_to_angle(FloatType Sg, const miller::index<>& h,
      const cart_t& K, FloatType sweep_angle = 3) const
    {
#ifdef _DEBUG
      SMTBX_ASSERT(std::abs(K.length() - std::abs(K[2])) < 1e-6);
#endif
      std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K, sweep_angle);
      return alpha + (Sg - k.first) / k.second;
    }

    /* returns Sg or the given angle.
    This could be used to firther remove reflections from VF where
    excitation andle is too hight at the edge of the frame
    */
    FloatType angle_to_Sg(FloatType ang, const miller::index<>& h,
      const cart_t& K, FloatType sweep_angle)
    {
#ifdef _DEBUG
      SMTBX_ASSERT(std::abs(K.length() - std::abs(K[2])) < 1e-6);
#endif
      FloatType Sg1 = utils<FloatType>::calc_Sg(RMf * h, K);
      FloatType ang_diff = scitbx::deg_as_rad(sweep_angle);
      std::pair<mat3_t, cart_t> r = compute_RMf_N(alpha + ang_diff);
      FloatType Sg2 = utils<FloatType>::calc_Sg(r.first * h, K);
      FloatType k = (Sg2 - Sg1) / ang_diff;
      //FloatType a = Sg1 - k*alpha;
      //return k * ang + a = k*(ang-alpha) + Sg1;
      return Sg1 + k*(ang-alpha);
    }

    FloatType PL_correctionROD(const miller::index<>& h) const {
      return utils<FloatType>::PL_correctionROD(RMf * h);
    }

    bool is_excited_index(const miller::index<> &h,
      const cart_t& K,
      FloatType MaxSg,
      FloatType MaxG
      ) const;
  
    bool is_excited_beam(const BeamInfo<FloatType>& beam,
      const cart_t& K,
      FloatType MaxSg,
      FloatType MaxG
    ) const;

    void top_up(cart_t const& K,
      size_t num, FloatType min_d,
      FloatType MaxSg, FloatType MaxG,
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      sgtbx::space_group const& uniq_sg, bool force_frame_normal
    );

    void add_indices(const af::shared<miller::index<> >& indices);

    void add_beam(const miller::index<>& index,
      FloatType I, FloatType sig);
    // resets indices to match the beams
    void set_beams(af::shared<BeamInfo<FloatType> > const& beams);
    /* removes symmetry equivalents beams (and corresponding indices) and
    return the number of removed elements
    */
    size_t unify(sgtbx::space_group const& space_group,
      sgtbx::space_group const& uniq_sg, bool anomalous, bool exclude_sys_abs);

    void analyse_strength(
      af::shared<complex_t> const& Fcs_k,
      cart_t const& K,
      FloatType MaxG,
      FloatType MaxSg,
      FloatType MinP,
      bool force_Sg);
    /* computes integration angles making sure that the diffraction angles are
    in the list. Using threshold to merge near-by points
    * */
    af::shared<FloatType> get_int_angles(const cart_t& K, FloatType span,
      FloatType step, size_t N, bool use_Sg) const;
    /* returns all angles for span with step */
    static af::shared<FloatType> get_angles(FloatType  ang,
      FloatType span, FloatType step);
    af::shared<FloatType> get_angles_Sg(const miller::index<> &h,
      const cart_t& K, FloatType Sg_span, FloatType Sg_step) const;

    int id, tag;
    cart_t normal, original_normal;
    mat3_t UB, RMf;
    FloatType alpha, beta, omega, angle, scale;
    size_t offset; // for internal bookeeping
    // experimental data
    af::shared<BeamInfo<FloatType> > beams;
    // populated by analyse_strength
    af::shared<cart_t> gs;
    // all indices, first ones are strong_beams
    af::shared<miller::index<> > indices;
    af::shared<size_t> strong_beams, weak_beams,
      strong_measured_beams, weak_measured_beams;
    // 2 * Kl * Sg, populated by analyse_strength
    af::shared<FloatType> excitation_errors;
  };

  template <typename FloatType>
  struct BeamInfo {
    BeamInfo() {}

    BeamInfo(const miller::index<>& index,
      FloatType I, FloatType sig)
      : index(index),
      I(I), sig(sig)
    {}
    miller::index<> index;
    FloatType I, sig, diffraction_angle;
  };

  template <typename FloatType>
  bool FrameInfo<FloatType>::is_excited_index(const miller::index<>& h,
    const cart_t& K, FloatType MaxSg, FloatType MaxG) const
  {
    return utils<FloatType>::is_excited_h(h, RMf, K, MaxSg, MaxG, angle);
  }

  template <typename FloatType>
  bool FrameInfo<FloatType>::is_excited_beam(const BeamInfo<FloatType>& beam,
    const cart_t& K, FloatType MaxSg, FloatType MaxG) const
  {
    return is_excited_index(beam.index, K, MaxSg, MaxG);
  }

  template <typename FloatType>
  void FrameInfo<FloatType>::add_beam(
    const miller::index<>& index,
    FloatType I, FloatType sig)
  {
    beams.push_back(BeamInfo<FloatType>(index, I, sig));
    indices.push_back(index);
  }

  template <typename FloatType>
  void FrameInfo<FloatType>::set_beams(af::shared<BeamInfo<FloatType> > const& beams) {
    this->beams = beams.deep_copy();
    indices.clear();
    indices.reserve(beams.size());
    for (size_t i = 0; i < beams.size(); i++) {
      indices.push_back(beams[i].index);
    }
  }

  template <typename FloatType>
  void FrameInfo<FloatType>::add_indices(
    const af::shared<miller::index<> >& indices)
  {
    this->indices.reserve(this->indices.size() + indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
      this->indices.push_back(indices[i]);
    }
  }

  template <typename FloatType>
  void FrameInfo<FloatType>::top_up(
    cart_t const& K,
    size_t num, FloatType min_d,
    FloatType MaxSg, FloatType MaxG,
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const &space_group,
    sgtbx::space_group const& uniq_sg, bool force_frame_normal)
  {
    if (indices.size() >= num) {
      return;
    }
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;

    af::shared<typename utils<FloatType>::ExcitedBeam> ebeams;
    if (force_frame_normal) {
      ebeams = utils<FloatType>::generate_index_set(RMf, K, min_d,
        MaxG, MaxSg, unit_cell);
    }
    else {
      af::shared<mat3_t> RMfs(af::reserve(beams.size()));
      for (size_t i = 0; i < beams.size(); i++) {
        FloatType da = this->get_diffraction_angle(beams[i].index, K);
        RMfs.push_back(this->compute_RMf_N(da).first);
      }
      ebeams = utils<FloatType>::generate_index_set(RMfs, K, min_d,
        MaxG, MaxSg, unit_cell);
    }

    lookup_t existing = lookup_t(
      indices.const_ref(),
      uniq_sg,
      true);

    for (size_t i = 0; i < ebeams.size(); i++) {
      if (space_group.is_sys_absent(ebeams[i].h)) {
        continue;
      }
      if (existing.add_hkl(ebeams[i].h)) {
        indices.push_back(ebeams[i].h);
        if (indices.size() >= num) {
          break;
        }
      }
    }
  }

  // J.M. Zuo, A.L. Weickenmeier / Ultramicroscopy 57 (1995) 375-383
  template <typename FloatType>
  void FrameInfo<FloatType>::analyse_strength(
    af::shared<complex_t> const& Ugs,
    cart_t const& K,
    FloatType MaxG,
    FloatType MaxSg,
    FloatType MinP,
    bool force_Sg)
  {
    SMTBX_ASSERT(indices.size() == Ugs.size());
    strong_beams.clear();
    strong_measured_beams.clear();
    weak_beams.clear();
    weak_measured_beams.clear();
    excitation_errors.resize(indices.size());
    gs.resize(indices.size());

    FloatType Kl_sq = K.length_sq(),
      Kl = std::sqrt(Kl_sq);
    const FloatType max_f_sq = MaxG * MaxG,
      p_ang_sq = angle * angle;
    for (size_t i = 0; i < indices.size(); i++) {
      miller::index<> const& h = indices[i];
      cart_t g = RMf * cart_t(h[0], h[1], h[2]);
      gs[i] = g;
      //FloatType g_sq = g.length_sq();
      //if (g_sq > max_f_sq) { // || g_sq * p_ang_sq > 1.0) {
      //  continue;
      //}
      cart_t K_g = K + g;
      FloatType s = Kl_sq - K_g.length_sq();
      excitation_errors[i] = s;
      if (s < 0) {
        s = -s;
      }
      FloatType Sg = s / (2 * Kl);
      if (!force_Sg) {
        strong_beams.push_back(i);
        if (i < beams.size() && Sg < MaxSg) {
          strong_measured_beams.push_back(i);
        }
        continue;
      }
      if (Sg < MaxSg) {
        strong_beams.push_back(i);
        if (i < beams.size()) {
          strong_measured_beams.push_back(i);
        }
        continue;
      }
      FloatType strength_p = std::abs(Ugs[i] / s);
      if (strength_p > MinP) { // add for perturbation
        weak_beams.push_back(i);
        if (i < beams.size()) {
          weak_measured_beams.push_back(i);
        }
      }
    }
  }

  template <typename FloatType>
  size_t FrameInfo<FloatType>::unify(sgtbx::space_group const& space_group,
    sgtbx::space_group const& uniq_sg,
    bool anomalous, bool exclude_sys_abs)
  {
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;
    lookup_t bm = lookup_t(uniq_sg, anomalous);
    size_t cnt = 0;
    for (size_t i = 0; i < beams.size(); i++) {
      if ((exclude_sys_abs && space_group.is_sys_absent(beams[i].index)) ||
        !bm.add_hkl(beams[i].index))
      {
        beams.erase(&beams[i]);
        indices.erase(&indices[i]);
        cnt++;
        i--;
      }
    }
    return cnt;

  }
  template <typename FloatType>
  af::shared<FloatType> FrameInfo<FloatType>::get_angles_Sg(
    const miller::index<>& h,
    const cart_t& K, FloatType Sg_span, FloatType Sg_step) const
  {
    std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K);
    af::shared<FloatType> rv(af::reserve(std::abs(Sg_span * 2 / Sg_step) + 1));
    for (FloatType Sg = -Sg_span; Sg <= Sg_span; Sg += Sg_step) {
      FloatType ang = alpha + (Sg - k.first) / k.second;
      rv.push_back(ang);
    }
    return rv;
  }

  /* input span and step are in degrees */
  template <typename FloatType>
  af::shared<FloatType> FrameInfo<FloatType>::get_angles(FloatType ang,
    FloatType sn, FloatType st)
  {
    FloatType angle = scitbx::deg_as_rad(sn),
      step = scitbx::deg_as_rad(st);
    int steps = round(angle / step);
    af::shared<FloatType> rv(af::reserve(std::abs(steps * 2) + 1));
    for (int st = -steps; st <= steps; st++) {
      rv.push_back(ang + st * step);
    }
    return rv;
  }

  /* input span and step are in degrees */
  template <typename FloatType>
  af::shared<FloatType> FrameInfo<FloatType>::get_int_angles(
    const cart_t& K, FloatType span_, FloatType step_, size_t N, bool use_Sg) const
  {
    af::shared<FloatType> angles, d_angles, res;
    for (size_t i = 0; i < strong_measured_beams.size(); i++) {
      const miller::index<>& h = indices[strong_measured_beams[i]];
      FloatType da = this->get_diffraction_angle(h, K);
      d_angles.push_back(da);
      af::shared<FloatType> b_angles;
      if (use_Sg) {
        b_angles = get_angles_Sg(h, K, span_, step_);
      }
      else {
        b_angles = get_angles(alpha, span_, step_);
      }
      angles.extend(b_angles.begin(), b_angles.end());
    }
    if (angles.size() < 3) {
      return angles;
    }
    std::sort(angles.begin(), angles.end());
    std::sort(d_angles.begin(), d_angles.end());
    FloatType ang_span = (angles[angles.size() - 1] - angles[0]);
    FloatType threshold = ang_span / N;
    res.push_back(angles[0]);
    FloatType crv = angles[0];
    bool last_in = false;
    for (size_t i = 1; i < angles.size(); i++) {
      if (angles[i]-crv > threshold) {
        res.push_back(angles[i]);
        crv = angles[i];
        if (i == angles.size() - 1) {
          last_in = true;
        }
      }
    }
    if (!last_in) {
      res.push_back(angles[angles.size()-1]);
    }
    res.extend(d_angles.begin(), d_angles.end());

    std::sort(res.begin(), res.end());
    return res;
  }


  template <typename FloatType>
  struct PeakProfilePoint {
    FloatType I, Sg, angle, g;
    PeakProfilePoint()
    {}
    PeakProfilePoint(FloatType I, FloatType Sg, FloatType angle,
      FloatType g)
      : I(I), Sg(Sg), angle(angle), g(g)
    {}
  };

  template <typename FloatType>
  struct RefinementParams {
    af::shared<FloatType> values;
    RefinementParams(const af::shared<FloatType> &values)
      : values(values)
    {
      SMTBX_ASSERT(values.size() >= 14);
    }
    RefinementParams(const RefinementParams &params)
      : values(params.values)
    {}

    FloatType getKl_vac() const { return values[0]; }
    FloatType getKl() const { return values[1]; }
    FloatType getFc2Ug() const { return values[2]; }
    FloatType getEpsilon() const { return values[3]; }
    int getMatrixType() const { return static_cast<int>(values[4]); }
    void setMatrixType(int v) { values[4] = v; }
    int getBeamN() const { return static_cast<int>(values[5]); }
    int getThreadN() const { return static_cast<int>(values[6]); }
    FloatType getIntSpan() const { return values[7]; }
    FloatType getIntStep() const { return values[8]; }
    size_t getIntPoints() const { return static_cast<size_t>(values[9]); }
    bool isAngleInt() const { return values[10] == 1; }
    bool useNBeamSg() const { return values[11] == 1; }
    // with useNBeamSg - maxSg, otherwise is used as weight in |Fc|/(Sg+weight) 
    FloatType getNBeamWght() const { return values[12]; }
    bool isNBeamFloating() const { return values[13] != 0; }
  };

}}