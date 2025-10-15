#pragma once
#include <smtbx/ED/utils.h>
#include <scitbx/constants.h>
#include <boost/shared_ptr.hpp>

namespace smtbx { namespace ED
{
  using namespace cctbx;

  template <typename FloatType> struct BeamInfo;
  template <typename FloatType> struct BeamGroup;
  template <typename FloatType>
  class EDData {
    ED_UTIL_TYPEDEFS;
    typedef typename utils<FloatType>::a_geometry geometry_t;
    struct Frame {
      cart_t normal;
      FloatType start, end, scale;
      Frame(const cart_t& n,
        FloatType start, FloatType end)
        : normal(n), start(start), end(end)
      {}
    };
    struct Beam {
      miller::index<> h;
      FloatType I, sig, scale, start, end;
      const Frame* frame;
      Beam(const Frame& f, const miller::index<>& h,
        FloatType I, FloatType sig, FloatType scale,
        FloatType start, FloatType end)
        : frame(&f),
        h(h), I(I), sig(sig), scale(scale),
        start(start), end(end)
      {}
    };
    af::shared<Beam> beams;
  public:
    EDData(boost::shared_ptr<geometry_t> geometry)
      : geometry(geometry)
    {}
    void add(const miller::index<>& h,
      FloatType start, FloatType middle, FloatType end)
    {
      beams.push_back(Beam(h, start, middle, end));
    }
    af::shared<BeamGroup<FloatType> > generate_groups(FloatType width) const;

    boost::shared_ptr<geometry_t> geometry;
  };

  template <typename FloatType>
  class BeamGroup {
  public:
    ED_UTIL_TYPEDEFS;
    typedef typename utils<FloatType>::a_geometry geometry_t;
    BeamGroup() {}

    BeamGroup(int id,
      boost::shared_ptr<geometry_t> geometry,
      FloatType angle, FloatType scale)
      : id(id), tag(-1), geometry(geometry),
      angle(angle), scale(scale),
      offset(~0), next(0), prev(0)
    {
      RMf = geometry->get_RMf(angle);
    }

    // convenience method
    mat3_t get_R(FloatType ang) const {
      return geometry->get_RMf(ang);
    }

    // convenience method
    const cart_t &get_N() const {
      return geometry->get_normal();
    }

    // returns angle in rads at which the excitation error is 0
    FloatType get_diffraction_angle(const miller::index<>& h,
      const cart_t& K) const
    {
      return geometry->get_diffraction_angle(h, K);
      //FloatType ang = Sg_to_angle(0, h, K);
      //mat3_t m = geometry->get_RMf(ang);
      //std::pair<FloatType, FloatType> k = Sg_to_angle_k(m, ang, h, K, 0.3);
      //ang = ang - k.first / k.second;
      //m = geometry->get_RMf(ang);
      //k = Sg_to_angle_k(m, ang, h, K, 0.1);
      //return ang - k.first / k.second;
    }

    /* returns angle in rads at which the excitation error is Sg as
    angle = rv.first + Sg/rv.second
    */
    std::pair<FloatType, FloatType> Sg_to_angle_k(
      const miller::index<>& h,
      const cart_t& K, FloatType sweep_angle = 0.5) const
    {
      FloatType da = geometry->get_diffraction_angle(h, K);
      FloatType ang_diff = scitbx::deg_as_rad(sweep_angle);
      mat3_t m = geometry->get_RMf(da + ang_diff);
      FloatType Sg2 = utils<FloatType>::calc_Sg(m * h, K);
      FloatType k = Sg2 / ang_diff;
      //FloatType a = Sg1 - k*angle;
      //return (Sg - a) / k = (Sg + k*angle - sg1)/k = angle + (Sg - Sg1)/k;
      return std::make_pair(da, k);
    }

    /* returns angle in rads at which the excitation error is Sg */
    FloatType Sg_to_angle(FloatType Sg, const miller::index<>& h,
      const cart_t& K) const
    {
      std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K);
      return k.first + Sg / k.second;
    }

    FloatType angle_to_Sg(FloatType ang, const miller::index<>& h,
      const cart_t& K, FloatType sweep_angle = 0.5) const
    {
      std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K, sweep_angle);
      return k.second*(ang - k.first);
    }

    FloatType calc_Sg(const miller::index<> &h,
      const cart_t& K) const
    {
      return utils<FloatType>::calc_Sg(
        RMf*cart_t(h[0], h[1], h[2]), K);
    }

    FloatType PL_correctionROD(const miller::index<>& h) const {
      return utils<FloatType>::PL_correctionROD(RMf * h);
    }

    bool is_excited_index(const miller::index<> &h,
      const cart_t& K,
      FloatType MaxSg,
      FloatType MaxG
      ) const;

    af::shared<miller::index<> > get_indices() const;

    bool contains_index(const miller::index<>& h) const;

    bool is_excited_beam(const BeamInfo<FloatType>& beam,
      const cart_t& K,
      FloatType MaxSg,
      FloatType MaxG
    ) const;

    void add_beam(const miller::index<>& index,
      FloatType I, FloatType sig);

    void set_beams(af::shared<BeamInfo<FloatType> > const& beams);
    /* removes symmetry equivalents beams and return the number of removed
    elements
    */
    size_t unify(sgtbx::space_group const& space_group,
      sgtbx::space_group const& uniq_sg, bool anomalous, bool exclude_sys_abs);

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
    af::shared<FloatType> get_angles_Sg_N(const miller::index<>& h,
      const cart_t& K, FloatType Sg_span, size_t N) const;
    af::shared<FloatType> get_angles_Sg_for_angles(const miller::index<>& h,
      const cart_t& K, FloatType start_ang, FloatType end_ang, size_t N) const;

    static bool sort_groups(const BeamGroup& a, const BeamGroup& b)
    {
      return (a.angle < b.angle);
    }

    static void link_groups(const af::shared<BeamGroup> &groups_) {
      af::shared<BeamGroup> groups = groups_;
      std::sort(groups.begin(), groups.end(), sort_groups);
      for (size_t i = 1; i < groups.size(); i++) {
        groups[i - 1].next = &groups[i];
        groups[i].prev = &groups[i-1];
      }
    }

    int id, tag;
    boost::shared_ptr<geometry_t> geometry;
    mat3_t RMf;
    FloatType angle, scale;
    size_t offset; // for internal bookeeping
    BeamGroup* next, * prev;

    // experimental data
    af::shared<BeamInfo<FloatType> > beams;
  };

  template <typename FloatType>
  struct BeamInfo {
    BeamInfo() {}

    BeamInfo(const miller::index<>& h,
      FloatType I, FloatType sig)
      : h(h),
      I(I), sig(sig)
    {}
    miller::index<> h;
    FloatType I, sig, diffraction_angle;
  };

  template <typename FloatType>
  bool BeamGroup<FloatType>::is_excited_index(const miller::index<>& h,
    const cart_t& K, FloatType MaxSg, FloatType MaxG) const
  {
    return utils<FloatType>::is_excited_h(h, RMf, K, MaxSg, MaxG, angle);
  }

  template <typename FloatType>
  bool BeamGroup<FloatType>::is_excited_beam(const BeamInfo<FloatType>& beam,
    const cart_t& K, FloatType MaxSg, FloatType MaxG) const
  {
    return is_excited_index(beam.h, K, MaxSg, MaxG);
  }

  template <typename FloatType>
  void BeamGroup<FloatType>::add_beam(
    const miller::index<>& h,
    FloatType I, FloatType sig)
  {
    beams.push_back(BeamInfo<FloatType>(h, I, sig));
  }

  template <typename FloatType>
  void BeamGroup<FloatType>::set_beams(af::shared<BeamInfo<FloatType> > const& beams) {
    this->beams = beams.deep_copy();
  }

  template <typename FloatType>
  size_t BeamGroup<FloatType>::unify(sgtbx::space_group const& space_group,
    sgtbx::space_group const& uniq_sg,
    bool anomalous, bool exclude_sys_abs)
  {
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;
    lookup_t bm = lookup_t(uniq_sg, anomalous);
    size_t cnt = 0;
    for (size_t i = 0; i < beams.size(); i++) {
      if ((exclude_sys_abs && space_group.is_sys_absent(beams[i].h)) ||
        !bm.add_hkl(beams[i].h))
      {
        beams.erase(&beams[i]);
        cnt++;
        i--;
      }
    }
    return cnt;

  }
  template <typename FloatType>
  af::shared<FloatType> BeamGroup<FloatType>::get_angles_Sg(
    const miller::index<>& h,
    const cart_t& K, FloatType Sg_span, FloatType Sg_step) const
  {
    std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K);
    af::shared<FloatType> rv(af::reserve(std::abs(Sg_span * 2 / Sg_step) + 1));
    for (FloatType Sg = -Sg_span; Sg <= Sg_span; Sg += Sg_step) {
      FloatType ang = k.first + Sg / k.second;
      rv.push_back(ang);
    }
    return rv;
  }

  template <typename FloatType>
  af::shared<FloatType> BeamGroup<FloatType>::get_angles_Sg_N(
    const miller::index<>& h,
    const cart_t& K, FloatType Sg_span, size_t N) const
  {
    std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K);
    af::shared<FloatType> rv(af::reserve(2*N+2));
    FloatType sg_step = Sg_span / (N - 1), ang;

    for (size_t i = N - 1; i >= 2; i--) {
      ang = k.first - sg_step*i / k.second;
      rv.push_back(ang);
    }
    for (int i = -2; i <= 2; i++) {
      ang = k.first + 0.5*sg_step * i / k.second;
      rv.push_back(ang);
    }
    for (size_t i = 2; i < N; i++) {
      ang = k.first + sg_step * i / k.second;
      rv.push_back(ang);
    }
    return rv;
  }

  template <typename FloatType>
  af::shared<FloatType> BeamGroup<FloatType>::get_angles_Sg_for_angles(
    const miller::index<>& h,
    const cart_t& K, FloatType start_ang, FloatType end_ang, size_t N) const
  {
    //std::pair<FloatType, FloatType> k = Sg_to_angle_k(h, K);
    //FloatType start_Sg = angle_to_Sg(start_ang, h, K),
    //  end_Sg = angle_to_Sg(end_ang, h, K);
    //if (start_Sg > end_Sg) {
    //  std::swap(start_Sg, end_Sg);
    //}
    //af::shared<FloatType> rv(af::reserve(N + 1));
    //FloatType Sg_step = (end_Sg - start_Sg) / N;
    //for (size_t i = 0; i <= N; i++) {
    //  FloatType ang = angle + (start_Sg + Sg_step*N - k.first) / k.second;
    //  rv.push_back(ang);
    //}

    //FloatType dang = get_diffraction_angle(h, K);
    //bool da_added = false;
    //af::shared<FloatType> rv(af::reserve(N + 1));
    //FloatType ang_step = (end_ang - start_ang) / N;
    //for (size_t i = 0; i <= N; i++) {
    //  FloatType ang = start_ang + ang_step*N;
    //  if (!da_added && ((ang_step > 0 && dang > ang) || (ang_step < 0 && dang < ang))) {
    //    rv.push_back(ang);
    //    da_added = true;
    //  }
    //  rv.push_back(ang);
    //}
    //return rv;

    FloatType start_Sg = angle_to_Sg(start_ang, h, K),
      end_Sg = angle_to_Sg(end_ang, h, K);
    FloatType Sg_span = std::max(std::abs(start_Sg), std::abs(end_Sg));
    FloatType Sg_step = Sg_span / N;
    return get_angles_Sg(h, K, Sg_span, Sg_step);

  }

  /* input span and step are in degrees */
  template <typename FloatType>
  af::shared<FloatType> BeamGroup<FloatType>::get_angles(FloatType ang,
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
  af::shared<FloatType> BeamGroup<FloatType>::get_int_angles(
    const cart_t& K, FloatType span_, FloatType step_, size_t N, bool use_Sg) const
  {
    af::shared<FloatType> angles, d_angles, res;
    for (size_t i = 0; i < beams.size(); i++) {
      const miller::index<>& h = beams[i].h;
      FloatType da = this->get_diffraction_angle(h, K);
      d_angles.push_back(da);
      af::shared<FloatType> b_angles;
      if (use_Sg) {
        b_angles = get_angles_Sg(h, K, span_, step_);
      }
      else {
        b_angles = get_angles(angle, span_, step_);
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
  af::shared<miller::index<> > BeamGroup<FloatType>::get_indices() const {
    af::shared<miller::index<> > rv(af::reserve(beams.size()));
    for (size_t bi = 0; bi < beams.size(); bi++) {
      rv.push_back(beams[bi].h);
    }
    return rv;
  }

  template <typename FloatType>
  bool BeamGroup<FloatType>::contains_index(const miller::index<>& h) const {
    for (size_t bi = 0; bi < beams.size(); bi++) {
      if (beams[bi].h == h) {
        return true;
      }
    }
    return false;
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
      SMTBX_ASSERT(values.size() > 20);
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
    /* The following three params are used in reflection profile 'detection'
    threshold is a fraction of the intensity to the intensity range to determine
    when the reflection 'starts'
    */
    FloatType getIntProfileStartTh() const { return values[14]; }
    FloatType getIntProfileSpan_Sg() const { return values[15]; }
    size_t getIntProfilePoints() const { return static_cast<size_t>(values[16]); }
    /* if scales ar not flat - neigbouring scales mixed through Sg to scale any
    particular intensity after the integration
    */
    bool useFlatScales() const { return static_cast<size_t>(values[17]); }
    FloatType getTopUpD() const { return values[18]; }
    FloatType getTopUpMaxSg() const { return values[19]; }
    FloatType getGroupWidth() const { return values[20]; }
  };

}}
