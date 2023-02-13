#pragma once
#include <smtbx/ED/utils.h>
#include <scitbx/constants.h>

namespace smtbx { namespace ED {
using namespace cctbx;

template <typename FloatType> struct BeamInfo;

template <typename FloatType>
class FrameInfo {
public:
  typedef scitbx::vec3<FloatType> cart_t;
  typedef scitbx::mat3<FloatType> mat3_t;
  typedef std::complex<FloatType> complex_t;

  FrameInfo() {}

  FrameInfo(int id, const cart_t &f_normal,
    FloatType alpha, FloatType beta, FloatType omega,
    FloatType angle, FloatType scale, mat3_t const& UB)
    : id(id), tag(-1),
    alpha(alpha), beta(beta), omega(omega),
    angle(angle), scale(scale),
    offset(~0)
  {
    FloatType ca = std::cos(alpha), sa = std::sin(alpha),
      cb = std::cos(beta), sb = std::sin(beta),
      co = std::cos(omega), so = std::sin(omega);
    mat3_t rxa(1, 0, 0, 0, ca, -sa, 0, sa, ca),
      ryb(cb, 0, sb, 0, 1, 0, -sb, 0, cb),
      rzo(co, -so, 0, so, co, 0, 0, 0, 1);
    RM = rzo * rxa * ryb;
    RMf = RM * UB;
    normal = RM * f_normal;
    normal /= normal.length();
  }

  bool is_excited_index(const miller::index<> &h,
    FloatType Kl,
    FloatType MaxSg,
    FloatType MaxG
    ) const;
  
  bool is_excited_beam(const BeamInfo<FloatType>& beam,
    FloatType Kl,
    FloatType MaxSg,
    FloatType MaxG
  ) const;

  void top_up(FloatType Kl,
    size_t num, FloatType min_d,
    FloatType MaxSg, FloatType MaxG,
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const& space_group, bool anomalous
  );

  void add_indices(const af::shared<miller::index<> >& indices);

  void add_beam(const miller::index<>& index,
    FloatType I, FloatType sig);
  // resets indices to match the beams
  void set_beams(af::shared<BeamInfo<FloatType> > const& beams);
  /* removes symmetry equivalents beams (and corresponding indices) and
  return the number of removed elements
  */
  size_t unify(sgtbx::space_group const& space_group, bool anomalous);

  void analyse_strength(
    af::shared<complex_t> const& Fcs_k,
    cart_t const& K,
    FloatType MaxG,
    FloatType MaxSg,
    FloatType MinP);

  int id, tag;
  cart_t normal;
  mat3_t RM, RMf;
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
  FloatType I, sig;
};

template <typename FloatType>
bool FrameInfo<FloatType>::is_excited_index(const miller::index<>& h,
  FloatType Kl, FloatType MaxSg, FloatType MaxG) const
{
  return utils<FloatType>::is_excited_h(h, RMf, Kl, MaxSg, MaxG, angle);
}

template <typename FloatType>
bool FrameInfo<FloatType>::is_excited_beam(const BeamInfo<FloatType>& beam,
  FloatType Kl, FloatType MaxSg, FloatType MaxG) const
{
  return is_excited_index(beam.index, Kl, MaxSg, MaxG);
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
  FloatType Kl,
  size_t num, FloatType min_d,
  FloatType MaxSg, FloatType MaxG,
  uctbx::unit_cell const& unit_cell,
  sgtbx::space_group const &space_group, bool anomalous)
{
  if (indices.size() >= num) {
    return;
  }
  typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;

  std::vector<typename utils<FloatType>::ExcitedBeam> ebeams =
    utils<FloatType>::generate_index_set(RMf, Kl, min_d,
      MaxG, unit_cell);

  lookup_t existing = lookup_t(
    indices.const_ref(),
    space_group,
    anomalous);

  for (size_t i = 0; i < ebeams.size(); i++) {
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
  FloatType MinP)
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
    FloatType g_sq = g.length_sq();
    if (g_sq > max_f_sq) { // || g_sq * p_ang_sq > 1.0) {
      continue;
    }
    cart_t K_g = K + g;
    FloatType s = Kl_sq - K_g.length_sq();
    excitation_errors[i] = s;
    if (s < 0) {
      s = -s;
    }
    if (s / (2 * Kl) < MaxSg) {
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
size_t FrameInfo<FloatType>::unify(sgtbx::space_group const& space_group, bool anomalous)
{
  typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;
  lookup_t bm = lookup_t(space_group, anomalous);
  size_t cnt = 0;
  for (size_t i = 0; i < beams.size(); i++) {
    if (!bm.add_hkl(beams[i].index)) {
      beams.erase(&beams[i]);
      indices.erase(&indices[i]);
      cnt++;
      i--;
    }
  }
  return cnt;

}


}}