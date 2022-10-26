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
    normal = RMf * f_normal;
    normal /= normal.length();
  }
  bool is_excited(const BeamInfo<FloatType> &beam,
    FloatType Kl,
    FloatType MaxSg,
    FloatType MaxG
    ) const;
  
  void top_up(FloatType Kl,
    size_t num, FloatType min_d,
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const& space_group, bool anomalous
  );

  void add_beam(const miller::index<>& index,
    FloatType I, FloatType sig);

  int id, tag;
  cart_t normal;
  mat3_t RM, RMf;
  FloatType alpha, beta, omega, angle, scale;
  size_t offset; // for internal bookeeping
  // first indices correspond to measured beams
  af::shared<miller::index<> > indices;
  af::shared<BeamInfo<FloatType> > beams;
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
bool FrameInfo<FloatType>::is_excited(const BeamInfo<FloatType>& beam,
  FloatType Kl, FloatType MaxSg, FloatType MaxG) const
{
  return utils<FloatType>::is_excited_h(beam.index, RMf, Kl, MaxSg, MaxG, angle);
}

template <typename FloatType>
void FrameInfo<FloatType>::add_beam(
  const miller::index<>& index,
  FloatType I, FloatType sig)
{
  beams.push_back(BeamInfo<FloatType>(index, I, sig));
}

template <typename FloatType>
void FrameInfo<FloatType>::top_up(
  FloatType Kl,
  size_t num, FloatType min_d,
  uctbx::unit_cell const& unit_cell,
  sgtbx::space_group const &space_group, bool anomalous)
{
  if (indices.size() >= num) {
    return;
  }
  typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;

  af::shared<typename utils<FloatType>::ExcitedBeam> ebeams =
    utils<FloatType>::generate_index_set(RMf, Kl, min_d, unit_cell, space_group, anomalous);

  lookup_t existing = lookup_t(
    indices.const_ref(),
    space_group,
    anomalous);

  for (size_t i = 0; i < ebeams.size(); i++) {
    if (existing.find_hkl(ebeams[i].h) >= 0) {
      continue;
    }
    indices.push_back(ebeams[i].h);
    if (indices.size() >= num) {
      break;
    }
  }
}


}}