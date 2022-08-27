#pragma once
#include <boost/shared_ptr.hpp>
#include <cctbx/miller.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <smtbx/error.h>

namespace smtbx { namespace ED {
using namespace cctbx;

template <typename FloatType> struct BeamInfo;

template <typename FloatType>
struct FrameInfo {
  typedef scitbx::vec3<FloatType> cart_t;
  typedef scitbx::mat3<FloatType> mat3_t;

  FrameInfo() {}

  FrameInfo(int id, const cart_t &f_normal,
    FloatType alpha, FloatType beta, FloatType omega,
    FloatType angle, FloatType scale, mat3_t const& UB)
    : id(id), tag(-1),
    alpha(alpha), beta(beta), omega(omega),
    angle(angle), scale(scale)
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
    FloatType MaxSg
    ) const;
  
  int id, tag;
  cart_t normal;
  mat3_t RM, RMf;
  FloatType alpha, beta, omega, angle, scale;
};

template <typename FloatType>
struct BeamInfo {
  typedef FrameInfo<FloatType> parent_t;

  BeamInfo() {}

  BeamInfo(parent_t *parent, const miller::index<>& index,
    FloatType I, FloatType sig)
    : parent(parent),
    index(index),
    I(I), sig(sig)
  {
    SMTBX_ASSERT(parent != 0);
  }
  
  FrameInfo <FloatType>& get_parent() { return *parent; }
  FrameInfo<FloatType>* parent;
  miller::index<> index;
  FloatType I, sig;
};
template <typename FloatType>
bool FrameInfo<FloatType>::is_excited(const BeamInfo<FloatType>& beam,
  FloatType Kl, FloatType MaxSg) const
{
  typedef scitbx::vec3<FloatType> cart_t;
  cart_t g = RMf * cart_t(beam.index[0], beam.index[1], beam.index[2]);
  g[2] += Kl;
  FloatType Sg = (Kl * Kl - g.length_sq()) / (2 * Kl);
  return std::abs(Sg) < MaxSg;
}

}}