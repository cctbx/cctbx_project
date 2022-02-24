#pragma once
#include <boost/shared_ptr.hpp>
#include <cctbx/miller.h>
#include <scitbx/vec3.h>

namespace cctbx { namespace smtbx { namespace ED {

template <typename FloatType>
struct FrameInfo {
  typedef scitbx::vec3<FloatType> cart_t;

  FrameInfo(const cart_t &normal,
    FloatType alpha, FloatType beta, FloatType omega,
    FloatType angle, FloatType scale)
    : normal(normal),
    alpha(alpha), beta(beta), omega(omega),
    angle(angle), scale(scale)
  {}

  cart_t normal;
  FloatType alpha, beta, omega, angle, scale;
};

template <typename FloatType>
struct BeamInfo {
  typedef FrameInfo<FloatType> parent_t;
  typedef boost::shared_ptr<parent_t> parent_ptr_t;

  BeamInfo(parent_ptr_t parent, const miller::index<>& index,
    FloatType I, FloatType sig)
    : parent(parent),
    index(index),
    I(I), sig(sig)
  {}
  
  boost::shared_ptr<FrameInfo<FloatType> > parent;
  miller::index<> index;
  FloatType I, sig;
};
}}}