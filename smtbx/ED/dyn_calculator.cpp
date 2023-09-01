#include <smtbx/ED/dyn_calculator.h>
#include <smtbx/ED/dyn_calculator_def.h>
#include <smtbx/ED/dyn_calculator_2013.h>
#include <smtbx/ED/dyn_calculator_2015.h>

using namespace smtbx::ED;

template <typename FloatType>
smtbx::ED::dyn_calculator_factory<FloatType>::dyn_calculator_factory(int type)
  : type(type)
{
  SMTBX_ASSERT(type >0 && type <= DYN_CALCULATOR_LAST);
}

template <typename FloatType>
boost::shared_ptr<a_dyn_calculator<FloatType> > smtbx::ED::dyn_calculator_factory<FloatType>::make(
  const af::shared<miller::index<> >& indices,
  const cmat_t& mat_Ug,
  const cart_t& K,
  const mat3_t& RMf,
  const cart_t& N,
  FloatType thickness) const
{
  typedef boost::shared_ptr<a_dyn_calculator<FloatType> > ptr_t;
  switch (type) {
  case DYN_CALCULATOR_2013:
    return ptr_t(new dyn_calculator_2013<FloatType>(indices, mat_Ug, K, RMf, N, thickness));
  case DYN_CALCULATOR_2015:
    return ptr_t(new dyn_calculator_2015<FloatType>(indices, mat_Ug, K, RMf, N, thickness));
  }
  return ptr_t(new dyn_calculator_def<FloatType>(indices, mat_Ug, K, RMf, N, thickness));
}

template <typename FloatType>
boost::shared_ptr<a_dyn_calculator<FloatType> > smtbx::ED::dyn_calculator_factory<FloatType>::make(
  const af::shared<miller::index<> >& indices,
  const cart_t& K,
  FloatType thickness) const
{
  typedef boost::shared_ptr<a_dyn_calculator<FloatType> > ptr_t;
  switch (type) {
  case DYN_CALCULATOR_2013:
    return ptr_t(new dyn_calculator_2013<FloatType>(indices, K, thickness));
  case DYN_CALCULATOR_2015:
    return ptr_t(new dyn_calculator_2015<FloatType>(indices, K, thickness));
  }
  return ptr_t(new dyn_calculator_def<FloatType>(indices, K, thickness));
}

template class dyn_calculator_factory<double>;
