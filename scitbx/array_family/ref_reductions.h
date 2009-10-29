#ifndef SCITBX_ARRAY_FAMILY_REDUCTIONS_H
#define SCITBX_ARRAY_FAMILY_REDUCTIONS_H

#include <scitbx/error.h>
#include <scitbx/math/approx_equal.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/misc_functions.h>
#include <boost/optional.hpp>
#include <complex>

namespace scitbx { namespace af {

  template <typename ElementType1, typename AccessorType1,
            typename ElementType2, typename AccessorType2>
  int
  order(
    const_ref<ElementType1, AccessorType1> const& a1,
    const_ref<ElementType2, AccessorType2> const& a2)
  {
    std::size_t sz_min = (a1.size() < a2.size() ? a1.size() : a2.size());
    for(std::size_t i=0;i<sz_min;i++) {
      if (a1[i] < a2[i]) return -1;
      if (a2[i] < a1[i]) return 1;
    }
    if (a1.size() < a2.size()) return -1;
    if (a2.size() < a1.size()) return 1;
    return 0;
  }

  template <typename ElementType, typename AccessorType>
  std::size_t
  max_index(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("max_index() argument is an empty array");
    }
    std::size_t result = 0;
    for(std::size_t i=1;i<a.size();i++) {
      if (a[result] < a[i]) result = i;
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  std::size_t
  min_index(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("min_index() argument is an empty array");
    }
    std::size_t result = 0;
    for(std::size_t i=1;i<a.size();i++) {
      if (a[result] > a[i]) result = i;
    }
    return result;
  }

  template<typename ElementType, typename AccessorType, class PredicateType>
  boost::optional<std::size_t>
  first_index(const_ref<ElementType, AccessorType> const& a,
              PredicateType p)
  {
    typedef
      typename const_ref<ElementType, AccessorType>::const_iterator
      iter;
    boost::optional<std::size_t> result;
    iter i = std::find_if(a.begin(), a.end(), p);
    if (i != a.end()) result = i - a.begin();
    return result;
  }

  template<typename ElementType, typename AccessorType, class PredicateType>
  boost::optional<std::size_t>
  last_index(const_ref<ElementType, AccessorType> const& a,
             PredicateType p)
  {
    typedef
      typename const_ref<ElementType, AccessorType>::const_reverse_iterator
               iter;
    boost::optional<std::size_t> result;
    iter i = std::find_if(a.rbegin(), a.rend(), p);
    if (i != a.rend()) result = a.rend() - i - 1;
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  max(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("max() argument is an empty array");
    }
    ElementType result = a[0];
    for(std::size_t i=1;i<a.size();i++) {
      if (result < a[i]) result = a[i];
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  min(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("min() argument is an empty array");
    }
    ElementType result = a[0];
    for(std::size_t i=1;i<a.size();i++) {
      if (result > a[i]) result = a[i];
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  max_absolute(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("max_absolute() argument is an empty array");
    }
    ElementType result = fn::absolute(a[0]);
    for(std::size_t i=1;i<a.size();i++) {
      ElementType const& ai = a[i];
      if (ai > 0) {
        if (result <  ai) result =  ai;
      }
      else {
        if (result < -ai) result = -ai;
      }
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  sum(const_ref<ElementType, AccessorType> const& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  sum_sq(const_ref<ElementType, AccessorType> const& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += a[i] * a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  sum_sq(const_ref<std::complex<ElementType>, AccessorType> const& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += std::norm(a[i]);
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  norm(const_ref<ElementType, AccessorType> const& a)
  {
    return std::sqrt(sum_sq(a));
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  product(const_ref<ElementType, AccessorType> const& a)
  {
    std::size_t sz = a.size();
    if (sz == 0) return 0;
    ElementType result = 1;
    for(std::size_t i=0;i<sz;i++) result *= a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("mean() argument is an empty array");
    }
    ElementType result = a[0];
    for(std::size_t i=1;i<a.size();i++) result += a[i];
    return result * (1. / a.size());
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean_sq(const_ref<ElementType, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("mean_sq() argument is an empty array");
    }
    ElementType result = a[0] * a[0];
    for(std::size_t i=1;i<a.size();i++) result += a[i] * a[i];
    return result * (1. / a.size());
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean_sq(const_ref<std::complex<ElementType>, AccessorType> const& a)
  {
    if (a.size() == 0) {
      throw std::runtime_error("mean_sq() argument is an empty array");
    }
    ElementType result = std::norm(a[0]);
    for(std::size_t i=1;i<a.size();i++) result += std::norm(a[i]);
    return result / a.size();
  }

  template <typename ElementTypeValues, typename AccessorTypeValues,
            typename ElementTypeWeights, typename AccessorTypeWeights>
  ElementTypeValues
  mean_weighted(
    const_ref<ElementTypeValues, AccessorTypeValues> const& values,
    const_ref<ElementTypeWeights, AccessorTypeWeights> const& weights)
  {
    if (values.size() != weights.size()) throw_range_error();
    if (values.size() == 0) {
      throw std::runtime_error("mean_weighted() argument is an empty array");
    }
    ElementTypeValues sum_vw = values[0] * weights[0];
    ElementTypeWeights sum_w = weights[0];
    for(std::size_t i=1;i<values.size();i++) {
      sum_vw += values[i] * weights[i];
      sum_w += weights[i];
    }
    return sum_vw / sum_w;
  }

  template <typename ElementTypeValues, typename AccessorTypeValues,
            typename ElementTypeWeights, typename AccessorTypeWeights>
  ElementTypeValues
  mean_sq_weighted(
    const_ref<ElementTypeValues, AccessorTypeValues> const& values,
    const_ref<ElementTypeWeights, AccessorTypeWeights> const& weights)
  {
    if (values.size() != weights.size()) throw_range_error();
    if (values.size() == 0) {
      throw std::runtime_error(
        "mean_sq_weighted() argument is an empty array");
    }
    ElementTypeValues sum_vvw = values[0] * values[0] * weights[0];
    ElementTypeWeights sum_w = weights[0];
    for(std::size_t i=1;i<values.size();i++) {
      sum_vvw += values[i] * values[i] * weights[i];
      sum_w += weights[i];
    }
    return sum_vvw / sum_w;
  }

  template <typename ElementType>
  struct min_max_mean
  {
    min_max_mean() {}

    template <typename AccessorType>
    min_max_mean(const_ref<ElementType, AccessorType> const& values)
    :
      n(values.size())
    {
      if (n == 0) return;
      ElementType min_v = values[0];
      ElementType max_v = min_v;
      ElementType sum_v = min_v;
      for(std::size_t i=1;i<values.size();i++) {
        ElementType const& vi = values[i];
        if (vi < min_v) min_v = vi;
        if (max_v < vi) max_v = vi;
        sum_v += vi;
      }
      min = min_v;
      max = max_v;
      sum = sum_v;
      mean = sum_v / values.size();
    }

    std::size_t n;
    boost::optional<ElementType> min;
    boost::optional<ElementType> max;
    boost::optional<ElementType> sum;
    boost::optional<ElementType> mean;
  };

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_eq(const_ref const& other) const
  {
    const ElementType* o = other.begin();
    const ElementType* t = begin();
    const ElementType* e = end();
    if (e-t != other.end()-o) return false;
    while (t != e) {
      if (!(*t++ == *o++)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_eq(ElementType const& other) const
  {
    const ElementType* t = begin();
    const ElementType* e = end();
    while (t != e) {
      if (!(*t++ == other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_ne(const_ref const& other) const
  {
    const ElementType* o = other.begin();
    const ElementType* t = begin();
    const ElementType* e = end();
    if (e-t != other.end()-o) return false;
    while (t != e) {
      if (!(*t++ != *o++)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_ne(ElementType const& other) const
  {
    const ElementType* t = begin();
    const ElementType* e = end();
    while (t != e) {
      if (!(*t++ != other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_lt(const_ref const& other) const
  {
    const ElementType* o = other.begin();
    const ElementType* t = begin();
    const ElementType* e = end();
    if (e-t != other.end()-o) throw_range_error();
    while (t != e) {
      if (!(*t++ < *o++)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_lt(ElementType const& other) const
  {
    const ElementType* t = begin();
    const ElementType* e = end();
    while (t != e) {
      if (!(*t++ < other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_gt(const_ref const& other) const
  {
    const ElementType* o = other.begin();
    const ElementType* t = begin();
    const ElementType* e = end();
    if (e-t != other.end()-o) throw_range_error();
    while (t != e) {
      if (!(*t++ > *o++)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_gt(ElementType const& other) const
  {
    const ElementType* t = begin();
    const ElementType* e = end();
    while (t != e) {
      if (!(*t++ > other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_le(const_ref const& other) const
  {
    const ElementType* o = other.begin();
    const ElementType* t = begin();
    const ElementType* e = end();
    if (e-t != other.end()-o) throw_range_error();
    while (t != e) {
      if (!(*t++ <= *o++)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_le(ElementType const& other) const
  {
    const ElementType* t = begin();
    const ElementType* e = end();
    while (t != e) {
      if (!(*t++ <= other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_ge(const_ref const& other) const
  {
    const ElementType* o = other.begin();
    const ElementType* t = begin();
    const ElementType* e = end();
    if (e-t != other.end()-o) throw_range_error();
    while (t != e) {
      if (!(*t++ >= *o++)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_ge(ElementType const& other) const
  {
    const ElementType* t = begin();
    const ElementType* e = end();
    while (t != e) {
      if (!(*t++ >= other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_approx_equal(
    const_ref const& other,
    typename scitbx::math::abs_traits<
      ElementType>::result_type const& tolerance) const
  {
    if (size() != other.size()) return false;
    scitbx::math::approx_equal_absolutely<ElementType> predicate(tolerance);
    for (std::size_t i=0; i<size(); ++i) {
      if (!predicate((*this)[i], other[i])) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_approx_equal(
    ElementType const& other,
    typename scitbx::math::abs_traits<
      ElementType>::result_type const& tolerance) const
  {
    scitbx::math::approx_equal_absolutely<ElementType> predicate(tolerance);
    for (std::size_t i=0; i<size(); ++i) {
      if (!predicate((*this)[i], other)) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_approx_equal_relatively(
    const_ref const& other,
    typename scitbx::math::abs_traits<
      ElementType>::result_type const& relative_error) const
  {
    if (size() != other.size()) return false;
    scitbx::math::approx_equal_relatively<ElementType> predicate(relative_error);
    for (std::size_t i=0; i<size(); ++i) {
      if (!predicate((*this)[i], other[i])) return false;
    }
    return true;
  }

  template <typename ElementType, typename AccessorType>
  bool
  const_ref<ElementType, AccessorType>
  ::all_approx_equal_relatively(
    ElementType const& other,
    typename scitbx::math::abs_traits<
      ElementType>::result_type const& relative_error) const
  {
    scitbx::math::approx_equal_relatively<ElementType> predicate(relative_error);
    for (std::size_t i=0; i<size(); ++i) {
      if (!predicate((*this)[i], other)) return false;
    }
    return true;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_REDUCTIONS_H
