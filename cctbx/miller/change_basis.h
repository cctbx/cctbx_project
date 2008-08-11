#ifndef CCTBX_MILLER_CHANGE_BASIS_H
#define CCTBX_MILLER_CHANGE_BASIS_H

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace miller {

  template <typename FloatType>
  struct change_basis_phase_policy
  {
    static
    FloatType
    get(
      sym_equiv_index const& hr_ht,
      FloatType const& datum_in,
      bool deg)
    {
      return hr_ht.phase_eq(datum_in, deg);
    }
  };

  template <typename FloatType>
  struct change_basis_complex_policy
  {
    static
    std::complex<FloatType>
    get(
      sym_equiv_index const& hr_ht,
      std::complex<FloatType> const& datum_in,
      bool)
    {
      return hr_ht.complex_eq(datum_in);
    }
  };

  template <typename FloatType>
  struct change_basis_hendrickson_lattman_policy
  {
    static
    hendrickson_lattman<FloatType>
    get(
      sym_equiv_index const& hr_ht,
      hendrickson_lattman<FloatType> const& datum_in,
      bool)
    {
      return hr_ht.hendrickson_lattman_eq(datum_in);
    }
  };

  template <typename DataType, typename ChangeBasisPolicy>
  struct change_basis
  {
    af::shared<index<> > indices;
    af::shared<DataType> data;

    change_basis(
      sgtbx::change_of_basis_op const& cb_op,
      af::const_ref<index<> > const& indices_in,
      af::const_ref<DataType> const& data_in,
      bool deg=false)
    {
      CCTBX_ASSERT(data_in.size() == indices_in.size());
      indices.reserve(indices_in.size());
      data.reserve(data_in.size());
      sgtbx::tr_vec const& t = cb_op.c_inv().t();
      for(std::size_t i=0;i<indices_in.size();i++) {
        index<> const& h = indices_in[i];
        index<> hr = cb_op.apply(h);
        indices.push_back(hr);
        data.push_back(ChangeBasisPolicy::get(
          sym_equiv_index(hr, h*t, t.den(), false), data_in[i], deg));
      }
    }
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_CHANGE_BASIS_H
