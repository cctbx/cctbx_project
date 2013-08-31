#ifndef CCTBX_MAPTBX_ASYMMETRIC_MAP_H
#define CCTBX_MAPTBX_ASYMMETRIC_MAP_H

#include <scitbx/array_family/accessors/c_interval_grid.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <cctbx/sgtbx/direct_space_asu/proto/asymmetric_unit.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/fftpack/real_to_complex_3d.h>
#include <mmtbx/masks/grid_symop.h>

namespace cctbx { namespace maptbx
{

//! @todo code duplication: mmtbx/masks/atom_mask.cpp
inline void translate_into_cell(scitbx::int3 &num, const scitbx::int3 &den)
{
  for(register unsigned char i=0; i<3; ++i)
  {
    register int tn = num[i];
    register const int td = den[i];
    tn %= td;
    if( tn < 0 )
      tn += td;
    num[i] = tn;
  }
}

//! @todo: code duplication, see atom_mask.h
inline unsigned short site_symmetry_order(
  const std::vector<cctbx::sgtbx::grid_symop> &symops,
  const scitbx::int3 &num, const scitbx::int3 &den )
{
  unsigned short nops = 0;
  // num must be  inside cell
  for(size_t i=0; i<symops.size(); ++i)
  {
    scitbx::int3 sv = symops[i].apply_to( num );
    translate_into_cell(sv, den);
    if( scitbx::eq_all(sv , num) )
      ++nops;
  }
  CCTBX_ASSERT( nops>0U );
  return nops;
}


namespace asu = cctbx::sgtbx::asu;

///// UNDER CONSTRUCTION
//  purporse: convert full unit_cell map into asymmetric unit sized map
//            convert processed aysmmetric map into form suitable for FFT
// motivation: performance improvment of various map algorithms
// code duplication: see atom_mask.h
class asymmetric_map
{
public:
  typedef scitbx::af::c_interval_grid<3> asu_grid_t;
  typedef scitbx::af::nested_loop<scitbx::int3> grid_iterator_t;
  typedef scitbx::af::c_grid_padded<3> fft_grid_t;
  typedef scitbx::af::versa<double, fft_grid_t> fft_map_t;

  // typedef scitbx::af::c_grid_periodic<3> unit_cell_grid_t;
  typedef scitbx::af::versa<double, asu_grid_t  > data_type;
  typedef scitbx::af::ref<double, asu_grid_t  >  data_ref_t;

  typedef scitbx::af::versa<double, fft_grid_t> unit_cell_map_t;

  const data_type &data() const { return data_; }

  data_ref_t data_ref() { return data_.ref(); }

  scitbx::int3 box_begin() const
  {
    return scitbx::int3(data_.accessor().origin());
  }

  scitbx::int3 box_end() const
  {
    return scitbx::int3(data_.accessor().last());
  }

  const asu::asymmetric_unit<asu::direct,asu::optimized> &optimized_asu() const
  {
    return optimized_asu_;
  }

  template<typename Grid>
  asymmetric_map(const sgtbx::space_group_type &group,
    scitbx::af::const_ref<double,Grid> cell_data) : asu_(group),
    optimized_asu_(asu_,cell_data.accessor().focus())
  {
    this->copy_to_asu_box(cell_data);
  }

  template<typename Grid>
  void copy_to_asu_box(scitbx::af::const_ref<double,Grid> cell_data)
  {
    scitbx::int3 grid( cell_data.accessor().focus() );
    CCTBX_ASSERT( this->unit_cell_grid_size().as_tiny().all_eq(grid) );
    asu::rvector3_t box_min, box_max;
    asu_.box_corners(box_min, box_max);
    scitbx::mul(box_min, grid);
    scitbx::mul(box_max, grid);
    auto ibox_min = scitbx::floor(box_min);
    auto ibox_max = scitbx::ceil(box_max);
    ibox_max += scitbx::int3(1,1,1); // open range: [box_min,box_max)
    // auto asu_size = ibox_max - ibox_min;
    data_.resize(asu_grid_t(ibox_min, ibox_max), 0.); //! @todo optimize
    for(grid_iterator_t loop3(ibox_min, ibox_max); !loop3.over(); loop3.incr() )
    {
      scitbx::int3 pos = loop3();
      if( optimized_asu_.where_is(pos)!=0 )
      {
        scitbx::int3 pos_in_cell(pos);
        translate_into_cell(pos_in_cell, grid);
        data_(pos) = cell_data(pos_in_cell);
      }
    }
  }

  const std::string &hall_symbol() const
  {
    return asu_.hall_symbol;
  }

  cctbx::sgtbx::space_group space_group() const
  {
    cctbx::sgtbx::space_group_symbols symbol("Hall: "+this->hall_symbol());
    cctbx::sgtbx::space_group group(symbol);
    CCTBX_ASSERT( group.type().hall_symbol() == this->hall_symbol() );
    return group;
  }


  scitbx::int3 unit_cell_grid_size() const
  {
    return optimized_asu_.grid_size();
  }

  fft_grid_t fft_grid() const
  {
    // @todo: size compatible with space group for fft ?
    // is it guaranteed by asymmetric_unit ?
    scitbx::af::tiny<std::size_t,3> n_real = this->unit_cell_grid_size();
    return fft_grid_t(scitbx::fftpack::m_real_from_n_real(n_real), n_real);
  }

  grid_iterator_t grid_begin() const
  {
    // open range : [begin, end)
    return grid_iterator_t(this->box_begin(), this->box_end());
  }

  //! Unit cell sized map suitable for fft
  fft_map_t map_for_fft() const;

  std::vector<cctbx::sgtbx::grid_symop> grid_symops() const;

  //! Structure factors from map
  scitbx::af::shared< std::complex<double> > structure_factors(
    scitbx::af::const_ref< cctbx::miller::index<> > indices) const;

  // do I really need this one ?
  unit_cell_map_t symmetry_expanded_unit_cell_map() const;

private:
  data_type data_;
  asu::direct_space_asu asu_;
  asu::asymmetric_unit<asu::direct, asu::optimized> optimized_asu_;
};

}} // cctbx::maptbx
#endif
