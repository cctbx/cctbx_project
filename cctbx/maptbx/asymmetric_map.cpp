#include <cctbx/maptbx/structure_factors.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <cctbx/maptbx/asymmetric_map.h>

namespace cctbx { namespace maptbx
{

scitbx::af::shared< std::complex<double > > asymmetric_map::structure_factors(
  scitbx::af::const_ref< cctbx::miller::index<> > indices) const
{
  typedef scitbx::af::c_grid_padded<3> grid_t;
  typedef std::complex<double> cmplx;
  scitbx::fftpack::real_to_complex_3d<double> fft(this->fft_grid().focus());
  auto fftmap = this->map_for_fft();
  CCTBX_ASSERT( fftmap.accessor().all().all_eq(fft.m_real()) );
  CCTBX_ASSERT( fftmap.accessor().focus().all_eq(fft.n_real()) );
  fft.forward(fftmap);
  grid_t cmplxpad(fft.n_complex(), fft.n_complex());
  scitbx::af::versa< cmplx, grid_t > cmap( fftmap.handle(),
    cmplxpad);
  bool anomalous_flag = false, conjugate_flag = true;
  cctbx::maptbx::structure_factors::from_map<double> frommap(
    this->space_group(),
    anomalous_flag,
    indices,
    cmap.const_ref(),
    conjugate_flag);
  scitbx::af::shared< cmplx > fcalc = frommap.data();
  double vol = 1. / fftmap.accessor().focus_size_1d();
  for(auto &fc : fcalc.ref())
    fc *= vol;
  CCTBX_ASSERT( indices.size() == fcalc.size() );
  return fcalc;
}

// prepare grid adapted integer symmetry operators
//! @todo code duplication, see atom_mask.h
std::vector<cctbx::sgtbx::grid_symop> asymmetric_map::grid_symops() const
{
  cctbx::sgtbx::space_group group = this->space_group();
  unsigned short order = group.order_z();
  CCTBX_ASSERT( order>0 );
  const scitbx::int3 n = this->unit_cell_grid_size();
  CCTBX_ASSERT( n[0]>0 && n[1]>0 && n[2] >0 );
  std::vector<cctbx::sgtbx::grid_symop> symops;
  symops.reserve(order);
  for(size_t i=0; i<order; ++i)
  {
    cctbx::sgtbx::grid_symop grsym( group(i), n );
    symops.push_back(grsym);
  }
  CCTBX_ASSERT( symops.size() == order );
  return symops;
}

asymmetric_map::fft_map_t asymmetric_map::map_for_fft() const
{
  scitbx::int3 map_size = this->unit_cell_grid_size();
  auto symops = this->grid_symops();
  unsigned short order = this->space_group().order_z();
  fft_map_t result(this->fft_grid(), 0.);
  for(grid_iterator_t i3=this->grid_begin(); !i3.over(); i3.incr() )
  {
    scitbx::int3 pos = i3();
    // if( has_enclosed_box )
    // {
    //   if( scitbx::ge_all(pos, emn) && scitbx::le_all(pos, emx) )
    //     nops = order;
    // }
    short ww = optimized_asu_.where_is(pos);
    if( ww!=0 )
    {
      scitbx::int3 pos_in_cell(pos);
      translate_into_cell(pos_in_cell, map_size);
      auto map_value = data_(pos);
      if( ww==-1 ) // inside on the face
      {
        unsigned short nops=site_symmetry_order(symops, pos_in_cell, map_size);
        CCTBX_ASSERT( nops>0 );
        CCTBX_ASSERT( order%nops == 0);
        map_value /= nops;
      }
      result(pos_in_cell) = map_value;
    }
    // if( nops!=0 )
    // {
    //   cell_volume += static_cast<size_t>(nops);
    // }
  }
  return result;
}

}}
