#include <boost/timer/timer.hpp>

#include <cctbx/maptbx/structure_factors.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <cctbx/maptbx/asymmetric_map.h>
#include <iotbx/xplor/map_writer.h>

namespace cctbx { namespace maptbx
{

scitbx::af::shared< std::complex<double > > asymmetric_map::structure_factors(
  scitbx::af::const_ref< cctbx::miller::index<> > indices) const
{
  typedef scitbx::af::c_grid_padded<3> grid_t;
  typedef std::complex<double> cmplx;
  scitbx::fftpack::real_to_complex_3d<double> fft(this->fft_grid().focus());
  boost::timer::cpu_timer timer;
  auto fftmap = this->map_for_fft();
  this->map_for_fft_times_ = timer.format();
  CCTBX_ASSERT( fftmap.accessor().all().all_eq(fft.m_real()) );
  CCTBX_ASSERT( fftmap.accessor().focus().all_eq(fft.n_real()) );
  fft.forward(fftmap);
  grid_t cmplxpad(fft.n_complex(), fft.n_complex());
  scitbx::af::versa< cmplx, grid_t > cmap( fftmap.handle(), cmplxpad);
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
  scitbx::int3 ebox_min, ebox_max;
  bool has_enclosed_box=asu_.enclosed_box_corners(ebox_min, ebox_max, map_size);
  bool is_asu_inside_cell = this->box_begin().as_tiny().all_ge(0)
    && this->box_end().as_tiny().all_le(map_size+1);
  // CCTBX_ASSERT( is_asu_inside_cell ); No: 48 fails
  // CCTBX_ASSERT( has_enclosed_box ); No: 99
  scitbx::int3 padded_grid_size( result.accessor().all() );
  boost::timer::cpu_timer timer;
  if( has_enclosed_box && is_asu_inside_cell )
  {
    auto result_ref = result.ref();
    const auto *d = data_.begin();
    for(auto i3=this->mapped_begin(padded_grid_size); !i3.over(); i3.advance(),
      ++d)
    {
      scitbx::int3 pos = i3();
      auto t = pos.as_tiny();
      if( t.all_ge(ebox_min) && t.all_le(ebox_max) )
        result_ref[i3.mapped_index_1d()] = *d;
      else
      {
        short ww = optimized_asu_.where_is(pos);
        if( ww!=0 )
        {
          if( ww==-1 ) // inside on the face
          {
            scitbx::int3 pos_in_cell(pos);
            translate_into_cell(pos_in_cell, map_size);
            unsigned short ns=site_symmetry_order(symops,pos_in_cell,map_size);
            result_ref(pos_in_cell) = *d / ns;
          }
          else
            result_ref[i3.mapped_index_1d()] = *d;
        }
      }
    }
    CCTBX_ASSERT( d==data_.end() );
  }
  else
  {
    // slower generic version
    for(auto i3=this->grid_begin(); !i3.over(); i3.incr() )
    {
      scitbx::int3 pos = i3();
      short ww = optimized_asu_.where_is(pos);
      if( ww!=0 )
      {
        scitbx::int3 pos_in_cell(pos);
        translate_into_cell(pos_in_cell, map_size);
        auto map_value = data_(pos);
        if( ww==-1 ) // inside on the face
        {
          unsigned short ns=site_symmetry_order(symops, pos_in_cell, map_size);
          map_value /= ns;
        }
        result(pos_in_cell) = map_value;
      }
    }
  }
  this->fill_fft_times_ = timer.format();
  return result;
}

void asymmetric_map::copy_to_asu_box(const scitbx::int3 &map_size,
  const scitbx::int3 &padded_map_size, const double *cell_data)
{
  //! @todo: code duplication, see atom_mask
  //! @todo: test grid for compatiblity with the troup, see
  //     crystal_form::determin_grid
  scitbx::int3 grid = map_size;
  CCTBX_ASSERT( this->unit_cell_grid_size().as_tiny().all_eq(grid) );
  asu::rvector3_t box_min, box_max;
  asu_.box_corners(box_min, box_max);
  scitbx::int3 ebox_min, ebox_max;
  bool has_enclosed_box = asu_.enclosed_box_corners(ebox_min, ebox_max, grid);
  scitbx::mul(box_min, grid);
  scitbx::mul(box_max, grid);
  auto ibox_min = scitbx::floor(box_min);
  auto ibox_max = scitbx::ceil(box_max);
  if( has_enclosed_box )
  {
    CCTBX_ASSERT( scitbx::gt_all( ebox_min, ibox_min ) );
    CCTBX_ASSERT( scitbx::lt_all( ebox_max, ibox_max ) );
  }
  bool is_asu_inside_cell = ibox_min.as_tiny().all_ge(0)
    && ibox_max.as_tiny().all_le(grid);
  ibox_max += scitbx::int3(1,1,1); // open range: [box_min,box_max)
  data_.resize(asu_grid_t(ibox_min, ibox_max), 0.); //! @todo optimize
  boost::timer::cpu_timer timer;
  scitbx::int3 padded_grid_size = padded_map_size;
  CCTBX_ASSERT( padded_grid_size.as_tiny().all_ge(grid) );
  if( has_enclosed_box && is_asu_inside_cell )
  {
    auto *d = data_.begin();
    for(auto l3=this->mapped_begin(padded_grid_size); !l3.over(); l3.advance(),
      ++d)
    {
      scitbx::int3 pos = l3();
      if( scitbx::ge_all(pos,ebox_min) && scitbx::le_all(pos,ebox_max) )
          *d = cell_data[l3.mapped_index_1d()];
      else
      {
        short ww = optimized_asu_.where_is(pos);
        if( ww==-1 ) // inside on the face
        {
          scitbx::int3 pos_in_cell(pos);
          translate_into_cell(pos_in_cell, grid);
          *d = cell_data[ l3.mapped_index_1d(pos_in_cell) ];
        }
        else if( ww!=0 ) // inside
          *d = cell_data[l3.mapped_index_1d()];
      }
    }
    CCTBX_ASSERT( d==data_.end() );
  }
  else
  {
    // slower generic version
    auto *d = data_.begin();
    auto indexer = this->mapped_begin(padded_grid_size);
    for(grid_iterator_t l3=this->grid_begin(); !l3.over(); l3.incr(), ++d)
    {
      scitbx::int3 pos = l3();
      if(optimized_asu_.where_is(pos)!=0 )
      {
        scitbx::int3 pos_in_cell(pos);
        translate_into_cell(pos_in_cell, grid);
        *d = cell_data[ indexer.mapped_index_1d(pos_in_cell) ];
      }
    }
    CCTBX_ASSERT( d==data_.end() );
  }
  this->fill_density_times_ = timer.format();
}

void asymmetric_map::save(const std::string &file_name,
    const uctbx::unit_cell &unit_cell, format f) const
{
  if( f!=format::xplor )
    throw error("unsupported file format");
  af::const_ref<double, af::flex_grid<> > flex_data(data_.begin(),
    data_.accessor().as_flex_grid());
  iotbx::xplor::map_writer(file_name, unit_cell, flex_data,
    this->unit_cell_grid_size());
}

}}
