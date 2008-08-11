#ifndef CCTBX_MAPTBX_ABSTRACT_MAP_H
#define CCTBX_MAPTBX_ABSTRACT_MAP_H

#include<scitbx/array_family/flex_types.h>
#include<cctbx/uctbx.h>
#include<cctbx/maptbx/coordinate_transformers.h>
#include<cctbx/maptbx/mapper.h>
#include<cctbx/maptbx/generic_grid.h>
#include<cctbx/maptbx/map_out_of_bounds.h>
#include<vector>
#include<cmath>

#include<iostream>

namespace cctbx {

namespace maptbx {

template < typename FloatType, typename IntType >
class basic_map {
public:
  static const std::size_t                              corners_8 = 8;
  typedef generic_grid<void,FloatType,IntType>          grid_intf_type;
  typedef chiltbx::handle::handle<grid_intf_type>       grid_handle_type;
  typedef fractional<FloatType>                         frac_type;
  typedef cartesian<FloatType>                          cart_type;
  typedef grid_point<IntType>                           grid_type;
  typedef scitbx::vec3<FloatType>                       vec3_type;
  typedef af::tiny<IntType,dimension_3>                 tin3_type;
  typedef scitbx::mat3<FloatType>                       mat3_type;
  typedef maptbx::transform<frac_type,grid_type>        f2g_type;
  typedef maptbx::transform<frac_type,cart_type>        f2c_type;
  typedef maptbx::transform<grid_type,frac_type>        g2f_type;
  typedef maptbx::transform<cart_type,frac_type>        c2f_type;
  typedef maptbx::transform<cart_type,grid_type>        c2g_type;
  typedef maptbx::transform<grid_type,cart_type>        g2c_type;
  typedef out_of_bounds<void,FloatType,IntType>         out_of_bounds_type;
  typedef chiltbx::handle::handle<out_of_bounds_type>   out_of_handle_type;

  typedef std::vector<grid_type>                        grid_list_type;
  typedef std::vector<frac_type>                        frac_list_type;
  typedef af::tiny<FloatType,2>                         weight_pair_type;
  typedef af::tiny<weight_pair_type,3>                  weight_pairs_type;
  typedef maptbx::basic_mapper<void,FloatType,IntType>  basic_mapper;
  typedef af::tiny<grid_type,2>                         grid_pair_type;

  typedef af::ref<FloatType, af::flex_grid<> >          flex_grid_ref;
  typedef af::const_ref<FloatType, af::flex_grid<> >    const_flex_grid_ref;
  typedef af::c_grid_padded<dimension_3>                c_grid_3;
  typedef af::ref<FloatType, c_grid_3>                  c_grid_ref;
  typedef af::const_ref<FloatType, c_grid_3 >           const_c_grid_ref;

  typedef cctbx::uctbx::unit_cell                       tbx_unit_cell;

  typedef cctbx::maptbx::unit_cell                      unit_cell_tag;

  basic_map (
    asu const&,
    af::versa<FloatType, af::flex_grid<> > const& data,
    tbx_space_group const& space_group,
    cdsa::float_asu<FloatType> const& fasu,
    tin3_type const& extents,
    mat3_type const& matrix,
    out_of_handle_type const& h_out_of,
    tbx_unit_cell const& unitcell,
    FloatType const& min_distance_sym_equiv=0.5,
    bool assert_min_distance_sym_equiv=true  )
  : out_of_handle_(h_out_of) {
    this->set_grid_handle( data,
                           space_group,
                           fasu,
                           extents,
                           min_distance_sym_equiv,
                           assert_min_distance_sym_equiv);
    this->rebuild_transformers(extents,matrix);
    this->unit_cell_ = unitcell;
  }
  basic_map (
    unit_cell_tag const&,
    af::versa<FloatType, af::flex_grid<> > const& data,
    tin3_type const& extents,
    mat3_type const& matrix,
    out_of_handle_type const& h_out_of,
    tbx_unit_cell const& unitcell )
  : out_of_handle_(h_out_of) {
    this->set_grid_handle(data);
    this->rebuild_transformers(extents,matrix);
    this->unit_cell_ = unitcell;
  }
  basic_map (
    unit_cell_tag const&,
    af::versa<FloatType, af::c_grid_padded<dimension_3> > const& data,
    tin3_type const& extents,
    mat3_type const& matrix,
    out_of_handle_type const& h_out_of,
    tbx_unit_cell const& unitcell )
  : out_of_handle_(h_out_of) {
    this->set_grid_handle(data);
    this->rebuild_transformers(extents,matrix);
    this->unit_cell_ = unitcell;
  }
  basic_map (
    non_symmetric const&,
    af::versa<FloatType, af::flex_grid<> > const& data,
    tin3_type const& extents,
    mat3_type const& matrix,
    out_of_handle_type const& h_out_of,
    tbx_unit_cell const& unitcell )
  : out_of_handle_(h_out_of) {
    this->set_grid_handle(data,extents);
    this->rebuild_transformers(extents,matrix);
    this->unit_cell_ = unitcell;
  }
  basic_map ( basic_map const& bm ) {
    this->grid_handle_    = bm.grid_handle_;
    this->out_of_handle_  = bm.out_of_handle_;
    this->rebuild_transformers(bm.extents_,bm.matrix_);
    this->unit_cell_      = bm.unit_cell_;
  }
  basic_map& operator = ( basic_map const& bm ) {
    if ( this == &bm )
      return *this;
    this->grid_handle_    = bm.grid_handle_;
    this->out_of_handle_  = bm.out_of_handle_;
    this->rebuild_transformers(bm.extents_,bm.matrix_);
    this->unit_cell_      = bm.unit_cell_;
    return *this;
  }
  void rebuild_transformers ( tin3_type const& extents, mat3_type const& matrix ) {
    this->extents_    = extents;
    this->matrix_    = matrix;
    this->frac2grid_  = f2g_type(this->extents_);
    this->grid2frac_  = this->frac2grid_.inverse();
    this->frac2cart_  = f2c_type(this->matrix_);
    this->cart2frac_  = this->frac2cart_.inverse();
    this->cart2grid_  = c2g_type(this->cart2frac_,this->frac2grid_);
    this->grid2cart_  = g2c_type(this->grid2frac_,this->frac2cart_);
  }
  af::tiny<IntType,dimension_3> origin () const {
    return this->grid_handle_->origin();
  }
  af::tiny<IntType,dimension_3> last () const {
    return this->grid_handle_->last();
  }
  af::tiny<IntType,dimension_3> focus () const {
    return this->grid_handle_->focus();
  }
  af::tiny<IntType,dimension_3> all () const {
    return this->grid_handle_->all();
  }
  FloatType grid_value ( grid_point<IntType> const& coordinate ) const {
    // does not remap, check, etc. etc.
    return this->grid_handle_->get_value(coordinate);
  }
  void set_grid_value ( grid_point<IntType> const& coordinate, FloatType const& F ) {
    this->grid_handle_->set_value(coordinate,F);
  }
  FloatType get_cart_value ( cart_type const& coordinate ) const {
    return this->base_get_value( this->cart2frac_(coordinate) );
  }
  af::shared<FloatType>
  get_cart_values ( af::const_ref<scitbx::vec3<FloatType> > const& coordinates ) const {
    af::shared<FloatType> result(coordinates.size());
    for ( std::size_t i=0; i<coordinates.size(); ++i )
      result[i] = this->get_cart_value(coordinates[i]);
    return result;
  }
  FloatType get_frac_value ( frac_type const& coordinate ) const {
    return this->base_get_value(coordinate);
  }
  af::shared<FloatType>
  get_frac_values ( af::const_ref<scitbx::vec3<FloatType> > const& coordinates ) const {
    af::shared<FloatType> result(coordinates.size());
    for ( std::size_t i=0; i<coordinates.size(); ++i )
      result[i] = this->get_frac_value(coordinates[i]);
    return result;
  }
  FloatType get_grid_value ( grid_type const& coordinate ) const {
    grid_type gcoord = this->remap_coordinate(coordinate);
    if ( this->grid_handle_->is_inside(gcoord) )
      return this->grid_handle_->get_value(gcoord);
    else
      try {
        return this->out_of_handle_->retry(
            *this->grid_handle_.const_get_raw(),
            this->grid2frac_,
            this->grid2frac_(gcoord));
      } catch ( FloatType const& result ) {
        return result;
      }
  }
  af::shared<FloatType>
  get_grid_values ( af::const_ref<scitbx::vec3<IntType> > const& coordinates ) const {
    af::shared<FloatType> result(coordinates.size());
    for ( std::size_t i=0; i<result.size(); ++i )
      result[i] = this->get_grid_value(coordinates[i]);
    return result;
  }
  FloatType get_value ( cart_type const& coordinate ) const {
    return this->get_cart_value(coordinate);
  }
  FloatType get_value ( frac_type const& coordinate ) const {
    return this->get_frac_value(coordinate);
  }
  FloatType get_value ( grid_type const& coordinate ) const {
    return this->get_grid_value(coordinate);
  }
  template < typename CoordType >
  FloatType operator [] ( CoordType const& coordinate ) const {
    return this->get_value(coordinate);
  }
  template < typename CoordType >
  af::shared<FloatType> operator [] ( af::const_ref<CoordType> const& coordinates ) const {
    af::shared<FloatType> result(coordinates.size());
    for ( std::size_t i=0; i<result.size(); ++i )
      result[i] = (*this)[coordinates[i]];
    return result;
  }
  frac_type remap_frac_coordinate ( frac_type const& coordinate ) const {
    frac_type result = this->grid_handle_->remap(coordinate).mapped_coordinate;
    return result;
  }
  cart_type remap_cart_coordinate ( cart_type const& coordinate ) const {
    frac_type fcoord = this->cart2frac_(coordinate);
    fcoord = this->grid_handle_->remap(fcoord).mapped_coordinate;
    return this->frac2cart_(fcoord);
  }
  grid_type remap_grid_coordinate ( grid_type const& coordinate ) const {
    frac_type fcoord = this->grid2frac_(coordinate);
    fcoord = this->grid_handle_->remap(fcoord).mapped_coordinate;
    return this->frac2grid_(fcoord);
  }
  frac_type remap_coordinate ( frac_type const& coordinate ) const {
    return this->remap_frac_coordinate(coordinate);
  }
  cart_type remap_coordinate ( cart_type const& coordinate ) const {
    return this->remap_cart_coordinate(coordinate);
  }
  grid_type remap_coordinate ( grid_type const& coordinate ) const {
    return this->remap_grid_coordinate(coordinate);
  }
  grid_type frac_nearest_grid_point ( frac_type const& coordinate ) const {
    return this->remap_coordinate(this->frac2grid_(coordinate));
  }
  grid_type cart_nearest_grid_point ( cart_type const& coordinate ) const {
    return this->remap_coordinate(this->cart2grid_(coordinate));
  }
  grid_type nearest_grid_point ( frac_type const& coordinate ) const {
    return this->frac_nearest_grid_point(coordinate);
  }
  grid_type nearest_grid_point ( cart_type const& coordinate ) const {
    return this->cart_nearest_grid_point(coordinate);
  }
  grid_type nearest_grid_point ( grid_type const& coordinate ) const {
    return this->remap_coordinate(coordinate);
  }
  FloatType frac_value_at_nearest_grid_point ( frac_type const& coordinate ) const {
    return this->value_at_nearest_grid_point(coordinate);
  }
  FloatType cart_value_at_nearest_grid_point ( cart_type const& coordinate ) const {
    return this->value_at_nearest_grid_point(coordinate);
  }
  template < typename CoordType >
  FloatType value_at_nearest_grid_point ( CoordType const& coordinate ) const {
    return this->grid_handle_->get_value( this->nearest_grid_point(coordinate) );
  }
  g2f_type grid_to_frac () const {
    return this->grid2frac_;
  }
  c2f_type cart_to_frac () const {
    return this->cart2frac_;
  }
  f2g_type frac_to_grid () const {
    return this->frac2grid_;
  }
  f2c_type frac_to_cart () const {
    return this->frac2cart_;
  }
  g2c_type grid_to_cart () const {
    return this->grid2cart_;
  }
  c2g_type cart_to_grid () const {
    return this->cart2grid_;
  }
  basic_mapper remap ( frac_type const& coordinate ) const {
    return this->grid_handle_->remap(coordinate);
  }
  bool is_inside ( grid_point<IntType> const& coordinate ) const {
    return this->grid_handle_->is_inside(coordinate);
  }
  bool is_inside ( fractional<FloatType> const& coordinate ) const {
    return this->grid_handle_->is_inside(coordinate);
  }
  bool is_inside ( cartesian<FloatType> const& coordinate ) const {
    return this->grid_handle_->is_inside( this->cart2frac_(coordinate) );
  }
  grid_type get_corner000 ( frac_type const& coordinate ) const {
    return this->frac2grid_.floor_transform(coordinate);
  }
  grid_list_type get_corners ( grid_type const& corner000 ) const {
    grid_list_type corners(corners_8,corner000);
    std::size_t index = 0;
    for ( std::size_t i=0; i<2; ++i )
      for ( std::size_t j=0; j<2; ++j )
        for ( std::size_t k=0; k<2; ++k, ++index ) {
          corners[index][0] += (1==i?0:1);
          corners[index][1] += (1==j?0:1);
          corners[index][2] += (1==k?0:1);
        }
    return corners;
  }
  weight_pairs_type
  get_weights ( grid_type const& corner000, frac_type const& coordinate ) const {

    scitbx::vec3<FloatType> g_coordinate = this->frac2grid_.strange_transform(coordinate);

    //weight_pairs_type result(3,weight_pair_type(2,0.0));
    weight_pairs_type result(weight_pair_type(0.0,0.0),
                             weight_pair_type(0.0,0.0),
                             weight_pair_type(0.0,0.0));

    for ( std::size_t i=0; i<dimension_3; ++i ) {
      result[i][0] = g_coordinate[i] - corner000[i];
      result[i][1] = 1 - result[i][0];
    }

    return result;
  }
  FloatType base_get_value ( frac_type const& coordinate ) const {

    // get the corners and map into fractional
    grid_type corner000 = this->get_corner000(coordinate);
    grid_list_type corners = this->get_corners(corner000);
    // generate the linear-alpha-weighting
    weight_pairs_type weights = this->get_weights(corner000,coordinate);

    CCTBX_ASSERT(corners.size()==corners_8);

    // get_grid_value remaps the coordinate
    // Build the interpolated value over the linear-alpha-weights.
    FloatType result = 0;
    std::size_t index = 0;
    for ( std::size_t i=0; i<2; ++i )
      for ( std::size_t j=0; j<2; ++j )
        for ( std::size_t k=0; k<2; ++k, ++index ) {
          result += this->get_grid_value(corners[index])
                     *weights[0][i]*weights[1][j]*weights[2][k];
        }

    return result;
  }
  void set_grid_handle (
    af::versa<FloatType, af::flex_grid<> > const& data,
    tbx_space_group const& space_group,
    cdsa::float_asu<FloatType> const& fasu,
    af::tiny<IntType,dimension_3> const& grid_len,
    FloatType const& min_distance_sym_equiv=0.5,
    bool assert_min_distance_sym_equiv=true  ) {
    this->grid_handle_
      = generic_grid<asu,FloatType,IntType>(
        data,space_group,fasu,grid_len,
        min_distance_sym_equiv,
        assert_min_distance_sym_equiv);
  }
  void set_grid_handle ( af::versa<FloatType,af::flex_grid<> > const& data ) {
    this->grid_handle_ = generic_grid<unit_cell_tag,FloatType,IntType>(data);
  }
  void set_grid_handle (
    af::versa<FloatType, af::c_grid_padded<dimension_3> > const& data ) {
    this->grid_handle_ = generic_grid<unit_cell_tag,FloatType,IntType>(data);
  }
  void set_grid_handle (
    af::versa<FloatType,af::flex_grid<> > const& data,
    af::tiny<IntType,dimension_3> const& grid_len ) {
    this->grid_handle_
      = generic_grid<non_symmetric,FloatType,IntType>(data,grid_len);
  }
  void set_grid_handle ( grid_handle_type const& handle ) {
    this->grid_handle_ = handle;
  }
  void set_out_of_bounds_handle ( out_of_handle_type const& handle ) {
    this->out_of_handle_ = handle;
  }
  tin3_type extents () const {
    return this->extents_;
  }
  tbx_unit_cell unit_cell () const {
    return this->unit_cell_;
  }
private:
  f2g_type            frac2grid_;
  f2c_type            frac2cart_;
  g2f_type            grid2frac_;
  c2f_type            cart2frac_;
  c2g_type            cart2grid_;
  g2c_type            grid2cart_;
  tin3_type           extents_;
  mat3_type           matrix_;
  grid_handle_type    grid_handle_;
  out_of_handle_type  out_of_handle_;
  tbx_unit_cell       unit_cell_;
};

}

}

#endif//CCTBX_MAPTBX_ABSTRACT_MAP_H
