#ifndef CCTBX_MAPTBX_GENERIC_GRID_H
#define CCTBX_MAPTBX_GENERIC_GRID_H

#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>

namespace cctbx {

namespace maptbx {

template <  typename MapType,
      typename FloatType,
      typename IntType > class generic_grid;

// generic interface to the map
template < typename FloatType, typename IntType >
class generic_grid<void,FloatType,IntType> {
public:
  typedef basic_mapper<void,FloatType,IntType>      mapper_type;
  typedef generic_grid<void,FloatType,IntType>      grid_intf_type;
  typedef chiltbx::handle::handle<grid_intf_type>      handle_type;
  typedef af::ref<FloatType, af::flex_grid<> >      flex_grid_ref;
  typedef af::const_ref<FloatType, af::flex_grid<> >    const_flex_grid_ref;
  typedef af::c_grid_padded<dimension_3>          c_grid_padded_3;
  typedef af::ref<FloatType, c_grid_padded_3>        c_grid_ref;
  typedef af::const_ref<FloatType, c_grid_padded_3>    const_c_grid_ref;
  virtual handle_type as_handle () const = 0;
  virtual mapper_type remap ( fractional<FloatType> const& ) const = 0;
  virtual FloatType get_value ( grid_point<IntType> const& ) const = 0;
  virtual void set_value ( grid_point<IntType> const&, FloatType const& ) = 0;
  virtual bool is_inside ( grid_point<IntType> const& ) const = 0;
  virtual bool is_inside ( fractional<FloatType> const& ) const = 0;
  virtual af::tiny<IntType,dimension_3> origin () const = 0;
  virtual af::tiny<IntType,dimension_3> last () const = 0;
  virtual af::tiny<IntType,dimension_3> focus () const = 0;
  virtual af::tiny<IntType,dimension_3> all () const = 0;
};

namespace cdsa = crystal::direct_space_asu;
typedef sgtbx::space_group tbx_space_group;

template < typename FloatType, typename IntType >
class generic_grid<asu,FloatType,IntType>
: public generic_grid<void,FloatType,IntType> {
public:
  typedef generic_grid<void,FloatType,IntType>      grid_intf_type;
  typedef chiltbx::handle::handle<grid_intf_type>      handle_type;
  typedef basic_mapper<void,FloatType,IntType>      mapper_type;
  typedef af::ref<FloatType, af::flex_grid<> >      flex_grid_ref;
  typedef af::const_ref<FloatType, af::flex_grid<> >    const_flex_grid_ref;
  typedef af::c_grid_padded<dimension_3>          c_grid_3;
  typedef af::ref<FloatType, c_grid_3>          c_grid_ref;
  typedef af::const_ref<FloatType, c_grid_3 >        const_c_grid_ref;
  typedef typename af::flex_grid<>::index_type      map_index;
  generic_grid ( af::versa<FloatType, af::flex_grid<> > const& data,
    tbx_space_group const& space_group,
    cdsa::float_asu<FloatType> const& fasu,
    af::tiny<IntType,dimension_3> const& grid_len,
    FloatType const& min_distance_sym_equiv=0.5,
    bool assert_min_distance_sym_equiv=true )
  : space_group_(space_group)
  , float_asu_(fasu)
  , grid_length_(grid_len)
  , min_distance_sym_equiv_(min_distance_sym_equiv)
  , assert_min_distance_sym_equiv_(assert_min_distance_sym_equiv) {
    this->versa_data_ = data;
    this->versa_ref_ = this->versa_data_.ref();
    CCTBX_ASSERT(this->versa_ref_.accessor().nd() == 3);
  }
  generic_grid ( generic_grid const& gintf )
  : space_group_(gintf.space_group_)
  , float_asu_(gintf.float_asu_)
  , grid_length_(gintf.grid_length_) {
    this->versa_data_ = gintf.versa_data_;
    this->versa_ref_ = this->versa_data_.ref();
  }
  generic_grid& operator = ( generic_grid const& gintf ) {
    this->versa_data_ = gintf.versa_data_;
    this->versa_ref_ = this->versa_data_.ref();
    this->space_group_ = gintf.space_group_;
    this->float_asu_ = gintf.float_asu_;
    this->grid_length_ = gintf.grid_length_;
    return *this;
  }
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  virtual mapper_type remap ( fractional<FloatType> const& F ) const {
    return basic_mapper<asu,FloatType,IntType>(F,this->space_group_,this->float_asu_);
  }
  map_index grid_point_to_map_index ( grid_point<IntType> const& gpt ) const {
    map_index midx(3,0);
    midx[0] = gpt[0];
    midx[1] = gpt[1];
    midx[2] = gpt[2];
    return midx;
  }
  virtual FloatType get_value ( grid_point<IntType> const& coord ) const {
    return this->versa_ref_(this->grid_point_to_map_index(coord));
  }
  void basic_set_value ( grid_point<IntType> const& coord, FloatType const& F ) {
    this->versa_ref_.begin()[
      this->versa_ref_.accessor()(
        this->grid_point_to_map_index(coord))] = F;
  }
  virtual void set_value ( grid_point<IntType> const& coord, FloatType const& F ) {
    this->basic_set_value(coord,F);
  }
  virtual bool is_inside ( grid_point<IntType> const& coord ) const {
    return this->is_inside( grid_fractionalization<FloatType>(coord,this->grid_length_) );
  }
  virtual bool is_inside ( fractional<FloatType> const& coord ) const {
    return this->float_asu_.is_inside(coord);
  }
  virtual af::tiny<IntType,dimension_3> origin () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().origin();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> last () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().last();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> focus () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().focus();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> all () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().all();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
private:
  tbx_space_group                        space_group_;
  cdsa::float_asu<FloatType>             float_asu_;
  af::tiny<IntType,dimension_3>          grid_length_;
  af::versa<FloatType, af::flex_grid<> > versa_data_;
  af::ref<FloatType, af::flex_grid<> >   versa_ref_;
  FloatType                              min_distance_sym_equiv_;
  bool                                   assert_min_distance_sym_equiv_;
};

template < typename FloatType, typename IntType >
class generic_grid<unit_cell,FloatType,IntType>
: public generic_grid<void,FloatType,IntType> {
public:
  typedef generic_grid<void,FloatType,IntType>      grid_intf_type;
  typedef chiltbx::handle::handle<grid_intf_type>      handle_type;
  typedef basic_mapper<void,FloatType,IntType>      mapper_type;
  typedef af::ref<FloatType, af::flex_grid<> >      flex_grid_ref;
  typedef af::const_ref<FloatType, af::flex_grid<> >    const_flex_grid_ref;
  typedef af::c_grid_padded<dimension_3>          c_grid_3;
  typedef af::ref<FloatType, c_grid_3>          c_grid_ref;
  typedef af::const_ref<FloatType, c_grid_3 >        const_c_grid_ref;
  generic_grid ( af::versa<FloatType,af::flex_grid<> > const& data ) {
    this->c_grid_data_ = af::versa<FloatType, af::c_grid_padded<dimension_3> >(
                data,
                af::c_grid_padded<dimension_3>(
                  data.accessor()));
    this->c_grid_ref_ = this->c_grid_data_.ref();
  }
  generic_grid ( af::versa<FloatType, af::c_grid_padded<dimension_3> > const& data ) {
    this->c_grid_data_ = data;
    this->c_grid_ref_ = this->c_grid_data_.ref();
  }
  generic_grid ( generic_grid const& gintf ){
    this->versa_data_ = gintf.versa_data_;
    this->versa_ref_ = this->versa_data_.ref();
    this->c_grid_data_ = gintf.c_grid_data_;
    this->c_grid_ref_ = this->c_grid_data_.ref();
  }
  generic_grid& operator = ( generic_grid const& gintf ) {
    this->versa_data_ = gintf.versa_data_;
    this->versa_ref_ = this->versa_data_.ref();
    this->c_grid_data_ = gintf.c_grid_data_;
    this->c_grid_ref_ = this->c_grid_data_.ref();
    return *this;
  }
  virtual mapper_type remap ( fractional<FloatType> const& F ) const {
    return basic_mapper<unit_cell,FloatType,IntType>(F);
  }
  typedef typename af::c_grid_padded<dimension_3>::index_type  map_index;
  map_index grid_point_to_map_index ( grid_point<IntType> const& gpt ) const {
    map_index midx(3,0);
    midx[0] = gpt[0];
    midx[1] = gpt[1];
    midx[2] = gpt[2];
    return midx;
  }
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  virtual FloatType get_value ( grid_point<IntType> const& coord ) const {
    return this->c_grid_ref_(this->grid_point_to_map_index(coord));
  }
  void basic_set_value ( grid_point<IntType> const& coord, FloatType const& F ) {
    this->c_grid_ref_.begin()[
      this->c_grid_ref_.accessor()(
        this->grid_point_to_map_index(coord))] = F;
  }
  virtual void set_value ( grid_point<IntType> const& coord, FloatType const& F ) {
    this->basic_set_value(coord,F);
  }
  virtual bool is_inside ( grid_point<IntType> const& coord ) const {
    for ( std::size_t i=0; i<dimension_3; ++i )
      if ( coord[i] < 0 || coord[i] >= this->c_grid_ref_.accessor().focus()[i] )
        return false;
    return true;
  }
  virtual bool is_inside ( fractional<FloatType> const& coord ) const {
    // fractional space for unit_cell is the easiest of all!!!
    for ( std::size_t i=0; i<dimension_3; ++i )
      if ( coord[i] < 0 || coord[i] >= 1 )
        return false;
    return true;
  }
  virtual af::tiny<IntType,dimension_3> origin () const {
    return af::tiny<IntType,dimension_3>(0,0,0);
  }
  virtual af::tiny<IntType,dimension_3> last () const {
    af::tiny<IntType,dimension_3> result;
    typename af::c_grid_padded<dimension_3>::index_type idx
      = this->c_grid_ref_.accessor().focus();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> focus () const {
    af::tiny<IntType,dimension_3> result;
    typename af::c_grid_padded<dimension_3>::index_type idx
      = this->c_grid_ref_.accessor().focus();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> all () const {
    af::tiny<IntType,dimension_3> result;
    typename af::c_grid_padded<dimension_3>::index_type idx
      = this->c_grid_ref_.accessor().all();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
private:
  af::versa<FloatType, af::c_grid_padded<dimension_3> >  c_grid_data_;
  af::ref<FloatType, af::c_grid_padded<dimension_3> >    c_grid_ref_;
  af::versa<FloatType, af::flex_grid<> >          versa_data_;
  af::ref<FloatType, af::flex_grid<> >          versa_ref_;
};

template < typename FloatType, typename IntType >
class generic_grid<non_symmetric,FloatType,IntType>
: public generic_grid<void,FloatType,IntType> {
public:
  typedef generic_grid<void,FloatType,IntType>      grid_intf_type;
  typedef chiltbx::handle::handle<grid_intf_type>      handle_type;
  typedef basic_mapper<void,FloatType,IntType>      mapper_type;
  typedef af::ref<FloatType, af::flex_grid<> >      flex_grid_ref;
  typedef af::const_ref<FloatType, af::flex_grid<> >    const_flex_grid_ref;
  typedef af::c_grid_padded<dimension_3>          c_grid_3;
  typedef af::ref<FloatType, c_grid_3>          c_grid_ref;
  typedef af::const_ref<FloatType, c_grid_3 >        const_c_grid_ref;
  generic_grid ( af::versa<FloatType,af::flex_grid<> > const& data,
    af::tiny<IntType,dimension_3> const& grid_len )
  : grid_length_(grid_len) {
    this->versa_data_ = data;
    this->versa_ref_ = this->versa_data_.ref();
    CCTBX_ASSERT(this->versa_ref_.accessor().nd() == 3);
  }
  generic_grid ( generic_grid const& gintf )
  : grid_length_(gintf.grid_length_) {
    this->versa_data_ = gintf.versa_data_;
    this->versa_ref_ = this->versa_data_.ref();
  }
  generic_grid& operator = ( generic_grid const& gintf ) {
    this->grid_length_ = gintf.grid_length_;
    this->versa_data_ = gintf.versa_data_;
    this->versa_ref_ = this->versa_data_.ref();
    return *this;
  }
  virtual mapper_type remap ( fractional<FloatType> const& F ) const {
    return basic_mapper<non_symmetric,FloatType,IntType>(F);
  }
  typedef typename af::flex_grid<>::index_type  map_index;
  map_index grid_point_to_map_index ( grid_point<IntType> const& gpt ) const {
    map_index midx(3,0);
    midx[0] = gpt[0];
    midx[1] = gpt[1];
    midx[2] = gpt[2];
    return midx;
  }
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  virtual FloatType get_value ( grid_point<IntType> const& coord ) const {
    return this->versa_ref_(this->grid_point_to_map_index(coord));
  }
  void basic_set_value ( grid_point<IntType> const& coord, FloatType const& F ) {
    this->versa_ref_.begin()[
      this->versa_ref_.accessor()(
        this->grid_point_to_map_index(coord))] = F;
  }
  virtual void set_value ( grid_point<IntType> const& coord, FloatType const& F ) {
    this->basic_set_value(coord,F);
  }
  virtual bool is_inside ( grid_point<IntType> const& gpt ) const {
    for ( std::size_t i=0; i<dimension_3; ++i )
      if ( gpt[i] < this->versa_ref_.accessor().origin()[i]
        || gpt[i] >= this->versa_ref_.accessor().last()[i] )
        return false;
    return true;
  }
  virtual bool is_inside ( fractional<FloatType> const& coord ) const {
    grid_point<IntType> g_origin(this->versa_ref_.accessor().origin().begin());
    grid_point<IntType> g_last(this->versa_ref_.accessor().last().begin());
    fractional<FloatType> origin = grid_fractionalization<FloatType>(g_origin,this->grid_length_);
    fractional<FloatType> last = grid_fractionalization<FloatType>(g_last,this->grid_length_);
    for ( std::size_t i=0; i<dimension_3; ++i )
      if ( coord[i] < origin[i] || coord[i] >= last[i] )
        return false;
    return true;
  }
  virtual af::tiny<IntType,dimension_3> origin () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().origin();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> last () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().last();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> focus () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().focus();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
  virtual af::tiny<IntType,dimension_3> all () const {
    af::tiny<IntType,dimension_3> result;
    typename af::flex_grid<>::index_type idx = this->versa_ref_.accessor().all();
    for ( std::size_t i=0; i<dimension_3; ++i )
      result[i] = idx[i];
    return result;
  }
private:
  af::tiny<IntType,dimension_3>      grid_length_;
  af::versa<FloatType, af::flex_grid<> >  versa_data_;
  af::ref<FloatType, af::flex_grid<> >  versa_ref_;
};

}

}

#endif//CCTBX_MAPTBX_GENERIC_GRID_H
