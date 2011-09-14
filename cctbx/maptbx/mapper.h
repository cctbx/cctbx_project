// copyright (c) Jacob N. Smith; leave this here; use at your whim
#ifndef CCTBX_MAPTBX_MAPPER_H
#define CCTBX_MAPTBX_MAPPER_H

#include<chiltbx/handle.h>
#include<scitbx/math/utils.h>
#include<scitbx/array_family/misc_functions.h>
#include<cctbx/crystal/direct_space_asu.h>
#include<cctbx/sgtbx/space_group.h>
#include<cctbx/maptbx/coordinate_transformers.h>

namespace cctbx {

namespace maptbx {

typedef sgtbx::space_group tbx_space_group;
namespace cdsa = crystal::direct_space_asu;

struct non_symmetric {};
struct unit_cell {};
struct asu {};

// converts a fractional into the space
template <  typename Tag,
      typename FloatType,
      typename IntType > struct basic_mapper;

template < typename FloatType, typename IntType >
struct basic_mapper<void,FloatType,IntType> {
  basic_mapper ( fractional<FloatType> const& coord,
    grid_point<IntType> const& ushft,
    IntType const& symop )
  : mapped_coordinate(coord)
  , unit_shift(ushft)
  , symmetry_operation(symop) {
  }
  basic_mapper () {}
  fractional<FloatType>  mapped_coordinate;
  grid_point<IntType>    unit_shift;
  IntType          symmetry_operation;
};

// non_symmetric one-shot mapper; fractional -> fractional
template < typename FloatType, typename IntType >
struct basic_mapper<non_symmetric,FloatType,IntType>
: basic_mapper<void,FloatType,IntType> {
  basic_mapper ( fractional<FloatType> const& coord ) {
    this->mapped_coordinate    = coord;
    this->unit_shift      = grid_point<IntType>(0,0,0);
    this->symmetry_operation  = IntType(-1);
  }
};

// unit_cell one-shot mapper; fractional -> fractional
template < typename FloatType, typename IntType >
struct basic_mapper<unit_cell,FloatType,IntType>
: basic_mapper<void,FloatType,IntType> {
  typedef scitbx::math::float_int_conversions<FloatType,IntType> fic;
  basic_mapper ( fractional<FloatType> const& coord ) {
    for ( std::size_t i=0; i<dimension_3; ++i ) {
      this->mapped_coordinate[i] = scitbx::fn::fmod_positive(coord[i],FloatType(1));
      this->unit_shift[i] = IntType( fic::ifloor(coord[i]-this->mapped_coordinate[i]) );
    }
    this->symmetry_operation = -1;
  }
};

// asu one-shot mapper; fractional -> fractional
template < typename FloatType, typename IntType >
struct basic_mapper<asu,FloatType,IntType>
: basic_mapper<void,FloatType,IntType> {
  basic_mapper ( fractional<FloatType> const& coord,
    tbx_space_group const& space_group,
    cdsa::float_asu<FloatType> const& asu,
    FloatType const& min_distance_sym_equiv=0.5,
    bool assert_min_distance_sym_equiv=true ) {

    sgtbx::site_symmetry sitesym(
      asu.unit_cell(),
      space_group,
      coord,
      min_distance_sym_equiv,
      assert_min_distance_sym_equiv);

    sgtbx::sym_equiv_sites<FloatType> equiv_sites(sitesym);

    af::const_ref<typename sgtbx::sym_equiv_sites<FloatType>::coor_t>
      coordinates = equiv_sites.coordinates().const_ref();
    af::const_ref<std::size_t>
      sym_op_indices = equiv_sites.sym_op_indices().const_ref();

    for ( std::size_t i_sym_eq=0; i_sym_eq<coordinates.size(); ++i_sym_eq ) {

      scitbx::vec3<FloatType> const& site = coordinates[i_sym_eq];

      // find the minimum and maximum unit shifts for this symmetry operation
      scitbx::vec3<IntType> unit_shifts_min;
      scitbx::vec3<IntType> unit_shifts_max;
      for(std::size_t i=0;i<3;i++) {
        unit_shifts_min[i] = scitbx::math::iceil(
          asu.box_min()[i] - site[i]
          - 2*asu.is_inside_epsilon());
        unit_shifts_max[i] = scitbx::math::ifloor(
          asu.box_max()[i] - site[i]
          + 2*asu.is_inside_epsilon());
      }

      scitbx::vec3<IntType> u;
      scitbx::vec3<FloatType> mapped_site;

      // apply a combination of unit shifts to find a value inside
      // the asu
      for(u[0]=unit_shifts_min[0];u[0]<=unit_shifts_max[0];u[0]++) {
        mapped_site[0] = site[0] + u[0];
        for(u[1]=unit_shifts_min[1];u[1]<=unit_shifts_max[1];u[1]++) {
          mapped_site[1] = site[1] + u[1];
          for(u[2]=unit_shifts_min[2];u[2]<=unit_shifts_max[2];u[2]++) {
            mapped_site[2] = site[2] + u[2];
            if ( asu.is_inside(mapped_site) ) {
              this->symmetry_operation = sym_op_indices[i_sym_eq];
              this->unit_shift = u;
              this->mapped_coordinate = mapped_site;
              return;
            }
          }
        }
      }
    }

    // failure
    throw error("basic_mapper<asu>: asu symmetry mapping failed internally.");
  }
};

template <  typename SymmetryType,
      typename FloatType,
      typename IntType > struct mapper_factory;

template < typename FloatType, typename IntType >
struct mapper_factory<void,FloatType,IntType> {
  typedef basic_mapper<void,FloatType,IntType>  mapper_type;
  typedef mapper_factory<void,FloatType,IntType>  factory_type;
  typedef chiltbx::handle::handle<factory_type>  factory_handle;
  virtual factory_handle as_handle () const = 0;
  virtual mapper_type map ( fractional<FloatType> const& ) const = 0;
};

template < typename FactoryType > struct mapper_factory_types;
template < typename SymmetryType, typename FloatType, typename IntType >
struct mapper_factory_types< mapper_factory<SymmetryType,FloatType,IntType> > {
  typedef SymmetryType  symmetry_type;
  typedef FloatType    float_type;
  typedef IntType      int_type;
};

// unit_cell and non_symmetric factory; fractional -> fractional
template <  typename SymmetryType,
      typename FloatType,
      typename IntType >
struct mapper_factory : public mapper_factory<void,FloatType,IntType> {
  typedef basic_mapper<SymmetryType,FloatType,IntType>  basic_mapper_type;
  typedef basic_mapper<void,FloatType,IntType>      mapper_type;

  mapper_factory () {}

  virtual mapper_type map ( fractional<FloatType> const& coordinate ) const {
    basic_mapper_type bmt(coordinate);
    return mapper_type( bmt.mapped_coordinate,
              bmt.unit_shift,
              bmt.symmetry_operation);
  }

  typedef mapper_factory<void,FloatType,IntType>  factory_type;
  typedef chiltbx::handle::handle<factory_type>  factory_handle;
  virtual factory_handle as_handle () const {
    return factory_handle(*this);
  }
};

// asu factory; fractional -> fractional
template < typename FloatType, typename IntType >
struct mapper_factory<asu,FloatType,IntType>
: public mapper_factory<void,FloatType,IntType> {
  typedef basic_mapper<asu,FloatType,IntType>        basic_mapper_type;
  typedef basic_mapper<void,FloatType,IntType>      mapper_type;

  mapper_factory () {}

  mapper_factory ( tbx_space_group const& space_group,
    cdsa::float_asu<FloatType> const& asu,
    FloatType const& min_distance_sym_equiv=0.5,
    bool assert_min_distance_sym_equiv=true )
  : space_group_(space_group)
  , asu_(asu)
  , min_distance_sym_equiv_(min_distance_sym_equiv)
  , assert_min_distance_sym_equiv_(assert_min_distance_sym_equiv) {}

  virtual mapper_type map ( fractional<FloatType> const& coordinate ) const {
    basic_mapper_type bmt(coordinate,
                          this->space_group_,
                          this->asu_,
                          this->min_distance_sym_equiv_,
                          this->assert_min_distance_sym_equiv_);
    return mapper_type(  bmt.mapped_coordinate,
              bmt.unit_shift,
              bmt.symmetry_operation);
  }

  typedef mapper_factory<void,FloatType,IntType>  factory_type;
  typedef chiltbx::handle::handle<factory_type>  factory_handle;
  virtual factory_handle as_handle () const {
    return factory_handle(*this);
  }

  tbx_space_group space_group_;
  cdsa::float_asu<FloatType> asu_;
  FloatType min_distance_sym_equiv_;
  bool assert_min_distance_sym_equiv_;
};

// transformer-mapper pipeline
template < typename ToFracType, typename FactoryType, typename FromFracType >
typename transformer_types<FromFracType>::to_type
transformer_mapper ( ToFracType const& to_frac,
  FactoryType const& factory,
  FromFracType const& from_frac,
  typename transformer_types<ToFracType>::from_type const& coordinate ) {
  return from_frac(factory.map(to_frac(coordinate)).mapped_coordinate);
}

// this does not have an interface to python
// it is emulated in python
// allows mapping from any arbitrary to to any other arbitrary type
template < typename ToFracType, typename FactoryType, typename FromFracType >
struct transformer_mapper_factory {
  typedef typename transformer_types<ToFracType>::from_type  from_type;
  typedef typename transformer_types<ToFracType>::to_type    frac_type;
  typedef typename transformer_types<FromFracType>::to_type  to_type;

  transformer_mapper_factory () {}

  transformer_mapper_factory (
    ToFracType const& tofrac,
    FactoryType const& fac,
    FromFracType const& fromfrac )
  : to_frac_(tofrac)
  , factory_(fac)
  , from_frac_(fromfrac) {
  }

  to_type map ( from_type const& coordinate ) const {
    return this->from_frac_(
      this->factory_.map(
        this->to_frac_(coordinate)).mapped_coordinate);
  }

  ToFracType    to_frac_;
  FactoryType    factory_;
  FromFracType  from_frac_;
};

}// end maptbx

}// end cctbx

#endif//CCTBX_MAPTBX_MAPPER_H
