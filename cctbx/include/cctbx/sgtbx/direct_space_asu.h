#ifndef CCTBX_SGTBX_DIRECT_SPACE_ASU_H
#define CCTBX_SGTBX_DIRECT_SPACE_ASU_H

#include <cctbx/uctbx.h>
#include <cctbx/error.h>

namespace cctbx { namespace sgtbx { namespace direct_space_asu {

  //! Floating-point parameterization of a cut plane.
  template <typename FloatType=double>
  class float_cut_plane
  {
    public:
      //! Default constructor. Some data members are not initialized!
      float_cut_plane() {}

      //! Initialization with normal vector n and constant c.
      float_cut_plane(
        fractional<FloatType> const& n_,
        FloatType const& c_)
      :
        n(n_),
        c(c_)
      {}

      //! Normal vector.
      fractional<FloatType> n;

      //! Constant.
      FloatType c;

      //! Returns n * point + c.
      FloatType
      evaluate(fractional<FloatType> const& point) const
      {
        return n * point + c;
      }

      //! Equivalent to evaluate(point) >= 0.
      bool
      is_inside(fractional<FloatType> const& point) const
      {
        if (evaluate(point) < 0) return false;
        return true;
      }

      //! Returns a point in the cut plane.
      fractional<FloatType>
      get_point_in_plane() const
      {
        fractional<FloatType> result(0,0,0);
        for(std::size_t i=0;i<3;i++) {
          if (n[i] != 0) {
            result[i] = -c / n[i];
            return result;
          }
        }
        throw error("float_cut_plane normal vector is the null vector.");
      }

      //! Shifts the plane by the distance specified as thickness.
      float_cut_plane
      add_buffer(
        uctbx::unit_cell const& unit_cell,
        FloatType const& thickness) const
      {
        typedef fractional<FloatType> f_t;
        typedef cartesian<FloatType> c_t;
        f_t x_frac = get_point_in_plane();
        c_t x_cart = unit_cell.orthogonalize(x_frac);
        c_t n_cart = (n * unit_cell.fractionalization_matrix()).normalize();
        c_t y_cart = x_cart - n_cart * thickness;
        f_t y_frac = unit_cell.fractionalize(y_cart);
        return float_cut_plane(n, -n * y_frac);
      }
  };

  //! Floating-point parameterization of an asymmetric unit.
  template <typename FloatType=double>
  class float_asu
  {
    public:
      //! Array type for facets.
      typedef af::small<float_cut_plane<FloatType>, 9> facets_t;

      //! Default constructor. Some data members are not initialized!
      float_asu() {}

      //! Initialization with unit cell and list of facets.
      float_asu(uctbx::unit_cell const& unit_cell, facets_t const& facets)
      :
        unit_cell_(unit_cell),
        facets_(facets)
      {}

      //! Unit cell as passed to the constructor.
      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      //! Facets as passed to the constructor.
      facets_t const&
      facets() const { return facets_; }

      //! True if is-inside test is true for all facets().
      bool
      is_inside(fractional<FloatType> const& point) const
      {
        for(std::size_t i=0;i<facets_.size();i++) {
          if (!facets_[i].is_inside(point)) return false;
        }
        return true;
      }

      /*! New asymmetric unit with all facets shifted by the distance
          specified as thickness.
       */
      float_asu
      add_buffer(FloatType const& thickness) const
      {
        facets_t buffer_facets;
        for(std::size_t i=0;i<facets_.size();i++) {
          buffer_facets.push_back(
            facets_[i].add_buffer(unit_cell_, thickness));
        }
        return float_asu(unit_cell_, buffer_facets);
      }

    protected:
      uctbx::unit_cell unit_cell_;
      facets_t facets_;
  };

}}} // namespace cctbx::sgtbx::direct_space_asu

#endif // CCTBX_SGTBX_DIRECT_SPACE_ASU_H
