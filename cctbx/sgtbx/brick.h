#ifndef CCTBX_SGTBX_BRICK_H
#define CCTBX_SGTBX_BRICK_H

#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace sgtbx {

  //! A "point" of a brick.
  class brick_point
  {
    public:
      //! Default constructor. Some data members are not initialized.
      brick_point() {};

      //! The point as a rational number.
      rat
      value() const { return value_; }

      //! Flag to indicate if the point is inside the brick or just outside.
      /*! Examples for a min point:
          <ul>
          <li>value() == 1/8 and off() == false: 1/8 <= x.
          <li>value() == 1/8 and off() == true: 1/8 < x.
          </ul>
          For a max point the inequalities would be 1/8 >= x and 1/8 > x,
          respectively.
       */
      int
      off() const { return off_; }

    protected:
      friend class brick;
      explicit brick_point(int raw_point);
      rat value_;
      bool off_;
  };

  //! Parallelepiped that contains an asymmetric unit.
  /*! A "brick" is a parallelepiped that contains an asymmetric unit
      and is optimized for memory allocation.
      <p>
      Bricks for 530 conventional settings (see class
      space_group_symbols) and an additional 223 primitive settings of
      centred space groups were computed with sginfo2 and
      formatted as a table for C++ with a Python script. Currently, the
      algorithm for computing the bricks is not available in the
      %cctbx. However, given the large number of tabulated bricks,
      this should hardly ever be noticeable.
      <p>
      The bricks minimize the memory that has to be allocated for
      storing a part of a map (e.g. an electron density map) that
      covers an asymmetric unit. However, in general, a brick will
      contain more than exactly one asymmetric unit.
      <p>
      A method for refining a brick is to allocate a map that contains
      flags indicating whether or not a certain point inside the brick
      is in the asymmetric unit or is redundant.  It is straightforward
      to generate such a map of flags by looping over the symmetry
      operations for each grid point.  One point is marked as being in
      the asymmetric unit, and all symmetry mates are marked as being
      outside.
      <p>
      Unfortunately, this simple algorithm does not necessarily lead to
      a contiguous asymmetric unit. If this could be achieved, it would
      not be necessary to allocate an entire 3-dimensional map of
      flags, but a 2-dimensional grid could be used. Each point in the
      2-dimensional grid would contain the indices of the first and the
      last grid point in the third dimension that is inside the
      asymmetric unit. Attempts to devise an algorithm that produces
      bricks with contiguous asymmetric units have not been successful.
      Ideas are very welcome. Until then, the bricks provided by this
      class are the best solution available.
   */
  class brick
  {
    public:
      //! Default constructor. Some data members are not initialized.
      brick() {}

      //! Constructor.
      explicit
      brick(sgtbx::space_group_type const& space_group_type);

      //! Access to the six points of the brick.
      /*! i_axis refers to the basis vectors a,b,c of the unit cell.
          <p>
          i_min_max == 0 selects the min point along the selected basis
                         vector.<br>
          i_min_max == 1 selects the max point.<br>
          <p>
          See also: class brick_point
          <p>
          Not available in Python.
       */
      brick_point
      operator()(std::size_t i_axis, std::size_t i_min_max) const;

      //! Formats the information about the brick as a string.
      /*! Example: 0<=x<=1/8; -1/8<=y<=0; 1/8<z<7/8
       */
      std::string
      as_string() const;

      //! Tests if a given point is inside the brick.
      /*! Not available in Python.
       */
      bool
      is_inside(vec3_rat const& point) const;

      //! Tests if a given point is inside the brick.
      bool
      is_inside(tr_vec const& point) const;

    private:
      brick_point points_[3][2];
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_BRICK_H
