// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     Created: 12-Jul-2001 (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_ASYM_UNITS_H
#define CCTBX_SGTBX_ASYM_UNITS_H

#include <cctbx/sgtbx/groups.h>

namespace cctbx { namespace sgtbx {

  //! A "point" of a Brick.
  class BrickPoint {
    public:
      //! Default constructor. Some data members are not initialized.
      BrickPoint() {};
      //! The point as a rational number.
      rational Point() const { return m_Point; }
      //! Flag to indicate if the point is inside the brick or just outside.
      /*! Examples for a min point:
          <ul>
          <li>Point() == 1/8 and Off() == false: 1/8 <= x.
          <li>Point() == 1/8 and Off() == true: 1/8 < x.
          </ul>
          For a max point the inequalities would be 1/8 >= x and 1/8 > x,
          respectively.
       */
      int Off() const { return m_Off; }
    private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      friend class Brick;
#endif // DOXYGEN_SHOULD_SKIP_THIS
      explicit BrickPoint(int RawPoint);
      rational m_Point;
      bool m_Off;
  };

  //! Parallelepiped that contains an asymmetric unit.
  /*! A "Brick" is a parallelepiped that contains an asymmetric unit
      and is optimized for memory allocation.
      <p>
      A Brick is the result of SpaceGroup::getBrick().
      <p>
      Bricks for 530 conventional settings (see class
      SpaceGroupSymbols) and an additional 223 primitive settings of
      centred space groups were computed with sginfo2 and were
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
      <p>
      Remark: The 126 tabluated Hermann-Mauguin symbols for the
      "Multiple cell C or F" and "Triple cell H" settings were added in
      the sgtbx and are not available in the older sginfo2 program.
      Therefore the bricks for these settings are currently missing.
   */
  class Brick {
    public:
      //! Default constructor. Some data members are not initialized.
      Brick() {}
      //! Constructor.
      explicit
      Brick(const SpaceGroupInfo& SgInfo);
      //! Access to the six points of the brick.
      /*! iBasisVector refers to the basis vectors a,b,c.
          <p>
          iMinMax == 0 selects the min point along the selected basis
                       vector.<br>
          iMinMax == 1 selects the max point.<br>
          <p>
          See BrickPoint for more details.
       */
      BrickPoint
      operator()(std::size_t iBasisVector, std::size_t iMinMax) const;
      //! Format the information about the Brick as a string.
      /*! Example: 0<=x<=1/8; -1/8<=y<=0; 1/8<z<7/8
       */
      std::string format() const;
      //! Test if a given point is inside the brick.
      bool isInBrick(const carray<rational, 3>& P) const;
    private:
      BrickPoint m_Points[3][2];
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_ASYM_UNITS_H
