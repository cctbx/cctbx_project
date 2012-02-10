#ifndef SCITBX_MATH_POLYGON_INTERSECTION_H
#define SCITBX_MATH_POLYGON_INTERSECTION_H

#include <scitbx/array_family/shared.h>
#include <scitbx/vec2.h>

namespace scitbx { namespace math {

  template <typename DataType = double>
  struct point_to_line_distance {
    scitbx::vec2<DataType> avec,bvec;
    DataType ab_length;
    point_to_line_distance(const scitbx::vec2<DataType>& A, const scitbx::vec2<DataType>& B):
      avec(A),bvec(B){
      ab_length = (avec-bvec).length();
    }
    DataType signed_distance(const scitbx::vec2<DataType>& P)const{
      return
        ( (P[0]-avec[0])*(bvec[1]-avec[1])-(P[1]-avec[1])*(bvec[0]-avec[0]) )/ab_length;
    }

  };

  //! Calculate whether two convex polygons intersect
  /*! Two dimensional case.
      It is the responsibility of the caller to assure that the
        polygons are convex.
      Also, the vectors describing the vertices must be given in
        sequential order, either or clockwise or counterclockwise.
      Algorithm inspired by Greg Magarshak (2000) Theory & Practice, issue 02,
        Collision Detection, part 2 (
        http://www.flipcode.com/archives/Theory_Practice-Issue_02_Collision_Detection_Part_2.shtml
       )
   */
  bool convex_polygons_intersect_2D (
    af::shared<scitbx::vec2<double> > polygon1,
    af::shared<scitbx::vec2<double> > polygon2
    )
  {
    int size1=polygon1.size();
    int size2=polygon2.size();

    //step 1.  try to make separation lines from the polygon1 edges; i.e., is it possible to find
    //        a line segment such that polygon1 is completely on one side, and polygon2 is on the other

    for (int pedge = 0; pedge< size1; ++pedge){
      int next_pt = (pedge == size1-1)? 0 : pedge+1;
      point_to_line_distance<double> candidate_separator(polygon1[pedge],polygon1[next_pt]);
      int ref_idx = (next_pt== size1-1)? 0 : next_pt + 1;
      double ref_distance = candidate_separator.signed_distance(polygon1[ref_idx]);
      bool found_separator_line = true;
      for (int itest = 0; itest < size2; ++itest){
        if (ref_distance > 0. ==
             (candidate_separator.signed_distance(polygon2[itest])>0.)) {
            found_separator_line = false; break;
        }
      }
      if (found_separator_line) return false;
    }

    //step 2.  same thing again with polygon2 edges
    for (int pedge = 0; pedge< size2; ++pedge){
      int next_pt = (pedge == size2-1)? 0 : pedge+1;
      point_to_line_distance<double> candidate_separator(polygon2[pedge],polygon2[next_pt]);
      int ref_idx = (next_pt== size2-1)? 0 : next_pt + 1;
      double ref_distance = candidate_separator.signed_distance(polygon2[ref_idx]);
      bool found_separator_line = true;
      for (int itest = 0; itest < size1; ++itest){
        if (ref_distance > 0. ==
             (candidate_separator.signed_distance(polygon1[itest])>0.)) {
            found_separator_line = false; break;
        }
      }
      if (found_separator_line) return false;
    }

    return true;
  }

}} //namespace

#endif // GUARD
