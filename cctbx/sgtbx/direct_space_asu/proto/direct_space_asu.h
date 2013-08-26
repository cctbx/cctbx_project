#ifndef CCTBX_SGTBX_DIRECT_SPACE_ASU_H
#define CCTBX_SGTBX_DIRECT_SPACE_ASU_H

#include <scitbx/array_family/shared.h>
#include <cctbx/crystal/direct_space_asu.h>

#include "small_vec_math.h"

#include "facet_collection.h"

namespace cctbx { namespace sgtbx { namespace asu {

  enum intersection_kind { none, fully, partially };

  //! Direct space asymmetric unit
  /*! \class direct_space_asu direct_space_asu.h cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h
   */
  class direct_space_asu
  {
  public:

    const facet_collection::pointer &get_faces() const {return faces;}

    //! Hall symbol
    std::string hall_symbol;

    //! Returns i'th face of the asymmetric unit
    void get_nth_plane(size_type i, cut &face) const
    {
      faces->get_nth_plane(i, face);
    }

    //! Returns 1D tolerance based on 3D
    double get_tolerance(const scitbx::af::double3 &tol3d) const
    {
      return faces->get_tolerance(tol3d);
    }

    void show_summary(std::ostream &os) const;

    std::string as_string() const;

    void show_comprehensive_summary(std::ostream &os) const
    {
      this->show_summary(os);
      os << '\n';
      faces->print(os);
    }

    void print_faces_as_xyz(std::ostream &os) const
    {
      this->faces->print_as_xyz(os);
    }

    //! Removes subexpressions from every face
    void shape_only()
    {
      faces = faces->new_shape_only();
    }

    // Experimental
    void shape_only_keep_inclusive_flag()
    {
      faces = faces->new_shape_only_keep_inclusive_flag();
    }


    //! Return set of planes which contain point
    void in_which_planes(const rvector3_t &point, std::vector<cut> &planes) const;

    scitbx::af::shared<cut> in_which_facets(const rvector3_t &point) const
    {
      // scitbx::af::shared<cut> planes;
      std::vector<cut> planes;
      this->in_which_planes(point,  planes);
      return scitbx::af::shared<cut> ( &*planes.begin(), &*planes.end() );
    }

    //! Tests if rational 3D point belongs to the asymmetric unit
    bool is_inside(const rvector3_t &p) const
    {
      return faces->is_inside(p);
    }

    //! Tests if num/den belongs to the asymmetric unit
    /*! This function should be faster than it's rational variant.
     * It imposes a constrain on the den. Approximately:
     * den[0]*den[1]*den[2]*4*8 < max_int.
     */
    bool is_inside(const scitbx::int3 &num, const scitbx::int3 &den) const
    {
      return faces->is_inside(num,den);
    }

    // this should go away?
    short where_is(const scitbx::int3 &num, const scitbx::int3 &den) const
    {
      return faces->where_is(num, den);
    }


    //! Tests if point belongs to the asymmetric unit, disregarding plane subexpressions
    bool is_inside_shape_only(const rvector3_t &point) const
    {
      const size_type sz = this->n_faces();
      cut plane;
      for(size_type i=0; i<sz; ++i )
      {
        this->get_nth_plane(i, plane);
        if( plane.evaluate(point)<0 )
          return false;
      }
      return true;
    }

    bool is_inside_shape_only(const scitbx::af::double3 &point, double tol) const
    {
      return faces->is_inside_shape_only(point, tol);
    }

    bool is_inside(const rvector3_t &point, bool vol_only) const
    {
      return vol_only ? is_inside_shape_only(point) : is_inside(point);
    }

    //! Returns number of the faces in the asu
    size_type n_faces() const {return faces->size(); }

    //! Changes space group basis of the asu
    void change_basis(const change_of_basis_op &op)
    {
      std::string new_hall;
      if( !hall_symbol.empty() )
        new_hall = space_group(hall_symbol).change_basis(op).type().hall_symbol();
      hall_symbol = new_hall;
      faces->change_basis(op);
    }

    //! Returns a set of all asu vertices
    void shape_vertices(set_rvector3_t &result) const;

    //! Returns bounding box for the asu. Not in python
    void box_corners(rvector3_t &mn, rvector3_t &mx) const;

    //! Computes boundaries of the box enclosed by the asu.
    /*! Boundaries are computed for a grid size 'grid'.
     * The enclosed box will be completely inside the asu,
     * and will not touch the surfaces of the asu.
     * Returns true if enclosed box found.
     * Currently will succeed only for parallelepiped like asu.
     * This function is not available in python.
     */
    bool enclosed_box_corners(scitbx::int3 &mn, scitbx::int3 &mx, const scitbx::int3 &grid) const;

    //! As box_corners for python
    rvector3_t box_max() const
    {
      rvector3_t mn, mx;
      this->box_corners(mn,mx);
      return mx;
    }

    //! As box_corners for python
    rvector3_t box_min() const
    {
      rvector3_t mn, mx;
      this->box_corners(mn,mx);
      return mn;
    }

    intersection_kind does_intersect(const scitbx::double3 &center, const scitbx::double3 &box) const;
    void get_adjacent_cells(std::vector<scitbx::tiny3> &cells) const;
    void get_cells(std::vector<scitbx::tiny3> &cells) const;
    rvector3_t move_inside(const cctbx::sgtbx::space_group &group, const rvector3_t &v) const;

    //! Converts to float_asu
    cctbx::crystal::direct_space_asu::float_asu<> as_float_asu(
      const cctbx::uctbx::unit_cell &cell,
      double epsilon=1.0E-6) const;

    direct_space_asu(const direct_space_asu &a) : hall_symbol(a.hall_symbol),
      faces(a.faces->new_copy()) {}
    direct_space_asu() : hall_symbol(), faces(NULL) {}

    //! Creates asymmetric unit from space group type
    explicit direct_space_asu(const space_group_type &group_type)
      : hall_symbol(group_type.hall_symbol()),
        faces(asu_table[group_type.number()-1]()) // build reference spacegroup asu
    {
      change_of_basis_op  op(  group_type.cb_op().inverse() );
      CCTBX_ASSERT( faces.get() != NULL );
      if( !op.is_identity_op() )
        faces->change_basis(op); // change to the real space group
    }

    // DO NOT USE!!!  For debugging purposes only!
    explicit direct_space_asu(const space_group_type &group_type,
      facet_collection::pointer &faces_)
      : hall_symbol(group_type.hall_symbol()),
        faces(faces_) // build custom asu
    {
      change_of_basis_op  op(  group_type.cb_op().inverse() );
      CCTBX_ASSERT( faces.get() != NULL );
      if( !op.is_identity_op() )
        faces->change_basis(op); // change to the real space group
    }

    // DO NOT USE!!!
    void add_face(const cut &face)
    {
      faces = faces->add_face(face);
    }

    //! Creates asymmetric unit from space group symbol
    explicit direct_space_asu(const std::string &group_symbol) :  hall_symbol(),
      faces(NULL)
    {
      // new(this) direct_space_asu( space_group_type(spgr) );  this fails
      *this =  direct_space_asu( space_group_type(group_symbol) );
    }

    direct_space_asu& operator= (const direct_space_asu &a)
    {
      faces = facet_collection::pointer(a.faces->new_copy());
      hall_symbol = a.hall_symbol;
      return *this;
    }

  private:

    // one template expression with all the faces
    facet_collection::pointer faces;

  }; // class direct_space_asu


}}}
#endif
