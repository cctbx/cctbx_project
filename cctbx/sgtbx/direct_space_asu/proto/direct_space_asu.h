#ifndef CCTBX_SGTBX_DIRECT_SPACE_ASU_H
#define CCTBX_SGTBX_DIRECT_SPACE_ASU_H

#include "facet_collection.h"

namespace cctbx { namespace sgtbx { namespace asu {

  //! Direct space asymmetric unit
  /*! \class direct_space_asu direct_space_asu.h cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h
   */
  class direct_space_asu
  {
  public:

    //! Hall symbol
    std::string hall_symbol;

    //! Returns i'th face of the asymmetric unit
    void get_nth_plane(size_type i, cut &face) const
    {
      faces->get_nth_plane(i, face);
    }

    void show_summary(std::ostream &os) const;

    std::string as_string() const;

    void show_comprehensive_summary(std::ostream &os) const
    {
      this->show_summary(os);
      os << '\n';
      faces->print(os);
    }

    //! Removes subexpressions from every face
    void volume_only()
    {
      faces = faces->new_volume_only();
    }

    //! Return set of planes which contain point
    void in_which_planes(const rvector3_t &point, std::vector<cut> &planes) const;

    //! Tests if rational 3D point belongs to the asymmetric unit
    bool is_inside(const rvector3_t &p) const
    {
      return faces->is_inside(p);
    }

    //! Tests if point belongs to the asymmetric unit, disregarding plane subexpressions
    bool is_inside_volume_only(const rvector3_t &point) const
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
    void volume_vertices(set_rvector3_t &result) const;

    //! Returns bounding box for the asu. Not in python
    void box_corners(rvector3_t &mn, rvector3_t &mx) const;

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

    direct_space_asu(const direct_space_asu &a) : hall_symbol(a.hall_symbol), faces(a.faces->new_copy()) {}
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

    //! Creates asymmetric unit from space group symbol
    explicit direct_space_asu(const std::string &group_symbol) :  hall_symbol(), faces(NULL)
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

