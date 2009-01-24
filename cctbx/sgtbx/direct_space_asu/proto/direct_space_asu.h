#ifndef CCTBX_SGTBX_DIRECT_SPACE_ASU_H 
#define CCTBX_SGTBX_DIRECT_SPACE_ASU_H

#include "abstract.h"

namespace cctbx { namespace sgtbx { namespace asu {

  //! Direct space asymmetric unit
  class direct_space_asu
  {
  public:

    //! one template expression with all the faces
    abstract::ptr faces;
    std::string hall_symbol;

    void get_nth_plane(size_type i, cut &face) const
    {
      faces->get_nth_plane(i, face);
    }

    void show_summary(std::ostream &os) const;

    void show_comprehensive_summary(std::ostream &os) const;

    //! Removes subexpressions from every face
    void volume_only()
    {
      faces = faces->new_volume_only();
    }

    void in_which_planes(const rvector3_t &p, std::vector<cut> &result) const;

    //! Tests if rational 3D point belongs to the asymmetric unit
    bool is_inside(const rvector3_t &p) const
    {
      return faces->is_inside(p);
    }

    bool is_inside_volume_only(const rvector3_t &p) const
    {
      const size_type sz = this->n_faces();
      cut plane;
      for(size_type i=0; i<sz; ++i )
      {
        this->get_nth_plane(i, plane);
        if( plane.evaluate(p)<0 )
          return false;
      }
      return true;
    }

    //! Returns number of the faces in the asu
    size_type n_faces() const {return faces->size(); }

    void change_basis(const change_of_basis_op &op)
    {
      std::string new_hall;
      if( !hall_symbol.empty() )
      {
        new_hall = cctbx::sgtbx::space_group("Hall: " + hall_symbol).change_basis(op).type().hall_symbol();
      }
      hall_symbol = new_hall;
      faces->change_basis(op);
    }

    //! Returns a set of all asu vertices
    void volume_vertices(set_rvector3_t &result) const;

    //! Returns bounding box for the asu
    bool box_corners(rvector3_t &mn, rvector3_t &mx) const;

    direct_space_asu(const direct_space_asu &a) : faces(a.faces->new_copy()), hall_symbol(a.hall_symbol) {}
    direct_space_asu() : faces(NULL), hall_symbol() {}

    explicit direct_space_asu(unsigned short spgr)
    {
      CCTBX_ASSERT( spgr>0 && spgr < 231 );
      hall_symbol = space_group_symbols(spgr).hall();
      faces = asu_table[spgr-1]();
    }

    explicit direct_space_asu(const std::string &spgr) : faces(NULL), hall_symbol()
    {
      // is this right?
      space_group_symbols s(spgr);
      new(this) direct_space_asu(s.number());
      std::string cbs = s.change_of_basis_symbol();
      if( !cbs.empty() )
        this->change_basis(change_of_basis_op(cbs));
    }

    direct_space_asu operator= (const direct_space_asu &a) 
    {
      faces= abstract::ptr(a.faces->new_copy()); return *this; 
      hall_symbol = a.hall_symbol;
    }

  }; // class direct_space_asu


}}}
#endif

