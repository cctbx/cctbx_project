#include <ostream>
#include <sstream>

#include "direct_space_asu.h"

namespace cctbx { namespace sgtbx { namespace asu {

  void direct_space_asu::show_summary(std::ostream &os) const
  {
    const size_type sz = this->n_faces();
    os << "Hall symbol: " << hall_symbol
       << "\nNumber of facets: " << sz;
  }

  std::string direct_space_asu::as_string() const
  {
    std::stringstream str;
    this->show_comprehensive_summary(str);
    return str.str();
  }

  void direct_space_asu::in_which_planes(const rvector3_t &p, std::vector<cut> &result) const
  {
    result.clear();
    const size_type sz = this->n_faces();
    cut plane;
    for(size_type i=0; i<sz; ++i)
    {
      this->get_nth_plane(i,plane);
      if( plane.evaluate(p)==0 )
        result.push_back(plane);
    }
  }

  void direct_space_asu::volume_vertices(set_rvector3_t &result) const
  {
    result.clear();
    size_type n_facets = this->n_faces();
    for(size_type i0=0; i0<n_facets-2; ++i0)  // in xrange(0,n_facets-2):
    {
      cut face0;
      faces->get_nth_plane(i0, face0);
      for(size_type i1=i0+1; i1<n_facets-1; ++i1)  // in xrange(i0+1,n_facets-1):
      {
        cut face1;
        faces->get_nth_plane(i1, face1);
        for( size_type i2=i1+1; i2<n_facets; ++i2)  // in xrange(i1+1,n_facets):
        {
          cut face2;
          faces->get_nth_plane(i2, face2);
          sg_mat3 m; //  = matrix.rec(facets[i0].n+facets[i1].n+facets[i2].n,(3,3))
          m.set_row(0, face0.n);
          m.set_row(1, face1.n);
          m.set_row(2, face2.n);
          int_type d  = m.determinant();
          if( d != 0 )
          {
            sg_mat3 c(  m.co_factor_matrix_transposed() ); //  / d );
            sg_vec3 b( -face0.c, -face1.c, -face2.c ); // matrix.col([-facets[i0].c,-facets[i1].c,-facets[i2].c])
            sg_vec3 iv = c * b;
            rvector3_t vertex( rvector3_t(iv) / rational_t(d) );
            if( this->is_inside_volume_only(vertex) ) // do not add if planes intersect outside of the asu
              result.insert(vertex);
          }
        }
      }
    }
    if( result.size()<4 )
      throw cctbx::error("Fewer than 4 vertices in asu ?");
  }

  void direct_space_asu::box_corners(rvector3_t &mn, rvector3_t &mx) const
  {
    set_rvector3_t vertices;
    volume_vertices(vertices);
    if( vertices.empty() )
      throw cctbx::error("No vertices in the asu");
    mn = *vertices.begin();
    mx = mn;
    for( set_rvector3_t::const_iterator v=vertices.begin(); v!=vertices.end(); ++v)
    {
      for( unsigned short i=0; i<3; ++i)
      {
        mn[i] = std::min(mn[i],(*v)[i]);
        mx[i] = std::max(mx[i],(*v)[i]);
      }
    }
  }


  inline scitbx::double3 conv_(const rvector3_t r)
  {
    return scitbx::double3( boost::rational_cast<double,int>(r[0]),
      boost::rational_cast<double,int>(r[1]), boost::rational_cast<double,int>(r[2]));
  }

  intersection_kind direct_space_asu::does_intersect(const scitbx::double3 &center,
    const scitbx::double3 &box) const
  {
    const signed char n_corners = 2;
    rvector3_t asu_box_r[n_corners];
    this->box_corners(asu_box_r[0], asu_box_r[1]);
    const scitbx::double3 asu_box[n_corners] = { conv_(asu_box_r[0]), conv_(asu_box_r[1])};
    CCTBX_ASSERT( scitbx::ge_all(asu_box[1], asu_box[0]) );
    const scitbx::double3 atom_box[n_corners] = { center - box, center + box };
    CCTBX_ASSERT( scitbx::ge_all(atom_box[1], atom_box[0]) );

    // need to test against enclosed_box for the 'fully' to work
    // if( scitbx::ge_all(atom_box[0], asu_enclosed_box[0]) &&
    // scitbx::le_all(atom_box[1], asu_enclosed_box[1]) )
    //  return fully;
    if( scitbx::ge_all(atom_box[1], asu_box[0]) && scitbx::le_all(atom_box[0], asu_box[1]) )
      return partially;

    return none;
  }

  void direct_space_asu::get_adjacent_cells(std::vector<scitbx::tiny3> &cells) const
  {
    cells.clear();
    BOOST_STATIC_ASSERT( sizeof(scitbx::tiny3::value_type) == sizeof(signed char) );
    for( signed char i=-1; i<=1; ++i)
      for( signed char j=-1; j<=1; ++j)
        for( signed char k=-1; k<=1; ++k)
          cells.push_back( scitbx::tiny3(i,j,k) );
  }

  void direct_space_asu::get_cells(std::vector<scitbx::tiny3> &cells) const
  {
    CCTBX_ASSERT( 0 ); // this function is not implemented yet
    cells.clear();
    cells.push_back( scitbx::tiny3(0,0,0) ); // TEMPORARY
    // this->get_adjacent_cells(cells);
  }


  rvector3_t direct_space_asu::move_inside(const cctbx::sgtbx::space_group &group, const rvector3_t &v) const
  {
    std::vector<scitbx::tiny3> cells;
    this->get_cells(cells);

    for(size_t i=0; i<group.order_z(); ++i)
    {
      const cctbx::sgtbx::rt_mx op = group(i);
      rvector3_t sv = op * v;
      sv -= scitbx::floor(sv);
      for(size_t icell=0; icell<cells.size(); ++icell)
      {
        rvector3_t sv_cell = sv + cells[icell];
        if( this->is_inside( sv_cell ) )
          return sv;
      }
    }
    CCTBX_ASSERT( 0 );
    return v;
  }

}}}

