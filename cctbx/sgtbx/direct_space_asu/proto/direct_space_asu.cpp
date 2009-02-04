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

}}}

