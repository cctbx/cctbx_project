#ifndef CCTBX_CRYSTL_ORIENT_H
#define CCTBX_CRYSTL_ORIENT_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/change_of_basis_op.h>

namespace cctbx {

  //! Shorthand for default vec3 type in orientation toolbox.
  typedef scitbx::vec3<double> oc_vec3;
  //! Shorthand for default mat3 type in orientation toolbox.
  typedef scitbx::mat3<double> oc_mat3;
  //! Shorthand for default vector list in orientation toolbox.
  typedef scitbx::af::shared<oc_vec3> oc_vec3_list;

  enum basis_type { direct=0, reciprocal=1 };

  //! Rotation of a vector
  /*!
      Right-hand rotation of a vector (first arg) about a laboratory unit vector
      (second arg) through an angle expressed in radians (third arg).
   */
  template <typename FloatType>
  oc_vec3
  vector_rotate_thru(oc_vec3 const& vec, oc_vec3 const& unit,
                     FloatType const& angle){
    //
    //Still need to add assertion that unit is indeed a unit vector
    //
    return unit*(unit*vec)
         +(vec-(unit*(unit*vec)))*std::cos(angle)
         -(vec.cross(unit))*std::sin(angle);
  }

  //! Class for the handling of crystal orientation matrix
  /*! All dimensions are in Angstroms (direct space) and inverse
      Angstroms (reciprocal space).
      <p>
      The reciprocal space matrix is a generalization of the fractionalization
      matrix:
      <pre>
                   / A*x B*x C*x \
      Matrix A* =  | A*y B*y C*y |
                   \ A*z B*z C*z /
      </pre>
      The direct space matrix is a generalization of the orthogonalization
      matrix:
      <pre>
                   / Ax  Ay  Az \
      Matrix A  =  | Bx  By  Bz |
                   \ Cx  Cy  Cz /
      </pre>
  */
  class crystal_orientation
  {
    public:
      //! Default constructor. Some data members are not initialized!
      /*!
          Not available in Python.
       */
      crystal_orientation(){}

      //! Constructor using a direct or reciprocal space matrix.
      /*!
          basis type flag:  false=direct space; true=reciprocal space
       */
      explicit
      crystal_orientation(oc_mat3 const&, bool const&);

      //! Access to direct space unit cell.
      cctbx::uctbx::unit_cell
      unit_cell() const;

      //! Access to reciprocal space unit cell.
      cctbx::uctbx::unit_cell
      unit_cell_inverse() const;

      //! Access to direct space orientation matrix.
      oc_mat3
      direct_matrix() const;

      //! Access to reciprocal space orientation matrix.
      oc_mat3
      reciprocal_matrix() const;

      //! Change of basis.  Return value gives mutated orientation.
      /*!
          Input is the change of basis operator
          for transforming direct space vectors from the current
          setting to a reference setting.
       */
      crystal_orientation
      change_basis(cctbx::sgtbx::change_of_basis_op const&) const;

      //! Change of basis.  Return value gives mutated orientation.
      /*!
          Input is a rotation matrix, assumed to be the matrix
          for transforming direct space vectors from the current
          setting to a reference setting.
       */
      crystal_orientation
      change_basis(oc_mat3 const&) const;

      //! Rotation mutator
      /*!
          Right-hand rotation of the crystal about a laboratory vector
          through an angle expressed in radians.
       */
      template <typename FloatType>
      crystal_orientation
      rotate_thru(
        oc_vec3 const& unit_axis,
        FloatType const& angle) const
      {
        //
        //still need to add an assertion that unit is a unit vector
        //
        oc_vec3_list abc; //astar,bstar,cstar
        for (int i=0; i<3; ++i){
           abc.push_back( Astar_.get_column(i) );
        }

        oc_vec3_list newabc;
        for (int i=0; i<3; ++i){
           newabc.push_back( vector_rotate_thru(abc[i],unit_axis,angle) );
        }

        oc_mat3 new_recip_matrix;
        for (int i=0; i<3; ++i){
           new_recip_matrix.set_column( i, newabc[i] );
        }

        return crystal_orientation(new_recip_matrix, cctbx::reciprocal);
      }

      //! Simple measure for the similarity of two orientatons.
      /*! The result is the mean of the squared differences between
          basis vectors. The basis vectors are taken in direct space.
          DEPRECATED 7/19/2009.
       */
      inline double
      direct_mean_square_difference(crystal_orientation const& other) const
      {
        return cctbx::uctbx::mean_square_difference(
                 direct_matrix(),other.direct_matrix())/3;
      }

      //! Simple measure for the similarity of two orientatons.
      /*! The result is a sum of mock Z-scores for the three basis vectors of
          the other orientation.  It assumes that the direct space vectors
          of this crystal orientation are distributed normally with a standard
          deviation of 1% of the basis vector length.
          Introduced 8/13/2009.
       */
      inline double
      difference_Z_score(crystal_orientation const& other) const
      {
        oc_mat3 diff = (direct_matrix() - other.direct_matrix());
        double Z = 0.0;
        for (int idx = 0; idx<3; ++idx){
          oc_vec3 this_basis_vector(direct_matrix().get_row(idx));
          oc_vec3 diff_vector(diff.get_row(idx));
          Z += diff_vector.length()/(0.01*this_basis_vector.length());
        }
        return Z;
      }

      //! Best change of basis for superimposing onto another orientation.
      /*! Used for aligning two orientations determined independently from the same
          crystal sample. The rotation part of an inverse cb_op is returned.
       */
      oc_mat3
      best_similarity_transformation(crystal_orientation const& other,
       double const& fractional_length_tolerance,
       int unimodular_generator_range=1) const;

    protected:
      //! Internal representation of the reciprocal matrix
      oc_mat3 Astar_;

  }; //class crystal_orientation
} //namespace cctbx

#endif // CCTBX_CRYSTL_ORIENT_H
