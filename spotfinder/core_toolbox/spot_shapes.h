#ifndef LIBDISTL_SHAPES_H
#define LIBDISTL_SHAPES_H

#include <boost/shared_ptr.hpp>
#include <scitbx/math/principal_axes_of_inertia.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/constants.h>
#include <iostream>

/* a pure-LABELIT source file; no distl equivalent */
namespace Distl {

typedef scitbx::math::principal_axes_of_inertia_2d<double> I_type_2d;

//calculate inertia data from experimental image
struct w_I_type: public I_type_2d {
      w_I_type(
        scitbx::af::const_ref<scitbx::vec2<double> > const& points,
        scitbx::af::const_ref<double> const& weights):
        I_type_2d(points,weights){}
      w_I_type():I_type_2d(){}
      inline virtual scitbx::vec2<double>
      center_of_mass(){
        return I_type_2d::center_of_mass();
      }
      inline virtual scitbx::vec2<double>
      eigenvector(const int& i){

        return scitbx::vec2<double>( eigensystem().vectors()(i,0) ,
          eigensystem().vectors()(i,1) );
      }
      inline virtual double
      eigenvalue(const int& i){
        return eigensystem().values()[i];
      }
      inline virtual ~w_I_type(){}
};

//contain unpickled inertia data
struct w2_I_type: public w_I_type {
  private:
      scitbx::vec2<double> center_of_mass_m;
      scitbx::vec2<double> values_m;
      scitbx::vec2<scitbx::vec2<double> > axis_m;
  public:
      w2_I_type(
        scitbx::vec2<double> const& center_of_mass,
        scitbx::vec2<double> const& values,
        scitbx::vec2<double> const& axis0,
        scitbx::vec2<double> const& axis1):
        center_of_mass_m(center_of_mass),
        values_m(values),
        axis_m(axis0,axis1)
        {
        }
      inline virtual scitbx::vec2<double>
      center_of_mass(){
        return center_of_mass_m;
      }
      inline virtual scitbx::vec2<double>
      eigenvector(const int& i){
        return axis_m[i];
      }
      inline virtual double
      eigenvalue(const int& i){
        return values_m[i];
      }
};

typedef boost::shared_ptr<w_I_type> ptr_I_type;

struct spot_shapes
{
    // instead of spot_shapes owning an I_type instance, it has a pointer
    // to a polymorphic wrapper type containing inertia data.  This complex
    // structure is intended to allow access to the spot inertia data
    // both before pickling (w_I_type) and after unpickling (w2_I_type).

    ptr_I_type model_m;
    double total_mass;

    spot_shapes():total_mass(0.0){} //default constructor

    inline void
    model_ellipse(scitbx::af::const_ref<scitbx::vec2<double> > const& points,
                scitbx::af::const_ref<double> const& weights) {
      total_mass = 0.0;
      for(std::size_t i_p=0;i_p<weights.size();i_p++) {
          total_mass += weights[i_p];
      }
      model_m = ptr_I_type(new w_I_type(points,weights));
      //show_axes();
    }

    inline void
    ss_setstate(scitbx::vec2<double> const& center_of_mass,
        scitbx::vec2<double> const& values,
        scitbx::vec2<double> const& axis0,
        scitbx::vec2<double> const& axis1,
        double const& mass) {
      total_mass = mass;
      model_m = ptr_I_type(new w2_I_type(
      center_of_mass,values,axis0,axis1));
    }

    inline scitbx::vec2<double>
    model_center() { return model_m->center_of_mass(); }

    inline double
    ctr_mass_x() const { return model_m->center_of_mass()[0]; }

    inline double
    ctr_mass_y() const { return model_m->center_of_mass()[1]; }

    inline scitbx::vec2<double>
      eigenvector(const int& i) const {
        return model_m->eigenvector(i);
      }

    inline double
      eigenvalue(const int& i){
        return model_m->eigenvalue(i);
      }

    inline double
    model_eccentricity() const {
      double b_squared = model_m->eigenvalue(1);
      double a_squared = model_m->eigenvalue(0);
      return std::sqrt( 1.0 - ( b_squared / a_squared ) );
    }

    inline bool
    com_valid(){ return total_mass > 0.; }

    inline void
    show_axes() {
      std::cout<<"mass "<<total_mass<<std::endl;
      std::cout<<"ctr mass "<<ctr_mass_x()<<" "<<ctr_mass_y()<<std::endl;
      std::cout<<"eigenvalues "<<model_m->eigenvalue(0)<<" "
                               <<model_m->eigenvalue(1)<<" "
      <<std::endl;
      std::cout<<"axis0 "<<model_m->eigenvector(0)[0]<<" "
                         <<model_m->eigenvector(0)[1]<<" "
      <<std::endl;
      std::cout<<"axis1 "<<model_m->eigenvector(1)[0]<<" "
                         <<model_m->eigenvector(1)[1]<<" "
      <<std::endl;
    }

    //model the spot as a uniform-density ellipse;
    // inertia tensor of an ellipse is I11=0.25*a*a*M; I22=0.25*b*b*M
    //semi-major axis
    inline double a() const {
      return std::sqrt(4.0*model_m->eigenvalue(0)/total_mass);
    }
    inline double b() const {
      //Problem:  when the points of a spot form a straight line segment,
      // the second eigenvalue is zero or a slightly negative number,
      // thus creating a floating point error.  What is reasonable in that
      // instance?  To set the eigenvalue to 1?
      return std::sqrt(4.0*std::max(1.,model_m->eigenvalue(1))/total_mass);
    }
};

}
#endif
