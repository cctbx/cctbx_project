#include <xfel/mono_simulation/parameters.h>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>

xfel::parameter::parameter_array::parameter_array(
  int const& ndata,int const& kdim,farray d):
  parameters(d),
  gradients(ndata*kdim,scitbx::af::init_functor_null<double>()),
  curvatures(ndata*kdim,scitbx::af::init_functor_null<double>()),
  ndata(ndata),
  local_flag(false)
{
  SCITBX_ASSERT(d.size()==ndata*kdim);
}

xfel::parameter::parameter_array
xfel::parameter::organizer_base::register_array(
  std::string const& tag, int const& ndata, int const& kdim, farray d){
  _P[ tag ] = parameter_array(ndata,kdim,d);
  return _P[ tag ];
}

xfel::parameter::parameter_array
xfel::parameter::organizer_base::register_local_array(
  std::string const& tag, int const& ndata, int const& kdim, farray d){
  _P[ tag ] = parameter_array(ndata,kdim,d);
  _P[ tag ].local_flag=true;
  return _P[ tag ];
}

xfel::farray
xfel::parameter::organizer_base::as_x_array()const{
  farray result;
  typedef std::map<std::string,parameter_array>::const_iterator citer;
  for (citer M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    const double* bare_ptr = &(M->second.parameters.const_ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      result.push_back(bare_ptr[count]);
    }
  }
  return result;
}

xfel::farray
xfel::parameter::organizer_base::get_gradient_array()const{
  farray result;
  typedef std::map<std::string,parameter_array>::const_iterator citer;
  for (citer M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    const double* bare_ptr = &(M->second.gradients.const_ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      //std::cout<<M->first<<" "<<bare_ptr[count]<<std::endl;
      result.push_back(bare_ptr[count]);
    }
  }
  return result;
}

void
xfel::parameter::organizer_base::set_gradient_array(
  std::string const& aname, farray d)
{
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    if (aname == M->first) {
      double* bare_ptr = &(M->second.gradients.ref()[0]);

      SCITBX_ASSERT(d.size()==M->second.size());
      for (int ix = 0; ix < d.size(); ++ix){
        *bare_ptr++ = d[ix];
      }
    }
  }
}

xfel::farray
xfel::parameter::organizer_base::get_curvature_array()const{
  farray result;
  typedef std::map<std::string,parameter_array>::const_iterator citer;
  for (citer M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    const double* bare_ptr = &(M->second.curvatures.const_ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      result.push_back(bare_ptr[count]);
    }
  }
  // rest of this is for debugging only and can be removed.
  bool is_negative=false;
  for (int i = 0; i<result.size(); ++i){
    if (result[i] < 0.){is_negative=true;break;}
  }
  if (is_negative){
    for (citer M = _P.begin(); M != _P.end(); ++M){
      if (M->second.local_flag) {continue;}
      std::cout<<M->first<<std::endl;
      const double* bare_ptr = &(M->second.curvatures.const_ref()[0]);
      for (int count = 0; count < M->second.size(); ++count){
        std::cout<<"item "<<count<<" value "<<(bare_ptr[count])<<std::endl;
      }
    }
  }
  return result;
}

void
xfel::parameter::organizer_base::from_x_array(farray const& newx){
  farray result;
  const double* bare_x = &(newx.const_ref()[0]);
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    double* bare_ptr = &(M->second.parameters.ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      *bare_ptr++ = *bare_x++;
    }
  }
}

void
xfel::parameter::organizer_base::initialize_gradients_curvatures(){
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    M->second.gradients = farray(M->second.parameters.size());
    M->second.curvatures = farray(M->second.parameters.size());
  }
}

void
xfel::parameter::organizer_base::rezero_gradients_curvatures(){
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    double* grad_ptr = &(M->second.gradients.ref()[0]);
    double* curv_ptr = &(M->second.curvatures.ref()[0]);
    for (int i = 0; i < M->second.gradients.size(); ++i){
      *grad_ptr++ = 0.;
      *curv_ptr++ = 0.;
    }
  }
}

