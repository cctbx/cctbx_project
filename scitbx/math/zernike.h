#ifndef SCITBX_MATH_ZERNIKE_H
#define SCITBX_MATH_ZERNIKE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <scitbx/math/gamma.h>
#include <scitbx/vec3.h>

#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>



namespace scitbx { namespace math {
namespace zernike{

  using boost::math::spherical_harmonic;

  template <typename FloatType = double>
  class log_factorial_generator
  {
    public:
    /* Default constructor */
    log_factorial_generator(){}
    /* Basic constructor */
    log_factorial_generator(int const& n_max)
    :
    n_max_(n_max)
    {
      build_log_factorial_lookup();
    }

    // precompute log factorials
    void build_log_factorial_lookup()
    {
      data_.push_back( 0.0 );
      data_.push_back( 0.0 );
      exp_data_.push_back( 1.0 );
      exp_data_.push_back( 1.0 );
      FloatType tmp;

      for (int ii=2;ii<n_max_*2+5;ii++){
        tmp = scitbx::math::gamma::log_complete<FloatType>(ii+1);
        data_.push_back( tmp  );
        //exp_data_.push_back( tmp );
      }


    }

    FloatType log_fact(int n)
    {
      SCITBX_ASSERT( n>=0 );
      return(data_[n]);
    }

    FloatType fact(int n)
    {
      SCITBX_ASSERT( n>=0 );
      SCITBX_ASSERT( n<n_max_ );
      return(std::exp(data_[n]));
    }

    private:
    int n_max_;
    scitbx::af::shared<FloatType> data_, exp_data_;

  };

  //--------------------------------------------------------------
  //                DOUBLE INTEGER INDEX
  //--------------------------------------------------------------


  template <typename IntType = int>
  class double_integer_index
  {
    public:
    double_integer_index(){}
    double_integer_index(IntType const& n, IntType const& l)
    {
       n_ = n;
       l_ = l;
    }

    int n(){ return(n_); }

    int l(){ return(l_); }

    scitbx::af::tiny<IntType,3> doublet()
    {
       scitbx::af::tiny<int,2> result;
       result[0]=n_;
       result[1]=l_;
       return( result );
    }


   int n_, l_;

  };

  template <typename IntType=int>
  class double_integer_index_fast_less_than
  {
    public:
      //! This fast comparison function is implemented as operator().
      bool operator()(double_integer_index<IntType> const& nl1, double_integer_index<IntType> const& nl2) const
      {

         if (nl1.n_ < nl2.n_ ) return true;
         if (nl1.n_ > nl2.n_ ) return false;

         if (nl1.l_ < nl2.l_ ) return true;
         if (nl1.l_ > nl2.l_ ) return false;



         return false;
      }
  };




  //--------------------------------------------------------------
  //                TRIPLE INTEGER INDEX
  //--------------------------------------------------------------


  template <typename IntType = int>
  class nlm_index
  {
    public:
    nlm_index(){}

    nlm_index(IntType n, IntType l, IntType m)
    {
      n_ = n;
      l_ = l;
      m_ = m;
    }

    IntType operator[](std::size_t const& index)
    {
      SCITBX_ASSERT( index <=2 );
      SCITBX_ASSERT( index >=0 );
      if (index==0){ return(n_); }
      if (index==1){ return(l_); }
      if (index==2){ return(m_); }
    }

    int n(){ return(n_); }
    int l(){ return(n_); }
    int m(){ return(n_); }


    scitbx::af::tiny<IntType,3> triple()
    {
       scitbx::af::tiny<IntType,3> result;
      result[0] = n_ ;
      result[1] = l_ ;
      result[2] = m_ ;
      return result;
    }

    IntType n_,l_,m_;
  };


  template <typename IntType=int>
  class nlm_fast_less_than
  {
    public:
      //! This fast comparison function is implemented as operator().
      bool operator()(nlm_index<IntType> const& nlm1, nlm_index<IntType> const& nlm2) const
      {

         if ( nlm1.n_ < nlm2.n_ ) return true;
         if ( nlm1.n_ > nlm2.n_ ) return false;

         if ( nlm1.l_ < nlm2.l_ ) return true;
         if ( nlm1.l_ > nlm2.l_ ) return false;

         if ( nlm1.m_ < nlm2.m_ ) return true;
         if ( nlm1.m_ > nlm2.m_ ) return false;

         return false;
      }
  };




  //--------------------------------------------------------------
  //                ZERNIKE INDEX ARRAY
  //--------------------------------------------------------------


  template <typename FloatType = double>
  class nlm_array
  {

    typedef std::map< nlm_index<int>,
                      std::size_t,
                      nlm_fast_less_than<int> > lookup_map_type;

    typedef std::map< double_integer_index<int>,
                      std::size_t,
                      double_integer_index_fast_less_than<int> > nl_lookup_map_type;

    public:
    /* Default constructor */
    nlm_array() {}
    /* Basic constructor, sets all coefs to zero */
    nlm_array(int const& n_max)
    {
      SCITBX_ASSERT (n_max>0);
      n_max_=n_max;
      int count=0, n_duplicates=0, nl_count=0;
      for (int nn=0; nn<=n_max_; nn++){
        for (int ll=0;ll<=nn;ll++){
          // restriction on even / odd
          if (is_even( nn-ll )){


            scitbx::af::shared<int> tmp2;
            // make a lookup table for nl
            double_integer_index<int> this_nl(nn,ll);
            nl_.push_back( this_nl );
            nl_lookup_map_type::const_iterator l = nl_lookup_.find( this_nl );

            if ( l == nl_lookup_.end() ) { // not in list
                nl_lookup_[ this_nl ] = nl_count;
            }
            nl_count++;



            for (int mm=0;mm<=ll;mm++){
              tmp2.push_back( count ); // this index has a specific n,l pair

              nlm_index<int> this_nlm(nn,ll,mm);

              indices_.push_back( this_nlm );
              coefs_.push_back( 0.0 );
              //
              lookup_map_type::const_iterator l = nlm_lookup_.find( this_nlm );
              if ( l == nlm_lookup_.end() ) { // not in list
                nlm_lookup_[ this_nlm ] = count;
              } else {
                n_duplicates++;
              }
              SCITBX_ASSERT( find_nlm(nn,ll,mm)==count );
              count++;
              //std::cout << nn << " " << ll << " " << mm << std::endl;
            }
            nl_index_.push_back( tmp2 );

          }
          // if odd, leave it be
        }
      }
    }




    int find_nlm(int const& n, int const& l, int const& m)
    {
       nlm_index<> this_nlm(n,l,m);
       return(find_nlm(this_nlm));
    }

    int find_nlm(nlm_index<> const& this_nlm )
    {
       long nlm_location;
       lookup_map_type::const_iterator l = nlm_lookup_.find( this_nlm );
       if (l == nlm_lookup_.end()) {
         nlm_location = -1; // !!! negative if not found !!!
       }
       else {
         nlm_location = l->second;
       }
       return (nlm_location);
    }

    int find_nl(int const& n, int const& l)
    {
       double_integer_index<int> this_nl(n,l);
       return(find_nl(this_nl));
    }


    int find_nl(double_integer_index<int> const& this_nl )
    {
       int nl_location;
       nl_lookup_map_type::const_iterator l = nl_lookup_.find( this_nl );
       if (l == nl_lookup_.end()) {
         nl_location = -1; // !!! negative if not found !!!
       }
       else {
         nl_location = l->second;
       }
       return (nl_location);
    }


    bool set_coef(int const& n, int const&l, int const&m, FloatType const&x )
    {
       int this_index = find_nlm(n,l,m);
       if (this_index>-1){
         coefs_[ this_index ] = x;
         return(true);
       }
       return(false);
    }

    FloatType get_coef(int const& n, int const& l, int const& m)
    {
       int this_index = find_nlm(n,l,m);
       if (this_index>-1){
         return(coefs_[ this_index ]);
       }
       return(0.0);
    }


    scitbx::af::shared< int > select_on_nl(int const& n, int const& l)
    {
       scitbx::af::shared< int > selection;
       int this_index;
       this_index = find_nl( double_integer_index<int>(n,l) );
       return ( nl_index_[ this_index ] );
    }


    scitbx::af::shared< scitbx::af::tiny<int,3> > nlm()
    {
      scitbx::af::shared< scitbx::af::tiny<int,3> > result;
      for (int ii=0; ii<indices_.size();ii++){
        result.push_back( indices_[ii].triple() );
      }
      return( result );
    }


    scitbx::af::shared< scitbx::af::tiny<int,2> > nl()
    {
      return( nl_ );
    }


    scitbx::af::shared< FloatType > coefs()
    {
      return( coefs_ );
    }

    bool load_coefs(scitbx::af::shared< scitbx::af::tiny<int,3> > nlm,
                    scitbx::af::const_ref< FloatType > const& coef)
    {

       SCITBX_ASSERT(nlm.size()==coef.size());
       SCITBX_ASSERT(nlm.size()>0 );
       int this_one;
       bool found_it, global_find=true;
       for (int ii=0;ii<nlm.size();ii++){
         found_it = set_coef(nlm[ii][0],nlm[ii][1],nlm[ii][2],coef[ii]);
         if (!found_it){
           global_find=false;
         }
       }
       return(global_find);
    }


    private:
      bool is_even(std::size_t value)
      {
        std::size_t res;
        res = 2*(value/2);
        if (res == value){
          return(true);
        }
        return(false);
      }
      lookup_map_type nlm_lookup_;
      nl_lookup_map_type nl_lookup_;

      int n_max_;
      scitbx::af::shared< nlm_index<int> > indices_;
      scitbx::af::shared< FloatType > coefs_;
      scitbx::af::shared< double_integer_index<int> > nl_;
      scitbx::af::shared< scitbx::af::shared<int> > nl_index_;

  } ;







  //--------------------------------------------------------------
  //                ZERNIKE FUNCTIONS
  //--------------------------------------------------------------



  /*
   *  Radial Zernike polynome
   */

  template <typename FloatType = double>
  class zernike_radial
  {
    /* This class implements the radial part of a 3D zernike polynome for a specified nl value
     *
     */
    public:
    /* Default constructor */
    zernike_radial(){}
    /* Basic constructor */
    zernike_radial(int const& n, int const& l, log_factorial_generator<FloatType> const& lgf)
    :
    n_(n),
    l_(l),
    eps_(1e-18)
    {
      lgf_ = lgf;
      SCITBX_ASSERT( (n-l)%2==0 );
      compute_Nnlk();
    }

    void compute_Nnlk()
    {
       FloatType top1,top2,bottom1,bottom2,bottom3, bottom4, tmp1,tmp2,tmp3;
       tmp2 = std::pow(2.0, l_-n_);
       tmp1 = std::sqrt( 2.0*n_ + 3.0 );

       for (int k=0; k<=(n_-l_)/2; k++) {
         top1    =  lgf_.log_fact(  2*(n_-k)+1    )  ;
         top2    =  lgf_.log_fact(  (n_+l_)/2-k   )  ;
         bottom1 =  lgf_.log_fact(  (n_-l_)/2-k   ) ;
         bottom2 =  lgf_.log_fact(  (n_+l_-2*k+1) ) ;
         bottom3 =  lgf_.log_fact(  (n_-k)        ) ;
         bottom4 =  lgf_.log_fact(  k             )  ;

         tmp3 = top1+top2-bottom1-bottom2-bottom3-bottom4;
         if (tmp3>1e45){
           tmp3 = 1e45;
         }
         tmp3 = std::exp( tmp3 );
         tmp3 = tmp1*tmp2*tmp3*std::pow(-1.0,k);
         Nnlk_.push_back( tmp3 );
       }
    }

    FloatType f( FloatType const& r)
    {
      FloatType rr;
      if (r<=eps_){
        rr = eps_;
      } else {
        rr = r;
      }

      FloatType result,tmp;
      result=0.0;

      for (int kk=0;kk<Nnlk_.size();kk++){
        result+=std::pow(rr, n_ - 2*kk) * Nnlk_[kk];
      }
      return (result);
    }

    scitbx::af::shared< FloatType> f( scitbx::af::const_ref< FloatType> const& r)
    {
      scitbx::af::shared< FloatType> result;
      for (int ii=0;ii<r.size();ii++){
        result.push_back( f(r[ii]) );
      }
      return (result);
    }


    scitbx::af::shared< FloatType > Nnlk()
    {
       return( Nnlk_ );
    }

    private:



    int n_;
    int l_;
    scitbx::af::shared< FloatType > Nnlk_;
    log_factorial_generator<FloatType> lgf_;
    FloatType eps_;

  };



  /*
   * A single Zernike polynome of index n,l,m
   */
  template <typename FloatType = double>
  class zernike_polynome
  {
    public:
    /* Default constructor */
    zernike_polynome(){}
    /* Basic constructor */
    zernike_polynome(int const& n, int const& l, int const& m, zernike_radial<FloatType> const& rnl)
    :
    n_(n),
    l_(l),
    m_(m)
    {
      rnl_ = rnl;
    }
    std::complex<FloatType> f(FloatType const& r, FloatType const& t, FloatType const& p)
    {
      //std::cout << r << " " << t << " " << p << std::endl;
      std::complex<FloatType> result;
      FloatType tmp;
      tmp = rnl_.f(r);
      result = spherical_harmonic(l_, m_, t, p)*tmp;
      return(result);
    }
    private:
    int n_, l_, m_;
    zernike_radial<FloatType> rnl_;
  };






  /*
   *  Zernike function on a grid
   */

  template <typename FloatType = double>
  class zernike_grid
  {
    public:
    /* Default constructor */
    zernike_grid(){}
    /* Basic constructor */
    zernike_grid(int const& m, int const& n_max)
    :
    m_(m),
    n_max_(n_max),
    eps_(1e-12),
    nlm_(n_max_)
    {
      delta_ = 1.0 / (m-1);
      FloatType x,y,z,r,t,p;
      scitbx::vec3<FloatType> xyz, rtp;
      for (int ix=-m;ix<=m;ix++){
        for (int iy=-m;iy<=m;iy++){
          for (int iz=-m;iz<=m;iz++){
            x = ix*delta_;
            y = iy*delta_;
            z = iz*delta_;

            r = std::sqrt(x*x+y*y+z*z);
            if (r>eps_){
              t = std::acos(z/r);
              p = std::atan2(y,x);
            } else {
              t = 0.0;
              p = 0.0;
            }

            xyz[0]=x;
            xyz[1]=y;
            xyz[2]=z;

            rtp[0]=r;
            rtp[1]=t;
            rtp[2]=p;

            xyz_.push_back( xyz );
            rtp_.push_back( rtp );
          }
        }
      }

    }

    scitbx::af::shared< scitbx::vec3<FloatType> > xyz()
    {
      return( xyz_ );
    }


    void compute_basis_functions()
    {

      scitbx::af::shared< scitbx::af::tiny<int,3> > nlm;
      nlm = nlm_.nlm();

      scitbx::af::shared< scitbx::af::tiny<int,2> > nl;
      nl = nlm_.nl();

      for (int ii=0;ii<nlm.size();ii++){

      }

    }

    private:
      int m_, n_max_;
      FloatType delta_, eps_;
      scitbx::af::shared< scitbx::vec3<FloatType> > xyz_;
      scitbx::af::shared< scitbx::vec3<FloatType> > rtp_;
      scitbx::af::shared< std::vector<FloatType> > partial_data_;
      nlm_array<FloatType> nlm_;
  };




}}} // namespace scitbx::math::zernike

#endif // SCITBX_MATH_ZERNIKE_H
