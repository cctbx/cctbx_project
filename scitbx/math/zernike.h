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
#include <scitbx/vec2.h>

#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>

using scitbx::constants::two_pi;

namespace scitbx { namespace math {
namespace zernike{

  using boost::math::spherical_harmonic;



  //--------------------------------------------------------------
  //           LOG FACTORIAL PRECOMPUTED LOOKUP TABLES
  //--------------------------------------------------------------


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

  template <typename IntType = int>
  class lm_array
  {

    typedef std::map< double_integer_index<int>,
                      std::size_t,
                      double_integer_index_fast_less_than<int> > lm_lookup_map_type;

    public:
    /* Default constructor */
    lm_array() {}
    /* Basic constructor, sets all coefs to zero */
    lm_array(int const& l_max)
    {
      SCITBX_ASSERT (l_max>0);
      l_max_=l_max;
      int count=0, n_duplicates=0, lm_count=0;
      for (int ll=0; ll<=l_max_; ll++){
        for (int mm=-ll;mm<=ll;mm++){
            scitbx::af::shared<int> tmp2;
            // make a lookup table for lm
            double_integer_index<int> this_lm(ll,mm);
            lm_.push_back( this_lm );
            lm_lookup_map_type::const_iterator l = lm_lookup_.find( this_lm );

            if ( l == lm_lookup_.end() ) { // not in list
                lm_lookup_[ this_lm ] = lm_count;
            }
            lm_count++;
        }
      }
    }

    int find_lm(int const& l, int const& m)
    {
       double_integer_index<int> this_lm(l,m);
       return(find_lm(this_lm));
    }


    int find_lm(double_integer_index<int> const& this_lm )
    {
       int lm_location;
       lm_lookup_map_type::const_iterator l = lm_lookup_.find( this_lm );
       if (l == lm_lookup_.end()) {
         lm_location = -1; // !!! negative if not found !!!
       }
       else {
         lm_location = l->second;
       }
       return (lm_location);
    }

    scitbx::af::shared< scitbx::af::tiny<int,2> > lm()
    {
      scitbx::af::shared< scitbx::af::tiny<int,2> > result;
      for(int ii=0;ii<lm_.size();ii++){
        scitbx::af::tiny<int,2> tmp( lm_[ii].n_, lm_[ii].l_);
        result.push_back( tmp );
      }
      return( result );
    }

    scitbx::af::shared< double_integer_index<int> > some_indices_in_legendre_recursion_order(int const& this_m)
    {
      scitbx::af::shared< double_integer_index<int> > result;
      for (int ii=this_m;ii<l_max_;ii++){
        double_integer_index<int> this_lm(ii,this_m);
        result.push_back( this_lm );
      }
      return ( result );
    }

    scitbx::af::shared< int > lut_of_some_indices_in_legendre_recursion_order(int const& this_m)
    {
      scitbx::af::shared< int > result;
      for (int ii=this_m;ii<l_max_;ii++){
        double_integer_index<int> this_lm(ii,this_m);
        int l = find_lm(this_lm);
        result.push_back( l );
      }
      return ( result );
    }

    private:

      lm_lookup_map_type lm_lookup_;

      int l_max_;
      scitbx::af::shared< double_integer_index<int> > lm_;
      scitbx::af::shared< scitbx::af::shared<int> > lm_index_;

  } ;




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
  //                ZERNIKE INDEX ARRAY of the NORM
  //--------------------------------------------------------------


  template <typename FloatType = double>
  class nl_array
  {

     typedef std::map< double_integer_index<int>,
                      std::size_t,
                      double_integer_index_fast_less_than<int> > nl_lookup_map_type;

    public:
    /* Default constructor */
    nl_array() {}
    /* Basic constructor, sets all coefs to zero */
    nl_array(int const& n_max)
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
            coefs_.push_back( 0.0 );
            nl_lookup_map_type::const_iterator l = nl_lookup_.find( this_nl );

            if ( l == nl_lookup_.end() ) { // not in list
                nl_lookup_[ this_nl ] = nl_count;
            }
            nl_count++;

          }
          // if odd, leave it be
        }
      }
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


    bool set_coef(int const& n, int const&l, FloatType const&x )
    {
       int this_index = find_nl(n,l);
       if (this_index>-1){
         coefs_[ this_index ] = x;
         return(true);
       }
       return(false);
    }

    FloatType get_coef(int const& n, int const& l)
    {
       int this_index = find_nl(n,l);
       if (this_index>-1){
         return(coefs_[ this_index ]);
       }
       return(0.0);
    }


    scitbx::af::shared< scitbx::af::tiny<int,2> > nl()
    {
      scitbx::af::shared< scitbx::af::tiny<int,2> > result;
      for(int ii=0;ii<nl_.size();ii++){
        scitbx::af::tiny<int,2> tmp( nl_[ii].n_, nl_[ii].l_);
        result.push_back( tmp );
      }
      return( result );
    }


    scitbx::af::shared< FloatType > coefs()
    {
      return( coefs_ );
    }

    bool load_coefs(scitbx::af::shared< scitbx::af::tiny<int,2> > nl,
                    scitbx::af::const_ref< FloatType > const& coef)
    {

       SCITBX_ASSERT(nl.size()==coef.size());
       SCITBX_ASSERT(nl.size()>0 );
       int this_one;
       bool found_it, global_find=true;
       for (int ii=0;ii<nl.size();ii++){
         found_it = set_coef(nl[ii][0],nl[ii][1],coef[ii]);
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

      nl_lookup_map_type nl_lookup_;

      int n_max_;
      scitbx::af::shared< FloatType > coefs_;
      scitbx::af::shared< double_integer_index<int> > nl_;
      scitbx::af::shared< scitbx::af::shared<int> > nl_index_;

  } ;


  //---------------------------------------------------
  // same as above, but now with complex numbers
  //
  template <typename FloatType = double>
  class nl_complex_array
  {

     typedef std::map< double_integer_index<int>,
                      std::size_t,
                      double_integer_index_fast_less_than<int> > nl_lookup_map_type;

    public:
    /* Default constructor */
    nl_complex_array() {}
    /* Basic constructor, sets all coefs to zero */
    nl_complex_array(int const& n_max)
    {
      SCITBX_ASSERT (n_max>0);
      n_max_=n_max;
      int count=0, n_duplicates=0, nl_count=0;
      for (int nn=0; nn<=n_max_; nn++){
        for (int ll=0;ll<=nn;ll++){
          // restriction on even / odd
          if (is_even( nn-ll )){
            //scitbx::af::shared<int> tmp2;
            // make a lookup table for nl
            double_integer_index<int> this_nl(nn,ll);
            nl_.push_back( this_nl );
            coefs_.push_back( 0.0 );
            nl_lookup_map_type::const_iterator l = nl_lookup_.find( this_nl );
            if ( l == nl_lookup_.end() ) { // not in list
                nl_lookup_[ this_nl ] = nl_count;
            }
            nl_count++;


          }
          // if odd, leave it be
        }
      }
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


    bool set_coef(int const& n, int const&l, std::complex<FloatType> const&x )
    {
       int this_index = find_nl(n,l);
       if (this_index>-1){
         coefs_[ this_index ] = x;
         return(true);
       }
       return(false);
    }

    std::complex<FloatType> get_coef(int const& n, int const& l)
    {
       int this_index = find_nl(n,l);
       if (this_index>-1){
         return(coefs_[ this_index ]);
       }
       return(0.0);
    }


    scitbx::af::shared< scitbx::af::tiny<int,2> > nl()
    {
      scitbx::af::shared< scitbx::af::tiny<int,2> > result;
      for(int ii=0;ii<nl_.size();ii++){
        scitbx::af::tiny<int,2> tmp( nl_[ii].n_, nl_[ii].l_);
        result.push_back( tmp );
      }
      return( result );
    }


    scitbx::af::shared< std::complex<FloatType> > coefs()
    {
      return( coefs_ );
    }

    bool load_coefs(scitbx::af::shared< scitbx::af::tiny<int,2> > nl,
                    scitbx::af::const_ref< std::complex<FloatType> > const& coef)
    {

       SCITBX_ASSERT(nl.size()==coef.size());
       SCITBX_ASSERT(nl.size()>0 );
       int this_one;
       bool found_it, global_find=true;
       for (int ii=0;ii<nl.size();ii++){
         found_it = set_coef(nl[ii][0],nl[ii][1],coef[ii]);
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

      nl_lookup_map_type nl_lookup_;

      int n_max_;
      scitbx::af::shared< std::complex<FloatType> > coefs_;
      scitbx::af::shared< double_integer_index<int> > nl_;
      scitbx::af::shared< scitbx::af::shared<int> > nl_index_;

  } ;





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
            //std::cout << nn << " " << ll << " " << std::endl;
            if ( l == nl_lookup_.end() ) { // not in list
                nl_lookup_[ this_nl ] = nl_count;
            }
            nl_count++;

            for (int mm=-ll;mm<=ll;mm++){
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

    bool set_coef(int const& n, int const&l, int const&m, std::complex<FloatType> const&x )
    {
       int this_index = find_nlm(n,l,m);
       if (this_index>-1){
         coefs_[ this_index ] = x;
         return(true);
       }
       return(false);
    }

    std::complex<FloatType> get_coef(int const& n, int const& l, int const& m)
    {
       int this_index = find_nlm(n,l,m);
       if (this_index>-1){
         return(coefs_[ this_index ]);
       }
       std::complex<FloatType> tmp(0.0);
       return(tmp);
    }

    scitbx::af::shared< int > select_on_nl(int const& n, int const& l)
    {
       scitbx::af::shared< int > selection;
       int this_index;
       this_index = find_nl( double_integer_index<int>(n,l) );
       return ( nl_index_[ this_index ] );
    }

    scitbx::af::shared< scitbx::af::shared<int> > nl_indices()
    {
      return(nl_index_);
    }

    scitbx::af::shared< scitbx::af::tiny<int,3> > nlm()
    {
      scitbx::af::shared< scitbx::af::tiny<int,3> > result;
      for (int ii=0; ii<indices_.size();ii++){
        result.push_back( indices_[ii].triple() );
      }
      return( result );
    }

    scitbx::af::shared< double_integer_index<int> > nl()
    {
      return( nl_ );
    }


    scitbx::af::shared< std::complex<FloatType> > coefs()
    {
      return( coefs_ );
    }


    bool load_coefs(scitbx::af::shared< scitbx::af::tiny<int,3> > nlm,
                    scitbx::af::const_ref< std::complex<FloatType> > const& coef)
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
      scitbx::af::shared< std::complex<FloatType> > coefs_;
      scitbx::af::shared< double_integer_index<int> > nl_;
      scitbx::af::shared< scitbx::af::shared<int> > nl_index_;

  } ;


  //--------------------------------------------------------------
  //                2D ZERNIKE Coefficient ARRAY B_nmk
  //--------------------------------------------------------------


  template <typename FloatType = double>
  class nmk_array
  {

    typedef std::map< nlm_index<int>,
                      std::size_t,
                      nlm_fast_less_than<int> > lookup_map_type;

    typedef std::map< double_integer_index<int>,
                      std::size_t,
                      double_integer_index_fast_less_than<int> > nl_lookup_map_type;

    public:
    /* Default constructor */
    nmk_array() {}
    /* Basic constructor, sets all coefs to zero */
    nmk_array(int const& n_max)
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
            //std::cout << nn << " " << ll << " " << std::endl;
            if ( l == nl_lookup_.end() ) { // not in list
                nl_lookup_[ this_nl ] = nl_count;
            }
            nl_count++;

            for (int kk=0;kk<=nn;kk++){
              if(is_even(nn-kk)) {
                tmp2.push_back( count ); // this index has a specific n,l pair

                nlm_index<int> this_nlm(nn,ll,kk);

                indices_.push_back( this_nlm );
                coefs_.push_back( 0.0 );
                //
                lookup_map_type::const_iterator l = nlm_lookup_.find( this_nlm );
                if ( l == nlm_lookup_.end() ) { // not in list
                  nlm_lookup_[ this_nlm ] = count;
                } else {
                n_duplicates++;
                }
                SCITBX_ASSERT( find_nlm(nn,ll,kk)==count );
                count++;
              //std::cout << nn << " " << ll << " " << kk << std::endl;
              }
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
       //std::cout<< nlm_location<<std::endl;
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
       FloatType tmp(0.0);
       return(tmp);
    }

    scitbx::af::shared< int > select_on_nl(int const& n, int const& l)
    {
       scitbx::af::shared< int > selection;
       int this_index;
       this_index = find_nl( double_integer_index<int>(n,l) );
       return ( nl_index_[ this_index ] );
    }

    scitbx::af::shared< scitbx::af::shared<int> > nl_indices()
    {
      return(nl_index_);
    }

    scitbx::af::shared< scitbx::af::tiny<int,3> > nlm()
    {
      scitbx::af::shared< scitbx::af::tiny<int,3> > result;
      for (int ii=0; ii<indices_.size();ii++){
        result.push_back( indices_[ii].triple() );
      }
      return( result );
    }

    scitbx::af::shared< double_integer_index<int> > nl()
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
  //                NOT-SO-SLOW SPHERICAL HARMONICS
  //--------------------------------------------------------------

  template <typename FloatType = double>
  class nss_spherical_harmonics
  {
    public:
    nss_spherical_harmonics();
    nss_spherical_harmonics(int const& max_l, int const& mangle, log_factorial_generator<FloatType> const& lgf)
    :
    max_l_(max_l),
    lm_engine_(max_l),
    mangle_(mangle),
    eps_(1e-18)
    {
      lgf_ = lgf;
      lm_indices_ = lm_engine_.lm();
      dtor_=3.14159265358979323846264338327950288/180.0;
      dangle_ = dtor_*360.0/mangle_;
      // first build precomputed sin and cos lookup tables
      // linear interpolation doesn't seem to slow things down much but does increase accurarcy of results.
      // it is however good to make sure that we use enough points to populate the table
      SCITBX_ASSERT( mangle> 3999 );
      dleg_ = 2.0/(mangle_-1.0);
      for(int ii=0;ii<mangle_;ii++){
        angles_.push_back( dangle_*ii );
        cos_.push_back( std::cos(dangle_*ii) );
        sin_.push_back( std::sin(dangle_*ii) );
        legendre_argument_.push_back( -1.0+ii*dleg_ );
      }

      // Now we need to build Legendre polynome lookup tables
      // again, we will not interpolate, but rely on nearest neighbour being okai for things up to order 30 or so
      //
      // We will build the lookuptables as using the recursion given by this relation
      //
      // P_{l+1,m} = [ (2*l+1)*x*P_{l,m} - (l-m+1)*P_{l-1,m} ] / [ l - m + 1]
      //
      //
      // we therefor have to compute
      // (m,m), (m+1,m) for each value of l
      // i.e, we do this
      // l,m = (0,0), (1,0) to get (2,0), (3,0), (4,0) etc etc
      // l,m = (1,1), (2,1) to get (3,1), (4,1), (5,1) etc etc
      // l,m = (2,2), (3,2) to get (4,2), (5,2), (6,2) etc etc
      // l,m = (3,3), (4,3) to get (5,3), (6,3), (7,3) etc etc
      // l,m = (4,4), (5,4) to get (6,4), (7,4), (8,4) etc etc
      // etc etc etc etc etc etc
      //
      //
      //
      for (int ii=0;ii<mangle_;ii++){
        std::vector< FloatType > plm( lm_indices_.size(), 0 ); // here we store the legendre functions
        FloatType x = legendre_argument_[ ii ];
        scitbx::af::shared<int> indices;
        for (int this_m=0; this_m < max_l_; this_m++){
          // get the indices please for this set of indices
          indices = lm_engine_.lut_of_some_indices_in_legendre_recursion_order( this_m );
          FloatType a,b,c;
          int this_l;
          this_l = lm_indices_[ indices[0] ][0] ; // the n/l l/m correspondence is annoying but a consequence of naming decisions made earlier
          //this_m = lm_indices_[ indices[0] ][1] ;
          a = boost::math::legendre_p(this_l, this_m, x);
          plm[ indices[0] ] = a;
          if (this_m>0){
                int neg_lm = lm_engine_.find_lm( this_l , -this_m );
                plm[ neg_lm ] = neg_legendre(this_l, this_m, b);
          }
          // make sure we can go further!
          if (indices.size() > 1){
            this_l = lm_indices_[ indices[1] ][0];
            b = boost::math::legendre_p(this_l, this_m, x);
            plm[ indices[1] ] = b;
            if (this_m>0){
                int neg_lm = lm_engine_.find_lm( this_l , -this_m );
                plm[ neg_lm ] = neg_legendre(this_l, this_m, b);
              }
            for (int uu=2;uu<indices.size();uu++){
              this_l = lm_indices_[ indices[uu] ][0] ;
              c = boost::math::legendre_next( this_l-1, this_m, x, b, a);
              plm[ indices[uu] ] = c;

              //---------------------------
              // now do (l,-m)
              if (this_m>0){
                int neg_lm = lm_engine_.find_lm( this_l , -this_m );
                plm[ neg_lm ] = neg_legendre(this_l, this_m, c);
              }
              //---------------------------
              //
              // swap vars
              a = b;
              b = c;
            }
          }
        }
        legendre_.push_back( plm );
      }

      // please compute the normalisation coefficients
      for (int ii=0;ii<lm_indices_.size();ii++){
        int l,m;
        l=lm_indices_[ii][0];
        m=lm_indices_[ii][1];
        FloatType tmp;
        tmp = std::sqrt( ( (2.0*l+1.0)/(4.0*3.14159265358979323846264338327950288) )*std::exp( lgf_.log_fact(l-m) - lgf_.log_fact(l+m) ) );
        norma_lm_.push_back( tmp );
      }


    }

    FloatType neg_legendre(int l, int m, FloatType plm)
    {
      FloatType tmp;
      tmp = std::pow(-1.0,m)*std::exp( lgf_.log_fact(l-m) - lgf_.log_fact(l+m) );
      return( plm*tmp );
    }

    scitbx::af::shared< FloatType > legendre_lm_pc(int const& l, int const& m)
    {
       // please find the index where to find this polynome
       int index = lm_engine_.find_lm(l,m);
       scitbx::af::shared< FloatType > result;
       for (int ii=0;ii<mangle_;ii++){
         result.push_back( legendre_[ii][index] );
       }
       return( result );
    }

    scitbx::af::shared< FloatType > legendre_lm(int const& l, int const& m)
    {
       // please find the index where to find this polynome
       int index = lm_engine_.find_lm(l,m);
       scitbx::af::shared< FloatType > result;
       for (int ii=0;ii<mangle_;ii++){
         result.push_back( boost::math::legendre_p(l,m,legendre_argument_[ii]) );
       }
       return( result );
    }

    std::complex<FloatType> spherical_harmonic_pc(int const& l, int const& m, FloatType const& theta, FloatType const& phi)
    {
       FloatType frac_theta, frac_phi, mphi, om_frac_phi, om_frac_theta, frac_cos, om_frac_cos;
       mphi = m*phi;
       int theta_index_l= static_cast<int>(theta/dangle_)%mangle_;
       int phi_index_l  = static_cast<int>(mphi/dangle_)%mangle_;
       int phi_index    = static_cast<int>(mphi/dangle_);
       int theta_index_h= (theta_index_l+1)%mangle_;
       int phi_index_h  = (phi_index_l+1)%mangle_;

       frac_theta    = (theta-theta_index_l*dangle_)/dangle_;
       om_frac_theta = 1.0-frac_theta;
       frac_phi      = (mphi-phi_index*dangle_)/dangle_;
       om_frac_phi   = 1.0-frac_phi;

       FloatType cos_phi, sin_phi, cos_theta;
       cos_phi   = cos_[phi_index_l  ]+(cos_[phi_index_h  ]-cos_[phi_index_l  ])*frac_phi   ;
       sin_phi   = sin_[phi_index_l  ]+(sin_[phi_index_h  ]-sin_[phi_index_l  ])*frac_phi   ;

       cos_theta = cos_[theta_index_l]+ (cos_[theta_index_h]-cos_[theta_index_l])*frac_theta;


       int i_cos_theta_l = static_cast<int>( (cos_theta+1.0)/dleg_ );
       int i_cos_theta_h;

       i_cos_theta_h = i_cos_theta_l + 1;
       frac_cos      = ( (cos_theta+1.0)-(i_cos_theta_l*dleg_) )/dleg_;
       om_frac_cos = 1.0-frac_cos;
       if (i_cos_theta_h > mangle_-1){
        i_cos_theta_h = mangle_-1;
       }


       FloatType part1, part2, part3;
       int lm_index = lm_engine_.find_lm(l,m);
       part1 = legendre_[i_cos_theta_l][lm_index]+(legendre_[i_cos_theta_h][lm_index]-legendre_[i_cos_theta_l][lm_index])*frac_cos;
       part1 = part1*norma_lm_[lm_index];
       part2 = cos_phi*part1;
       part3 = sin_phi*part1;

       std::complex<FloatType> result(part2,part3);
       return (result);
    }

    std::complex<FloatType> spherical_harmonic_direct(int const& l, int const& m, FloatType const& theta, FloatType const& phi)
    {
       std::complex<FloatType> result;
       result = spherical_harmonic(l, m, theta, phi);
       return (result);
    }


    private:

    std::vector< FloatType > cos_;
    std::vector< FloatType > sin_;
    std::vector< FloatType > angles_;
    std::vector< FloatType > legendre_argument_;
    std::vector< std::vector< FloatType> > legendre_;
    lm_array<int> lm_engine_;
    scitbx::af::shared< scitbx::af::tiny<int,2> > lm_indices_;
    scitbx::af::shared< FloatType > norma_lm_;
    int max_l_, mangle_;
    FloatType dangle_, dleg_, dtor_, eps_;
    log_factorial_generator<FloatType> lgf_;

  };



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
      n_terms_=Nnlk_.size();
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

      for (int kk=0;kk<n_terms_;kk++){
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

    int n(){ return(n_); }
    int l(){ return(l_); }


    private:
    int n_;
    int l_;
    int n_terms_;
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
      SCITBX_ASSERT( rnl_.n() == n_ );
      SCITBX_ASSERT( rnl_.l() == l_ );

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

    std::complex<FloatType> real_f(FloatType const& r, FloatType const& t, FloatType const& p) // the return value is not real, but it allows for faster computation of a real valued function
    {
      //std::cout << r << " " << t << " " << p << std::endl;
      std::complex<FloatType> result;
      result = f(r,t,p);
      if (m_!=0){ result = result*2.0; } // double it up so that we do not have to sum over -m
      return(result);
    }

    int n(){ return( n_ ); }
    int l(){ return( l_ ); }
    int m(){ return( m_ ); }


    private:
    int n_, l_, m_;
    zernike_radial<FloatType> rnl_;
  };


/*

  template <typename FloatType = double>
  class nss_zernike_polynome_engine
  {
    public:
    // Default constructor
    nss_zernike_polynome_engine(){}
    // Basic constructor
    nss_zernike_polynome_engine(int const& max_n, log_factorial_generator<FloatType> const& lgf)
    :
    n_max_(n_max),
    {
      rnl_ = rnl;
      SCITBX_ASSERT( rnl_.n() == n_ );
      SCITBX_ASSERT( rnl_.l() == l_ );

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

    std::complex<FloatType> real_f(FloatType const& r, FloatType const& t, FloatType const& p) // the return value is not real, but it allows for faster computation of a real valued function
    {
      //std::cout << r << " " << t << " " << p << std::endl;
      std::complex<FloatType> result;
      result = f(r,t,p);
      if (m_!=0){ result = result*2.0; } // double it up so that we do not have to sum over -m
      return(result);
    }
    private:
    int n_max_;

    zernike_radial<FloatType> rnl_;
  };

*/



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
    zernike_grid(int const& m, int const& n_max, bool hex=false)
    :
    m_(m),          // length of cube  = m*2+1
    n_max_(n_max),  // order of expansion
    hex_(hex),      // grid type, cubic or hexagonal
    eps_(1e-12),    // epsilon
    nlm_(n_max_),   // nlm index
    lgf_(n_max_*2+5)// factorial engine
    {
      delta_ = 1.0 / (m-1);
      FloatType x,y,z,r,t,p;
      scitbx::vec3<FloatType> xyz, rtp;
      scitbx::vec3<int> ijk;
      // get all indices please
      scitbx::af::shared< scitbx::af::tiny<int,3> > nlm_ind = nlm_.nlm();
      // make a log factorial lookup table
      log_factorial_generator<FloatType> lgf;
      scitbx::af::shared< zernike_radial<FloatType> > zr;
      scitbx::af::shared< zernike_polynome<FloatType> > zp;

      // build radial functions
      scitbx::af::shared< double_integer_index<int> > nl = nlm_.nl();
      for (int ii=0;ii<nl.size();ii++){
        zr.push_back( zernike_radial<FloatType>( nl[ii].n(), nl[ii].l(), lgf_ ) );
      }

      // build array of zernike polynomes
      for (int ii=0;ii<nlm_ind.size();ii++){
        int n,m,l,nl_index;
        n=nlm_ind[ii][0]; l=nlm_ind[ii][1]; m=nlm_ind[ii][2];
        nl_index = nlm_.find_nl(n,l);
        zernike_polynome<FloatType> this_zp(n,l,m,zr[nl_index]);
        zp.push_back( this_zp );
      }

      build_grid();
      int np_total = xyz_.size();
      for(int i=0; i< np_total; i++) {
        r=rtp_[i][0];
        t=rtp_[i][1];
        p=rtp_[i][2];

        // for each point, make a place holder for the zernike basis function
        std::vector< std::complex<FloatType> > tmp_result;
        // loop over all indices nlm and precompute all coefficients
        if (r<=1.0){
          for (int ii=0;ii<zp.size();ii++){
            tmp_result.push_back( zp[ii].f(r,t,p) );
            }
        } else {
          std::complex<FloatType> tmp(0,0);
          tmp_result.push_back( tmp ); // when radius bigger then 1
        }
        partial_data_.push_back( tmp_result );
      }
    }

    void build_grid()
    {
      FloatType x,y,z,r,t,p;
      scitbx::vec3<FloatType> xyz, rtp;
      if(hex_) {
        int yRow = 0;
        FloatType dr, dx, dy, dz, max_x, max_y, max_z;
        dr = 1.0/(2.0*m_);
        dx = 2*dr;
        dy = std::sqrt(3.0)*dr;
        dz = std::sqrt(6.0)*2.0/3.0*dr;
        max_x = 1.0;
        max_y = 1.0;
        max_z = 1.0;
        z = -1.0;
        bool is_plane_A = true; // A-B-A-B-A... planes

        while( z<=max_z ) {
          if(is_plane_A) {
            yRow = 0;
            x=-1.0;
            y=-1.0;
            while( y<=max_y) {
              while( x<max_x ) {
                xyz[0]=x;
                xyz[1]=y;
                xyz[2]=z;
                xyz_.push_back(xyz);
                x=x+dx;
              }
              yRow = yRow + 1;
              y=y+dy;
              if( yRow % 2 == 1 ) x=dr;
              else x = 0.0;
              x=x-1.0;
            }
          }
          else{
            yRow = 0;
            x=-1.0 + dr;
            y=-1.0 + dy/3.0;
            while( y<=max_y) {
              while( x<max_x ) {
                xyz[0]=x;
                xyz[1]=y;
                xyz[2]=z;
                xyz_.push_back(xyz);
                x=x+dx;
              }
              yRow = yRow + 1;
              y=y+dy;
              if( yRow % 2 == 0 ) x=dr;
              else x = 0.0;
              x=x-1.0;
            }
          }
          is_plane_A = !is_plane_A; // next plane
          z = z + dz;
        }

      }
      else {
        for (int ix=-m_;ix<=m_;ix++){
         for (int iy=-m_;iy<=m_;iy++){
          for (int iz=-m_;iz<=m_;iz++){

            x = ix*delta_;
            y = iy*delta_;
            z = iz*delta_;

            xyz[0]=x;
            xyz[1]=y;
            xyz[2]=z;
            xyz_.push_back( xyz );
          }
         }
        }
      }

     int np_total = xyz_.size();
     for(int i=0;i<np_total;i++)
     {
        x=xyz_[i][0];
        y=xyz_[i][1];
        z=xyz_[i][2];

        r = std::sqrt(x*x+y*y+z*z);
        if (r>eps_){
         t = std::acos(z/r);
         p = std::atan2(y,x);
         } else {
           t = 0.0;
           p = 0.0;
         }

        rtp[0]=r;
        rtp[1]=t;
        rtp[2]=p;

        rtp_.push_back( rtp );
     }

      return;
    }

    scitbx::af::shared< std::complex<FloatType> > slow_moments(
                                 scitbx::af::const_ref< FloatType> const& image
                                 )
    {
      scitbx::af::shared< std::complex<FloatType> > result;
      scitbx::af::shared< scitbx::af::tiny<int,3> > nlm_ind = nlm_.nlm();
      // first init the accumulator array
      for (int ii=0;ii<nlm_ind.size();ii++){
        std::complex<FloatType> tmp(0.0,0.0);
        result.push_back( tmp );
      }
      // now loop over all coordinates
      // and compute whatever we need to compute
      SCITBX_ASSERT( image.size() == xyz_.size() );

      for (int ii=0;ii<xyz_.size();ii++){
        for (int jj=0;jj<partial_data_[ii].size();jj++){
          result[jj] += image[ii]*std::conj(partial_data_[ii][jj]);
        }
      }
      // normalize
      for (int ii=0;ii<result.size();ii++){
        std::cout << result[ii]/static_cast<FloatType>(xyz_.size()) << std::endl;
        result[ii] = result[ii]/static_cast<FloatType>(xyz_.size());
      }
      return(result);
    }




    scitbx::af::shared< scitbx::vec3<FloatType> > xyz()
    {
      return( xyz_ );
    }
    scitbx::af::shared< scitbx::vec3<FloatType> > rtp()
    {
      return(rtp_ );
    }



    bool load_coefs(scitbx::af::shared< scitbx::af::tiny<int,3> > these_nlm,
                    scitbx::af::const_ref< std::complex<FloatType> > const& coef)
    {
       return( nlm_.load_coefs(these_nlm,coef) );
    }


    scitbx::af::shared< std::complex<FloatType> > f()
    {
       scitbx::af::shared< scitbx::af::tiny<int,3> > nlm_ind = nlm_.nlm();
       scitbx::af::shared< std::complex<FloatType> > coefs = nlm_.coefs();
       scitbx::af::shared< std::complex<FloatType> > result;
       std::complex<FloatType> tmp;
       for (int gg=0;gg<partial_data_.size();gg++){
         std::complex<FloatType> tmp_result=0.0;
         //std::cout << xyz_[gg][0] << " " << xyz_[gg][1] << " " << xyz_[gg][2] << "  --->  ";
         for (int ii=0;ii<partial_data_[gg].size();ii++){
           //std::cout <<  "(" << nlm_ind[ii][0] << " " << nlm_ind[ii][1] << " " << nlm_ind[ii][2]
           //          <<  " " << coefs[ii] << " )  ";
           tmp = partial_data_[gg][ii]*coefs[ii];
           tmp_result += tmp;//std::real(tmp) - std::imag(tmp);
         }
         //std::cout << std::endl;
         result.push_back( tmp_result );
       }
       return( result );
    }

    scitbx::af::shared< FloatType > f_real()
    {
       scitbx::af::shared< scitbx::af::tiny<int,3> > nlm_ind = nlm_.nlm();
       scitbx::af::shared< std::complex<FloatType> > coefs = nlm_.coefs();
       scitbx::af::shared< FloatType > result;
       std::complex<FloatType> tmp;
       for (int gg=0;gg<partial_data_.size();gg++){
         std::complex<FloatType> tmp_result=0.0;
         for (int ii=0;ii<partial_data_[gg].size();ii++){
           tmp = partial_data_[gg][ii]*coefs[ii];
           tmp_result += tmp;
         }
         result.push_back( std::real(tmp_result) );
       }
       return( result );
    }



    scitbx::af::shared< std::complex<FloatType> > coefs()
    {
      return( nlm_.coefs() );
    }

    scitbx::af::shared< scitbx::af::tiny<int,3> > nlm()
    {
      return( nlm_.nlm() );
    }


    private:
      int m_, n_max_;
      bool hex_;
      FloatType delta_, eps_;
      scitbx::af::shared< scitbx::vec3<FloatType> > xyz_;
      scitbx::af::shared< scitbx::vec3<FloatType> > rtp_;
      scitbx::af::shared< scitbx::vec3<int> > ijk_;

      scitbx::af::shared< std::vector< std::complex<FloatType> > > partial_data_;
      nlm_array<FloatType> nlm_;
      log_factorial_generator<FloatType> lgf_;
      scitbx::af::shared< FloatType > model_;

  };



  //------------------------------------------
  // Here we implement 2D Zernike functions
  // They are usefull for 2D type problems
  //------------------------------------------


  /*
   *  Radial Zernike polynome, 2D
   */

  template <typename FloatType = double>
  class zernike_2d_radial
  {
    /* This class implements the radial part of a 2D zernike polynome for a specified nl value
     *
     */
    public:
    /* Default constructor */
    zernike_2d_radial(){}
    /* Basic constructor */
    zernike_2d_radial(int const& n, int const& l, log_factorial_generator<FloatType> const& lgf)
    :
    n_(n),
    l_(l),
    eps_(1e-18)
    {
      lgf_ = lgf;
      SCITBX_ASSERT( (n-l)%2==0 );
      compute_Nnlk();
      n_terms_=Nnlk_.size();
    }

    void compute_Nnlk()
    {
       FloatType top1,bottom1,bottom2,bottom3,tmp3;
       for (int k=0; k<=(n_-l_)/2; k++) {
         top1    =  lgf_.log_fact(  (n_-k)        )  ;
         bottom1 =  lgf_.log_fact(  (n_+l_)/2-k   ) ;
         bottom2 =  lgf_.log_fact(  (n_-l_)/2-k   ) ;
         bottom3 =  lgf_.log_fact(  k             )  ;
         tmp3 = top1-bottom1-bottom2-bottom3;
         if (tmp3>1e45){
           tmp3 = 1e45;
         }
         tmp3 = std::exp( tmp3 );
         tmp3 = tmp3*std::pow(-1.0,k);
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
      for (int kk=0;kk<n_terms_;kk++){
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

    int n(){ return(n_); }
    int l(){ return(l_); }


    private:
    int n_;
    int l_;
    int n_terms_;
    scitbx::af::shared< FloatType > Nnlk_;
    log_factorial_generator<FloatType> lgf_;
    FloatType eps_;
  };


  /*
   * Implementation of 2d zernike radial function using discrete cosine expansion
   */
  template <typename FloatType = double>
  class zernike_2d_radial_dc
  {
    public:
    /* Default constructor */
    zernike_2d_radial_dc() {}
    /* Basic Constructor */
    zernike_2d_radial_dc(int const& n, int const& l)
    :
    n_(n),
    n_terms_(n*2+1),
    n_plus_1_(static_cast<FloatType>(n+1)),
    l_(l),
    eps_(1e-18)
    {
      SCITBX_ASSERT( (n_-l_)/2*2 ==(n_-l_) );
      two_pi_ = scitbx::constants::pi*2.0;
      if(n_>0) {
        two_pi_over_nterm_ = two_pi_/static_cast<FloatType>(n_terms_);
        two_pi_over_nterm_l_ = two_pi_over_nterm_*static_cast<FloatType>(l_);
//        build_array();
      }
    }

    void build_array()
    {
      FloatType n_terms_float( n_terms_ );
      step_size_ = 1.0/n_terms_float;
      r_array_.reserve(n_terms_+1);
      f_r_array_.reserve(n_terms_+1);
      for(int i=0;i<=n_terms_;i++)
      {
        x_=i/n_terms_float;
        r_array_.push_back( x_ );
        f_r_array_.push_back( f_exact(x_) );
      }
      return;
    }

    FloatType f( FloatType r )
    {
      if(n_==0) return 1.0;
      if(r==1.0) return 1.0;
      return f_exact(r);
/*      indx_=int(r*n_terms_);
      //std::cout<<r_array_[indx_]<<" "<<r_array_[indx_+1]<<" "<<f_r_array_[indx_]<<" "<<f_r_array_[indx_+1]<<r<<std::endl;
      value_=(r-r_array_[indx_])*f_r_array_[indx_+1]-(r-r_array_[indx_+1])*f_r_array_[indx_];
      value_/=step_size_;
      return value_;
 */
    }

    FloatType f_exact(FloatType r)
    {
      if(n_==0) return 1.0;
      if(r==1.0) return 1.0;
      value_=0.0;
      for(int k=0;k<n_terms_;k++)
      {
        x_=r*std::cos(two_pi_over_nterm_*k);
        mu_ = std::acos(x_);
        temp_ = std::sin(mu_*(n_plus_1_))/std::sin(mu_)*std::cos(two_pi_over_nterm_l_*k);
        value_+=temp_;
      }
      return value_/n_terms_;
    }

    int n(){ return(n_); }
    int l(){ return(l_); }

    private:
    int n_;
    int l_;
    int n_terms_;
    int indx_;
    FloatType eps_, value_, two_pi_, two_pi_over_nterm_, two_pi_over_nterm_l_;
    FloatType x_, temp_,mu_, n_plus_1_;
    FloatType step_size_;
    scitbx::af::shared< FloatType > f_r_array_;
    scitbx::af::shared< FloatType > r_array_;
  };

  /*
   * A single Zernike 2d polynome of index n,l
   */
  template <typename FloatType = double>
  class zernike_2d_polynome
  {
    public:
    /* Default constructor */
    zernike_2d_polynome(){}
    /* Basic constructor */
    zernike_2d_polynome(int const& n, int const& l):
    rnl_(n,l),
    n_(n),
    l_(l)
    {
      //rnl_ = rnl;
      SCITBX_ASSERT( rnl_.n() == n_ );
      SCITBX_ASSERT( rnl_.l() == l_ );
    }
    std::complex<FloatType> f(FloatType const& r, FloatType const& t)
    {
      //std::cout << r << " " << t << " " << p << std::endl;
      std::complex<FloatType> result;
      FloatType tmp;
      tmp = rnl_.f(r);
      result = expo_part(l_, t)*tmp;
      return(result);
    }

    std::complex<FloatType> expo_part(int l,FloatType t)
    {
      return (  std::complex<FloatType>( std::cos(l*t), std::sin(l*t) ) );
    }

    int n(){ return( n_ ); }
    int l(){ return( l_ ); }


    private:
    int n_, l_;
    zernike_2d_radial_dc<FloatType> rnl_;
  };




  template <typename FloatType = double>
  class zernike_grid_2d
  {
    public:
    /* Default constructor */
    zernike_grid_2d(){}
    /* Basic constructor */
    zernike_grid_2d(int const& m, int const& n_max)
    :
    m_(m),          // length of cube  = m*2+1
    n_max_(n_max),  // order of expansion
    eps_(1e-12),    // epsilon
    nl_(n_max_),    // nl index
    lgf_(n_max_*2+5)// factorial engine
    {
      delta_ = 1.0 / (m);
      FloatType x,y,r,t;
      scitbx::vec2<FloatType> xy, rt;
      scitbx::vec2<int> ij;
      // get all indices please
      // scitbx::af::shared< scitbx::af::tiny<int,2> > nl_indices_only_ = nl_.nl();
      nl_indices_only_ = nl_.nl();
      // make a log factorial lookup table
      log_factorial_generator<FloatType> lgf;
      scitbx::af::shared< zernike_2d_polynome<FloatType> > zp;
      for (int ii=0;ii<nl_indices_only_.size();ii++){
        zp.push_back( zernike_2d_polynome<FloatType>( nl_indices_only_[ii][0],
                                                      nl_indices_only_[ii][1] ) );
      }

      build_grid();
      int np_total = xy_.size();
      for(int i=0; i< np_total; i++) {
        r=rt_[i][0];
        t=rt_[i][1];

        // for each point, make a place holder for the zernike basis function
        std::vector< std::complex<FloatType> > tmp_result;
        // loop over all indices nlm and precompute all coefficients
        if (r<=1.0){
          for (int ii=0;ii<zp.size();ii++){
            tmp_result.push_back( zp[ii].f(r,t) );
            }
        } else {
          std::complex<FloatType> tmp(0,0);
          tmp_result.push_back( tmp ); // when radius bigger then 1
        }
        partial_data_.push_back( tmp_result );
      }
    }

    void build_grid()
    {
      FloatType x,y,r,t;
      scitbx::vec2<FloatType> xy, rt;

        for (int ix=-m_;ix<=m_;ix++){
         for (int iy=-m_;iy<=m_;iy++){
            x = ix*delta_;
            y = iy*delta_;
            xy[0]=x;
            xy[1]=y;
            xy_.push_back( xy );
          }
         }



     int np_total = xy_.size();
     for(int i=0;i<np_total;i++)
     {
        x=xy_[i][0];
        y=xy_[i][1];
        r = std::sqrt(x*x+y*y);
        if (r>eps_){
         t = std::atan2(y,x);
         } else {
           t = 0.0;
         }

        rt[0]=r;
        rt[1]=t;
        rt_.push_back( rt );
     }

      return;
    }




    scitbx::af::shared< scitbx::vec2<FloatType> > xy()
    {
      return( xy_ );
    }
    scitbx::af::shared< scitbx::vec2<FloatType> > rt()
    {
      return(rt_ );
    }



    bool load_coefs(scitbx::af::shared< scitbx::af::tiny<int,2> > these_nl,
                    scitbx::af::const_ref< std::complex<FloatType> > const& coef)
    {
       return( nl_.load_coefs(these_nl,coef) );
    }


    scitbx::af::shared< std::complex<FloatType> > f()
    {
       //scitbx::af::shared< scitbx::af::tiny<int,2> > nl_ind = nl_.nl();
       scitbx::af::shared< std::complex<FloatType> > coefs = nl_.coefs();
       scitbx::af::shared< std::complex<FloatType> > result;
       std::complex<FloatType> tmp;
       for (int gg=0;gg<partial_data_.size();gg++){
         std::complex<FloatType> tmp_result=0.0;
         for (int ii=0;ii<partial_data_[gg].size();ii++){
           tmp = partial_data_[gg][ii]*coefs[ii];
           if (nl_indices_only_[ii][1]!=0){
             tmp = tmp+std::conj(tmp);
           }
           tmp_result += tmp;
         }
         result.push_back( tmp_result );
       }
       return( result );
    }

    scitbx::af::shared< std::complex<FloatType> > coefs()
    {
      return( nl_.coefs() );
    }

    scitbx::af::shared< scitbx::af::tiny<int,2> > nl()
    {
      return( nl_.nl() );
    }


    private:
      int m_, n_max_;
      bool hex_;
      FloatType delta_, eps_;
      scitbx::af::shared< scitbx::vec2<FloatType> > xy_;
      scitbx::af::shared< scitbx::vec2<FloatType> > rt_;
      scitbx::af::shared< scitbx::vec2<int> > ij_;

      scitbx::af::shared< scitbx::af::tiny<int,2> > nl_indices_only_;

      scitbx::af::shared< std::vector< std::complex<FloatType> > > partial_data_;
      nl_complex_array<FloatType> nl_;
      log_factorial_generator<FloatType> lgf_;
      scitbx::af::shared< FloatType > model_;
  };







}}} // namespace scitbx::math::zernike

#endif // SCITBX_MATH_ZERNIKE_H
