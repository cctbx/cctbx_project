#ifndef SCITBX_ROTATE_H
#define SCITBX_ROTATE_H

#include<scitbx/constants.h>
#include<cmath>
#include<cfloat>
#include<scitbx/array_family/shared.h>
#include<scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include<scitbx/math/zernike.h>

namespace af=scitbx::af;

namespace scitbx { namespace math
{

  template <typename FloatType>
  class dmatrix
  {
    public:
      dmatrix() {}
      dmatrix( int max_L,
               FloatType beta
             ):
              max_L_(max_L), beta_(beta),
              lgf_(max_L), logBIGNUM_( std::log( DBL_MAX/1.0e15 ) )
      {
        for(int j=0;j<=max_L_;j++)
          dmatrix_.push_back( dj_table( j ));
      }

      FloatType djmn( int j, int m, int n) {
        return dmatrix_[j][j+m][j+n];
      }

      af::shared< af::shared<FloatType> > dj_table( int& j)
      {
        af::shared< af::shared< FloatType > > dj_t;

        for(int i=0;i<2*j+1;i++) {
          af::shared< FloatType > tmp_t( 2*j+1, 0.0);
          dj_t.push_back( tmp_t );
        }

        if( beta_ == 0 ) {
          for(int i=0;i<2*j+1;i++) {
          dj_t[i][i]=1.0;
          }
        return dj_t;
        } //end beta == 0

        FloatType djm_nm, djm_n, djm_np, csb, tmp;
        FloatType cos_half_beta( std::cos( beta_/2.0 )), sin_half_beta( std::sin(beta_/2.0));
        FloatType cosbeta( std::cos( beta_)), sinbeta( std::sin(beta_));
        FloatType powcosbeta(0.0), powsinbeta(0.0);
        for(int m=-j; m<=j; m++ ) {
          djm_np = 0.0;
          if (std::abs(cosbeta) > 0.001 || std::abs( std::log(cosbeta)*(j+m))<logBIGNUM_)
             powcosbeta = std::pow( cos_half_beta, j+m );
          if (std::abs(sinbeta) > 0.001 || std::abs( std::log(sinbeta)*(j+m))<logBIGNUM_)
             powsinbeta = std::pow( sin_half_beta, j-m );

          csb = powsinbeta*powcosbeta;
          tmp = std::exp( 0.5*(lgf_.log_fact(2*j)-lgf_.log_fact(j+m)-lgf_.log_fact(j-m) ) );

          djm_n = tmp * csb;  //djnm = (-1)^(m-n) djmn = dj-m-n
          dj_t[j+m][j+j] = djm_n;
          dj_t[j-m][j-j] = djm_n*pow_1(j-m);
          dj_t[j+j][j+m] = djm_n*pow_1(j-m);
          dj_t[j-j][j-m] = djm_n;

          for(int l=j-1; l>=scitbx::fn::absolute(m);l--) {
            djm_nm = -(std::sqrt( static_cast< FloatType > ((j+l+2)*(j-l-1)))*
                       djm_np + 2.0*(m-(l+1)*cosbeta)*djm_n/sinbeta)
                       /std::sqrt( static_cast< FloatType> ((j-l)*(j+l+1)));
            djm_np = djm_n;
            djm_n = djm_nm;
            dj_t[j+m][j+l] = djm_n;
            dj_t[j-m][j-l] = djm_n*pow_1(m-l);
            dj_t[j+l][j+m] = djm_n*pow_1(m-l);
            dj_t[j-l][j-m] = djm_n;
          }
        }
        //print_table( dj_t );
        return dj_t;
      }

      void print_table( af::shared< af::shared< FloatType> > table) {
        int j = table.size();
        for(int m=0;m<j;m++) {
          for(int m2=0;m2<j;m2++)
            std::cout<< table[m][m2]<<" ";
          std::cout<<std::endl;
          }
      }

      int pow_1( int n) {
        if( n == n/2*2 ) return 1;
        else return -1;
      }

    private:
      int max_L_;
      FloatType beta_;
      FloatType logBIGNUM_;
      scitbx::math::zernike::log_factorial_generator< FloatType > lgf_;
      af::shared< af::shared < af::shared< FloatType > > > dmatrix_;

  }; //end dmatrix

  template <typename FloatType>
  class correlation
  {
   public:
    correlation(
      scitbx::math::zernike::nlm_array<FloatType> const& F_nlm, //Source Coefs
      scitbx::math::zernike::nlm_array<FloatType> const& M_nlm, //Target Coefs
      int const& nmax,
      FloatType const& beta
      ): Fnlm_(F_nlm), Mnlm_(M_nlm), nmax_(nmax), beta_(beta), dM_(nmax, beta),
         complexI_(0.0,1.0), size_(2*nmax_+1), mm_grid_(size_, size_), mhm_grid_(size_, size_, size_),
         mm_(mm_grid_, 0.0), mhm_(mhm_grid_, 0.0)
    {
      calc_FMlmm();
    }

    void calc_FMlmm() {
      for(int l=0;l<=nmax_;l++) {
         af::shared< af::shared<std::complex< FloatType> > > mm;
         for(int m=-l;m<=l;m++) {
           af::shared< std::complex<FloatType> > m_array( 2*l+1, 0.0);
           mm.push_back( m_array );
           }
         FM_lmm_.push_back( mm );
      }

      std::complex<FloatType> s, t;
      for(int l=0;l<=nmax_;l++) {
        for(int m1=-l;m1<=l;m1++) {
          for(int m2=-l;m2<=l;m2++) {
            std::complex<FloatType> tmp_lmm(0,0);
            for(int n=l;n<=nmax_;n+=2) {
              s = Fnlm_.get_coef(n,l,m1);
              t = Mnlm_.get_coef(n,l,m2);
              tmp_lmm += std::conj(s) * t;
            }
            FM_lmm_[l][l+m1][l+m2] = tmp_lmm;
          }
        }
      }
      return;
    }

    void slow_calc_FMlmm() {
      for(int l=0;l<=nmax_;l++) {
         af::shared< af::shared<std::complex< FloatType> > > mm;
         for(int m=-l;m<=l;m++) {
           af::shared< std::complex<FloatType> > m_array( 2*l+1, 0.0);
           mm.push_back( m_array );
           }
         FM_lmm2_.push_back( mm );
      }
      std::complex<FloatType> s, t;
      for(int n=0;n<=nmax_;n++)
      for(int l=(n-n/2*2);l<=n;l+=2) {
        for(int m1=-l;m1<=l;m1++) {
          for(int m2=-l;m2<=l;m2++) {
              s = Fnlm_.get_coef(n,l,m1);
              t = Mnlm_.get_coef(n,l,m2);
              FM_lmm2_[l][l+m1][l+m2] += std::conj(s) * t;
          }
        }
      }
      return;
    }

    bool compare_FM() {
      slow_calc_FMlmm();
      for(int l=0;l<=nmax_;l++)
        for(int m1=-l;m1<=l;m1++)
          for(int m2=-l;m2<=l;m2++) {
            if(FM_lmm_[l][m1+l][m2+l] != FM_lmm2_[l][m1+l][m2+l] )
              return false;
            }

      return true;
    }

    af::versa<  std::complex<FloatType>, af::c_grid<3>  > mhm_coef(int border) {
      if( border == 0 ) return mhm();
      int new_size = 2*(nmax_+ border ) + 1;
      af::c_grid< 3 >mhm_grid(new_size, new_size, new_size);
      af::versa< std::complex< FloatType >, af::c_grid<3>  >  mhm( mhm_grid, 0.0 );
      for(int m1=0; m1< size_; m1++ )
        for(int h=0; h< size_; h++ )
          for(int m2=0; m2< size_; m2++ )
            mhm( m1 + border, h+border,  m2+border ) = mhm_( m1, h, m2 );
      return mhm;
    }

    af::versa<  std::complex<FloatType>, af::c_grid<3>  > mhm() {
      for(int i=0;i<mhm_.size();i++) mhm_[i] = 0;  //initialization
      set_beta(scitbx::constants::pi/2.0);
      for(int l=0;l<=nmax_;l++){
       for(int m1=-l; m1<=l; m1++ ) {
         for(int h=-l; h<=l; h++ ) {
           for(int m2=-l; m2<=l; m2++)
             mhm_( m1+nmax_, h+nmax_, m2+nmax_ ) += FM_lmm_[l][m1+l][m2+l]*dM_.djmn(l,m1,h)*dM_.djmn(l,h,m2);
         }
       }
      }
     return mhm_;
   }

    af::versa< std::complex< FloatType> , af::c_grid<2> > mm_coef( int border ) {
      if(border == 0) return mm();
      int new_size = 2*(nmax_+ border ) + 1;
      af::c_grid< 2 >mm_grid(new_size, new_size);
      af::versa< std::complex< FloatType >, af::c_grid<2>  >  mm( mm_grid, 0.0 );
      for(int m1=0; m1< size_; m1++ )
        for(int m2=0; m2< size_; m2++ )
          mm( m1 + border, m2+ border ) = mm_( m1, m2 );
      return mm;
    }

    af::versa<  std::complex<FloatType>, af::c_grid<2>  > mm() {
      for(int i=0;i<mm_.size();i++) mm_[i] = 0;
      for(int l=0;l<=nmax_;l++){
       for(int m1=-l; m1<=l; m1++ ) {
        for(int m2=-l; m2<=l; m2++)
          mm_( m1+nmax_, m2+nmax_) += FM_lmm_[l][m1+l][m2+l]*dM_.djmn(l,m1,m2);
       }
      }

      return mm_;
   }


    scitbx::math::zernike::nlm_array<FloatType>
    rotate_moving_obj( FloatType alpha, FloatType beta, FloatType gama) {
      scitbx::math::zernike::nlm_array<FloatType> result( nmax_ );
      dmatrix<FloatType> small_d( nmax_, beta );
      std::complex<FloatType>  Dlmn, exp_alpha, exp_gama, tmp_coef;
      af::shared<std::complex<FloatType> > a_array;
      af::shared<std::complex<FloatType> > g_array;

      for(int m=-nmax_;m<=nmax_;m++) {
        a_array.push_back( std::exp(complexI_ * (FloatType)m * alpha) );
        g_array.push_back( std::exp(complexI_ * (FloatType)m * gama) );
      }

      for(int n=0; n<=nmax_; n++) {
        for(int l=(n-n/2*2); l<=n; l+=2 ) {
          for(int m1=-l;m1<=l;m1++) {
            tmp_coef = 0.0;
            exp_alpha = a_array[m1+nmax_];
            for(int m2=-l;m2<=l;m2++) {
              exp_gama  = g_array[m2+nmax_];
              Dlmn = exp_alpha * small_d.djmn(l,m1,m2) * exp_gama;
              tmp_coef += Mnlm_.get_coef(n,l,m2) * Dlmn;
            }
            result.set_coef(n,l,m1, tmp_coef);
          }
        }
      }
      return result;
    }

    std::complex<FloatType> calc_correlation( FloatType alpha, FloatType beta, FloatType gama)
    {
      if(beta != beta_) {
        set_beta( beta );
      }

      int l, m1, m2;
      std::complex<FloatType>  Dlmn, exp_alpha, exp_gama;
      cc_ = 0.0;
      af::shared<std::complex<FloatType> > a_array;
      af::shared<std::complex<FloatType> > g_array;
      for(int m=-nmax_;m<=nmax_;m++) {
        a_array.push_back( std::exp(complexI_ * (FloatType)m * alpha) );
        g_array.push_back( std::exp(complexI_ * (FloatType)m * gama) );
      }

      for(l=0;l<=nmax_;l++) {
        for(m1=-l; m1<=l; m1++) {
          for(m2=-l; m2<=l; m2++) {
           exp_alpha = a_array[m1+nmax_];
           exp_gama  = g_array[m2+nmax_];
           Dlmn = exp_alpha * dM_.djmn(l,m1,m2) * exp_gama;
           cc_ += (FM_lmm_[l][m1+l][m2+l]*Dlmn );
          }
        }
      }
      return cc_;
    }

    void set_beta( FloatType beta) {
      beta_ = beta;
      dM_ = dmatrix<FloatType>( nmax_, beta_ );
      return;
    }



   private:
    int nmax_, size_;
    scitbx::math::zernike::nlm_array<FloatType> Fnlm_;
    scitbx::math::zernike::nlm_array<FloatType> Mnlm_;
    std::complex< FloatType> cc_;
    FloatType beta_;
    dmatrix<FloatType> dM_;
    af::shared< af::shared< af::shared< std::complex< FloatType > > > > FM_lmm_;
    af::shared< af::shared< af::shared< std::complex< FloatType > > > > FM_lmm2_;
    std::complex<FloatType> complexI_;
    af::c_grid<2> mm_grid_;
    af::c_grid<3> mhm_grid_;
    af::versa< std::complex< FloatType>, af::c_grid<2> > mm_;
    af::versa< std::complex< FloatType>, af::c_grid<3> > mhm_;

  }; //end correlation


}}//endof scitbx namespace



#endif
