#ifndef SCITBX_WIGNER3J_H
#define SCITBX_WIGNER3J_H

#include <scitbx/constants.h>
#include <scitbx/array_family/shared.h>

namespace scitbx {
/* the routine will calculate the wigner 3j symbol value for integer parameters
 * See: http://mathworld.wolfram.com/Wigner3j-Symbol.html
 *
 *
 * Due to the range of double value, (j1+j2+j3+1)<150, ie, the max<150,
 * is required to get reliable values
 * ATTEN:  this module is for INTEGER l,m calculation only, yet   */

namespace wigner {

  template<typename FloatType>
  class wigner3j_fast
  {
    public:
    wigner3j_fast( int const& max ):
                 max_(max)
    {
      value_ = 0;
      build_factorial();
    }

    void check(int j1, int j2, int j3, int m1, int m2, int m3)
    {
      valid_=true;

      if (std::abs(m1) > j1)
        valid_= false;
      else if( std::abs(m2) > j2)
        valid_=false;
      else if( std::abs(m3) > j3)
        valid_=false;
      else if( m1 + m2 + m3 !=0)
        valid_=false;
      else if( j1+ j2 < j3)
        valid_=false;
      else if( j3 < std::abs( j1-j2) )
        valid_=false;
      else if( (m1==0) && (m2==0) && (m3==0) )
        {if( (j1+j2+j3) % 2 !=0)
          valid_=false; }
      return;
    }

    FloatType compute(int j1, int j2, int j3, int m1, int m2, int m3)
    {
      value_ = 0;
      check(j1,j2,j3,m1,m2,m3);
      if( ! valid_ ) return value_;

      t1 = j2 - m1 - j3;
      t2 = j1 + m2 - j3;
      t3 = j1 + j2 - j3;
      t4 = j1 - m1;
      t5 = j2 + m2;

      tmin = max(0, t1, t2 );
      tmax = min(t3,t4, t5 );
      for(int tt=tmin; tt<=tmax; tt++)
      {
          value_ += std::pow(-1.0,tt)/std::exp(
             fact_array_[tt]+fact_array_[tt-t1]
            +fact_array_[tt-t2]+fact_array_[t3-tt]
            +fact_array_[t4-tt]+fact_array_[t5-tt]);
      }
      triangle_coef = fact_array_[j1+j2-j3]+fact_array_[j1-j2+j3]
                      +fact_array_[-j1+j2+j3]-fact_array_[j1+j2+j3+1];
      value_ *= std::pow(-1.0,j1-j2-m3) *std::sqrt( std::exp( triangle_coef
        +fact_array_[j1+m1]+fact_array_[j1-m1]
        +fact_array_[j2+m2] + fact_array_[j2-m2]
        +fact_array_[j3+m3] + fact_array_[j3-m3]));
      return value_;
    }

     FloatType get_value()
     {
        return value_;
     }

     void build_factorial()
     {
       FloatType fac = 0.0;
       fact_array_.push_back( fac );
       for(int i=1; i<=max_; i++)
       {
         fac += std::log(static_cast<FloatType>(i));
         fact_array_.push_back( fac );
       }
       return;
     }

     void update_factorial( int new_max )
     {
       FloatType fac = fact_array_[max_-1];
       for(int i = max_; i<= new_max; i++)
       {
         fac += std::log(static_cast<FloatType>(i));
         fact_array_.push_back( fac );
       }
       max_ = new_max;
       return;
     }

   private:
     FloatType fact( int n)
     {
      if(n==0)
        return 1;
      else
        return n*fact(n-1);
     }

     int max(int a, int b, int c)
     {
       int m=a;
       if(b>m) m=b;
       if(c>m) m=c;
       return m;
     }

     int min(int a, int b, int c)
     {
       int m=a;
       if(b<m) m=b;
       if(c<m) m=c;
       return m;
     }

     FloatType value_, max_;
     scitbx::af::shared<FloatType> fact_array_;
     bool valid_;
     int t1, t2, t3, t4, t5, tmax, tmin;
     double triangle_coef;

   }; //end of wigner class


   template<typename FloatType>
   class wigner3j
   {
     public:
       wigner3j( int const& j1,
                 int const& j2,
                 int const& j3,
                 int const& m1,
                 int const& m2,
                 int const& m3 )
       {
         check(j1,j2,j3,m1,m2,m3);
         value_=0.0;
         if(valid_)
           compute(j1,j2,j3,m1,m2,m3);
       }

       void check(int j1, int j2, int j3, int m1, int m2, int m3)
       {
         if (std::abs(m1) > j1)
           valid_= false;
         else if( std::abs(m2) > j2)
           valid_=false;
         else if( std::abs(m3) > j3)
           valid_=false;
         else if( m1 + m2 + m3 !=0)
           valid_=false;
         else if( j1+ j2 < j3)
           valid_=false;
         else if( j3 < std::abs( j1-j2) )
           valid_=false;
         else
           valid_=true;
         return;
       }

       void compute(int j1, int j2, int j3, int m1, int m2, int m3)
       {
         int t1, t2, t3, t4, t5, tmax, tmin;
         double triangle_coef;
         t1 = j2 - m1 - j3;
         t2 = j1 + m2 - j3;
         t3 = j1 + j2 - j3;
         t4 = j1 - m1;
         t5 = j2 + m2;

         tmin = max(0, t1, t2 );
         tmax = min(t3,t4, t5 );
         for(int tt=tmin; tt<=tmax; tt++)
         {
           value_ += std::pow(-1.0,tt)/( fact(tt)*fact(tt-t1)*fact(tt-t2)
                                   *fact(t3-tt)*fact(t4-tt)*fact(t5-tt)
                                  );
         }
         triangle_coef = fact(j1+j2-j3)*fact(j1-j2+j3)*fact(-j1+j2+j3)
                        /fact(j1+j2+j3+1);

         value_ *= std::pow(-1.0,j1-j2-m3)*sqrt( triangle_coef * fact(j1+m1)
                                          * fact(j1-m1)*fact(j2+m2)*fact(j2-m2)
                                          * fact(j3+m3) * fact(j3-m3)
                                          );
         return;
       }

       FloatType get_value()
       {
         return value_;
       }

     private:
       FloatType fact( int n)
       {
         if(n==0)
           return 1;
         else
           return n*fact(n-1);
       }

       int max(int a, int b, int c)
       {
         int m=a;
         if(b>m) m=b;
         if(c>m) m=c;
         return m;
       }

       int min(int a, int b, int c)
       {
         int m=a;
         if(b<m) m=b;
         if(c<m) m=c;
         return m;
       }

       FloatType value_;
       bool valid_;

   }; //end of wigner class


} //namepsace wigner
} //namespace scitbx

#endif
