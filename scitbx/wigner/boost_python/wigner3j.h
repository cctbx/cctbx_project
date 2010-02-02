#ifndef SCITBX_WIGNER3J_H
#define SCITBX_WIGNER3J_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <string>
#include <iomanip>

/*   this module is for INTEGER l,m calculation only, yet   */
namespace scitbx { namespace wigner {

   template<typename FloatType>
   class wigner3j_fast
   {
        public:
	  wigner3j_fast( int const& max ):
	  	max_(max)
          {
	    value_ = 0;
	    build_factorial();
	    //std::cout<<compute(0,0,0,0,0,0)<<std::endl; the output is 1, and this is verified
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
		value_=value_+pow(-1,tt)/(fact_array_[tt]*fact_array_[tt-t1]*fact_array_[tt-t2]*fact_array_[t3-tt]*fact_array_[t4-tt]*fact_array_[t5-tt]);
	    }
	    triangle_coef = fact_array_[j1+j2-j3]*fact_array_[j1-j2+j3]*fact_array_[-j1+j2+j3]/fact_array_[j1+j2+j3+1];
	    value_=value_*pow(-1,j1-j2-m3)*sqrt( triangle_coef * fact_array_[j1+m1]* fact_array_[j1-m1] * fact_array_[j2+m2] * fact_array_[j2-m2] * fact_array_[j3+m3] * fact_array_[j3-m3] );
	  return value_;	
	  }

	  FloatType get_value()
	  {
		return value_;
	  }

	  void build_factorial()
	  {
	    FloatType fac = 1;
	    fact_array_.push_back( fac );
	    for(int i=1; i<=max_; i++)
	    {
	      fac *= i;
	      fact_array_.push_back( fac );
	    }
	    return;
	  }

	  void update_factorial( int new_max )
	  {
	    FloatType fac = fact_array_[max_-1];
	    for(int i = max_; i<= new_max; i++)
	    {
		fac *= i;
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
		value_=value_+pow(-1,tt)/(fact(tt)*fact(tt-t1)*fact(tt-t2)*fact(t3-tt)*fact(t4-tt)*fact(t5-tt));
	    }
	    triangle_coef = fact(j1+j2-j3)*fact(j1-j2+j3)*fact(-j1+j2+j3)/fact(j1+j2+j3+1);
	    value_=value_*pow(-1,j1-j2-m3)*sqrt( triangle_coef * fact(j1+m1)* fact(j1-m1) * fact(j2+m2) * fact(j2-m2) * fact(j3+m3) * fact(j3-m3) );

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


   template<typename FloatType>
   class wigner3j_zero
   {
        public:
	  wigner3j_zero( int const& j1,
		    int const& j2,
		    int const& j3)
          {
	      compute(j1,j2,j3);
	  }
          
          void compute(int j1, int j2, int j3)
	  {
	    int t1, t2, t3, t4, t5, tmax, tmin;
	    double triangle_coef;
   	    if(( j1+ j2 < j3) ||( j3 < std::abs( j1-j2) )||j1<0||j2<0||j3<0)
	      {value_=0; return;}

	    t1 = j2 - j3;
	    t2 = j1 - j3;
	    t3 = j1 + j2 - j3;
	    t4 = j1 ;
	    t5 = j2 ;
	    
	    tmin = max(0, t1, t2 );	
	    tmax = min(t3,t4, t5 );
	    for(int tt=tmin; tt<=tmax; tt++)
	    {
		value_=value_+pow(-1,tt)/(fact(tt)*fact(tt-t1)*fact(tt-t2)*fact(t3-tt)*fact(t4-tt)*fact(t5-tt));
	    }
	    triangle_coef = fact(j1+j2-j3)*fact(j1-j2+j3)*fact(-j1+j2+j3)/fact(j1+j2+j3+1);
	    //value_=value_*pow(-1,j1-j2)*sqrt( triangle_coef )* fact(j1) * fact(j2) * fact(j3);
            value_=value_*pow(-1,j1-j2)*sqrt( triangle_coef * fact(j1)* fact(j1) * fact(j2) * fact(j2) * fact(j3) * fact(j3) );


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

   }; //end of wigner3j_zero class


} //namepsace wigner
} //namespace scitbx

#endif
