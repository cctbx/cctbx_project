
#include "abstract.h"
#include "expressions.h"
#include "shortcuts.h"

namespace cctbx { namespace sgtbx { namespace asu {

  typedef ivector3_t vvv;

  template< typename T> class wrapper : public abstract
  {
    public:
      T obj;

    void change_basis(const change_of_basis_op &o)
    {
      obj.change_basis(o);
    }

    void write(std::ostream &os) const
    {
      obj.print(os, true);
    }

    bool is_inside(const rvector3_t &p) const
    {
      return obj.is_inside(p);
    }
    
    size_type size() const 
    {
      return n_faces<T>::value;
    }

    void get_nth_plane(size_type i, cut &plane) const
    {
      cctbx::sgtbx::asu::get_nth_plane(obj, i, plane);  // obj must be and_expression
    }

    wrapper(const T &o) : obj(o) { if( ! is_proper_t<T>::value ) throw "improper type"; }

    abstract::ptr new_volume_only() const
    {
      typedef typename strip<T>::return_type return_type;
      return abstract::ptr( new wrapper< return_type >( strip<T>::execute(obj) ) );
    }

    abstract::ptr new_copy() const { return abstract::ptr( new wrapper<T>(*this) );}
  };

  template<typename TL, typename TR> 
    abstract::ptr abstract_asu(const and_expression<TL,TR> &expr)
  {
      return wrapper< and_expression<TL,TR> > ( expr ).new_copy();
  }


  //
  // the rest of the file was generated
  //

/////   Hall:  P 1
abstract::ptr  asu_001 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0
    & +z1
  );
}

/////   Hall:  -P 1
abstract::ptr  asu_002 ()
{
  return abstract_asu(
      x0(y0(z2) & y2(z2))
    & x2(y0(z2) & y2(z2))
    & y0
    & +y1
    & z0
    & +z1
  );
}

/////   Hall:  P 2y
abstract::ptr  asu_003 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(x2)
    & z2(x2)
  );
}

/////   Hall:  P 2yb
abstract::ptr  asu_004 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(x0(+y2) & x2(+y2))
    & z2(x0(+y2) & x2(+y2))
  );
}

/////   Hall:  C 2y
abstract::ptr  asu_005 ()
{
  return abstract_asu(
      x0(z2)
    & x2(z2)
    & y0
    & +y2
    & z0
    & +z1
  );
}

/////   Hall:  P -2y
abstract::ptr  asu_006 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & y2
    & z0
    & +z1
  );
}

/////   Hall:  P -2yc
abstract::ptr  asu_007 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  C -2y
abstract::ptr  asu_008 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & y4(+x2)
    & z0
    & +z1
  );
}

/////   Hall:  C -2yc
abstract::ptr  asu_009 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0(+z2)
    & y4(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  -P 2y
abstract::ptr  asu_010 ()
{
  return abstract_asu(
      x0(z2)
    & x2(z2)
    & y0
    & y2
    & z0
    & +z1
  );
}

/////   Hall:  -P 2yb
abstract::ptr  asu_011 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0(z0(x2) & z2(x2))
    & y4
    & z0
    & +z1
  );
}

/////   Hall:  -C 2y
abstract::ptr  asu_012 ()
{
  return abstract_asu(
      x0(z2)
    & x2(z2)
    & y0
    & y4(x4(z2))
    & z0
    & +z1
  );
}

/////   Hall:  -P 2yc
abstract::ptr  asu_013 ()
{
  return abstract_asu(
      x0(z0(y2) & z4)
    & x2(z0(y2) & z4)
    & y0
    & +y1
    & z0
    & +z2
  );
}

/////   Hall:  -P 2ybc
abstract::ptr  asu_014 ()
{
  return abstract_asu(
      x0(y0(z2))
    & +x1
    & y0(x2(z2))
    & y4(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  -C 2yc
abstract::ptr  asu_015 ()
{
  return abstract_asu(
      x0(z4)
    & x2(z4)
    & y0
    & +y2
    & z0(y4(x4))
    & z2(-y4(x4))
  );
}

/////   Hall:  P 2 2
abstract::ptr  asu_016 ()
{
  return abstract_asu(
      x0(z2)
    & x2(z2)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  );
}

/////   Hall:  P 2c 2
abstract::ptr  asu_017 ()
{
  return abstract_asu(
      x0(-z4 & z1*3/4)
    & x2(-z4 & z1*3/4)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  );
}

/////   Hall:  P 2 2ab
abstract::ptr  asu_018 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  );
}

/////   Hall:  P 2ac 2ab
abstract::ptr  asu_019 ()
{
  return abstract_asu(
      x0
    & +x2
    & y0(-z2)
    & y2(z2)
    & z0(+y2)
    & +z1
  );
}

/////   Hall:  C 2c 2
abstract::ptr  asu_020 ()
{
  return abstract_asu(
      x0(z4)
    & x2(-z4)
    & y0
    & y2(-z0)
    & z0
    & +z2
  );
}

/////   Hall:  C 2 2
abstract::ptr  asu_021 ()
{
  return abstract_asu(
      x0(z2)
    & x4(y4)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  );
}

/////   Hall:  F 2 2
abstract::ptr  asu_022 ()
{
  return abstract_asu(
      x0(z2)
    & x4(-z4 & z1*3/4)
    & y0(z2)
    & y4(-z4 & z1*3/4)
    & z0
    & +z1
  );
}

/////   Hall:  I 2 2
abstract::ptr  asu_023 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & y2(-z0)
    & z0
    & z2(-x0)
  );
}

/////   Hall:  I 2b 2c
abstract::ptr  asu_024 ()
{
  return abstract_asu(
      x0(y4)
    & x2(y4)
    & y0(z4)
    & y2(z4)
    & z0(x4)
    & z2(x4)
  );
}

/////   Hall:  P 2 -2
abstract::ptr  asu_025 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z1
  );
}

/////   Hall:  P 2c -2
abstract::ptr  asu_026 ()
{
  return abstract_asu(
      x0
    & x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  P 2 -2c
abstract::ptr  asu_027 ()
{
  return abstract_asu(
      x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  P 2 -2a
abstract::ptr  asu_028 ()
{
  return abstract_asu(
      x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z1
  );
}

/////   Hall:  P 2c -2ac
abstract::ptr  asu_029 ()
{
  return abstract_asu(
      x0(+z2)
    & x4(+z2)
    & y0
    & +y1
    & z0
    & +z1
  );
}

/////   Hall:  P 2 -2bc
abstract::ptr  asu_030 ()
{
  return abstract_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z2
  );
}

/////   Hall:  P 2ac -2
abstract::ptr  asu_031 ()
{
  return abstract_asu(
      x0
    & x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  P 2 -2ab
abstract::ptr  asu_032 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  );
}

/////   Hall:  P 2c -2n
abstract::ptr  asu_033 ()
{
  return abstract_asu(
      x0
    & +x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  P 2 -2n
abstract::ptr  asu_034 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  );
}

/////   Hall:  C 2 -2
abstract::ptr  asu_035 ()
{
  return abstract_asu(
      x0
    & x4(y4)
    & y0
    & y2
    & z0
    & +z1
  );
}

/////   Hall:  C 2c -2
abstract::ptr  asu_036 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & +y2
    & z0
    & +z2
  );
}

/////   Hall:  C 2 -2c
abstract::ptr  asu_037 ()
{
  return abstract_asu(
      x0(+z2)
    & x4(y4)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  A 2 -2
abstract::ptr  asu_038 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  );
}

/////   Hall:  A 2 -2b
abstract::ptr  asu_039 ()
{
  return abstract_asu(
      x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y4
    & z0
    & +z1
  );
}

/////   Hall:  A 2 -2a
abstract::ptr  asu_040 ()
{
  return abstract_asu(
      x0(+z2)
    & x4
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  A 2 -2ab
abstract::ptr  asu_041 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z2
  );
}

/////   Hall:  F 2 -2
abstract::ptr  asu_042 ()
{
  return abstract_asu(
      x0
    & x4(+z2)
    & y0
    & y4(+z2)
    & z0
    & +z1
  );
}

/////   Hall:  F 2 -2d
abstract::ptr  asu_043 ()
{
  return abstract_asu(
      x0
    & x4(-y0(+z2))
    & y0
    & +y4
    & z0
    & +z1
  );
}

/////   Hall:  I 2 -2
abstract::ptr  asu_044 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  );
}

/////   Hall:  I 2 -2c
abstract::ptr  asu_045 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z2
  );
}

/////   Hall:  I 2 -2a
abstract::ptr  asu_046 ()
{
  return abstract_asu(
      x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z2
  );
}

/////   Hall:  -P 2 2
abstract::ptr  asu_047 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & z2
  );
}

/////   Hall:  -P 2ab 2bc
abstract::ptr  asu_048 ()
{
  return abstract_asu(
      x0(-y0(z2))
    & x4(-z4 & z1*3/4)
    & ~y4(-z4 & z1*3/4)
    & y4(-z4 & z1*3/4)
    & z0
    & +z1
  );
}

/////   Hall:  -P 2 2c
abstract::ptr  asu_049 ()
{
  return abstract_asu(
      x0(z4)
    & x2(z4)
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  );
}

/////   Hall:  -P 2ab 2b
abstract::ptr  asu_050 ()
{
  return abstract_asu(
      x0(-y2)
    & x4(-y4 & y1*3/4)
    & y0
    & +y1
    & z0(-y4 & y1*3/4)
    & z2(-y4 & y1*3/4)
  );
}

/////   Hall:  -P 2a 2a
abstract::ptr  asu_051 ()
{
  return abstract_asu(
      x0(z2)
    & x4
    & y0
    & y2
    & z0
    & +z1
  );
}

/////   Hall:  -P 2a 2bc
abstract::ptr  asu_052 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0(-x4 & x1*3/4)
    & y4(z4)
    & z0(-x2)
    & z2(-x2)
  );
}

/////   Hall:  -P 2ac 2
abstract::ptr  asu_053 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & +y1
    & z0(y2)
    & z4(x4)
  );
}

/////   Hall:  -P 2a 2ac
abstract::ptr  asu_054 ()
{
  return abstract_asu(
      x0(-z4)
    & x2(z4)
    & y0(-x4)
    & y2(-x4)
    & z0
    & +z2
  );
}

/////   Hall:  -P 2 2ab
abstract::ptr  asu_055 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z2
  );
}

/////   Hall:  -P 2ab 2ac
abstract::ptr  asu_056 ()
{
  return abstract_asu(
      x0(y2(-z0))
    & x4(-y4 & y1*3/4)
    & y0
    & +y1
    & z0
    & +z2
  );
}

/////   Hall:  -P 2c 2b
abstract::ptr  asu_057 ()
{
  return abstract_asu(
      x0(-y2)
    & x2(-y2)
    & y0
    & +y1
    & z0(-y4 & y1*3/4)
    & z4
  );
}

/////   Hall:  -P 2 2n
abstract::ptr  asu_058 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z2
  );
}

/////   Hall:  -P 2ab 2a
abstract::ptr  asu_059 ()
{
  return abstract_asu(
      x0(-y0(z2))
    & x4
    & ~y4
    & y4
    & z0
    & +z1
  );
}

/////   Hall:  -P 2n 2ab
abstract::ptr  asu_060 ()
{
  return abstract_asu(
      x0(z4)
    & x2(-z4)
    & y0
    & y2(-x0(-z0))
    & z0
    & +z2
  );
}

/////   Hall:  -P 2ac 2ab
abstract::ptr  asu_061 ()
{
  return abstract_asu(
      x0
    & x2(-y0(-z0))
    & y0
    & +y2
    & z0
    & +z2
  );
}

/////   Hall:  -P 2ac 2n
abstract::ptr  asu_062 ()
{
  return abstract_asu(
      x0
    & x2(-y0(-z0))
    & y0(+z2)
    & y4
    & z0
    & +z1
  );
}

/////   Hall:  -C 2c 2
abstract::ptr  asu_063 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & +y2
    & z0(y4(x4))
    & z4
  );
}

/////   Hall:  -C 2ac 2
abstract::ptr  asu_064 ()
{
  return abstract_asu(
      x0
    & x4(z4)
    & y0
    & +y2
    & z0(y4)
    & z2(+y4)
  );
}

/////   Hall:  -C 2 2
abstract::ptr  asu_065 ()
{
  return abstract_asu(
      x0
    & x4(y4)
    & y0
    & y2
    & z0
    & z2
  );
}

/////   Hall:  -C 2 2c
abstract::ptr  asu_066 ()
{
  return abstract_asu(
      x0(z4)
    & x4(y4)
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  );
}

/////   Hall:  -C 2a 2
abstract::ptr  asu_067 ()
{
  return abstract_asu(
      x0
    & x2
    & y0(x4)
    & y4
    & z0(x4)
    & z2(x4)
  );
}

/////   Hall:  -C 2a 2ac
abstract::ptr  asu_068 ()
{
  return abstract_asu(
      x0(z4)
    & x2(z4)
    & y0(x4)
    & y4(z4)
    & z0(+x2 & y4(x4))
    & +z2
  );
}

/////   Hall:  -F 2 2
abstract::ptr  asu_069 ()
{
  return abstract_asu(
      x0
    & x4(z4)
    & y0
    & y4(z4)
    & z0
    & z2
  );
}

/////   Hall:  -F 2uv 2vw
abstract::ptr  asu_070 ()
{
  return abstract_asu(
      x0(-y0(z2))
    & x8(-z8 & z1*5/8)
    & ~y8(-z1*3/8 & z1*7/8)
    & y8(-z8 & z1*5/8)
    & z0
    & +z1
  );
}

/////   Hall:  -I 2 2
abstract::ptr  asu_071 ()
{
  return abstract_asu(
      x0
    & x4(y4(z4))
    & y0
    & y2
    & z0
    & z2
  );
}

/////   Hall:  -I 2 2c
abstract::ptr  asu_072 ()
{
  return abstract_asu(
      x0(z4)
    & x4(y4(z4))
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  );
}

/////   Hall:  -I 2b 2c
abstract::ptr  asu_073 ()
{
  return abstract_asu(
      x0(y4)
    & x4(z4(y4))
    & y0(z4)
    & y2(-z4)
    & z0
    & +z2
  );
}

/////   Hall:  -I 2b 2
abstract::ptr  asu_074 ()
{
  return abstract_asu(
      x0
    & x4(-z4 & z1*3/4)
    & y0(z2)
    & y4
    & z0
    & +z1
  );
}

/////   Hall:  P 4
abstract::ptr  asu_075 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z1
  );
}

/////   Hall:  P 4w
abstract::ptr  asu_076 ()
{
  return abstract_asu(
      x0(+z4)
    & x2(+z4)
    & y0(+z1*3/4)
    & y2(+z1*3/4)
    & z0
    & +z1
  );
}

/////   Hall:  P 4c
abstract::ptr  asu_077 ()
{
  return abstract_asu(
      x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

/////   Hall:   P 4cw
abstract::ptr  asu_078 ()
{
  return abstract_asu(
      x0(+-z1*3/4)
    & x2(+-z1*3/4)
    & y0(+-z4)
    & y2(+-z4)
    & z1
    & +z0
  );
}

/////   Hall:  I 4
abstract::ptr  asu_079 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  );
}

/////   Hall:  I 4bw
abstract::ptr  asu_080 ()
{
  return abstract_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z4
  );
}

/////   Hall:  P -4
abstract::ptr  asu_081 ()
{
  return abstract_asu(
      x0(-y0(z2))
    & x2
    & y0
    & y2(-x2(z2))
    & z0
    & +z1
  );
}

/////   Hall:  I -4
abstract::ptr  asu_082 ()
{
  return abstract_asu(
      x0(z0(-y0))
    & x2(-y0(z4))
    & y0
    & y2(-x0(z4))
    & z0
    & z2(-y0)
  );
}

/////   Hall:  -P 4
abstract::ptr  asu_083 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z2
  );
}

/////   Hall:  -P 4c
abstract::ptr  asu_084 ()
{
  return abstract_asu(
      x0(-y0(z4))
    & x2
    & y0
    & y2(-x2(z4))
    & z0
    & z2
  );
}

/////   Hall:  -P 4a
abstract::ptr  asu_085 ()
{
  return abstract_asu(
      ~x4(-~y4)
    & x4(z0(-~y4) & z2(-~y4))
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0))
    & z2(-y0(-x0))
  );
}

/////   Hall:  -P 4bc
abstract::ptr  asu_086 ()
{
  return abstract_asu(
      ~x4(-~y4(z4))
    & x4(z0(-~y4) & z2(+-~y4))
    & ~y4
    & y4(-x4(z4))
    & z0(-y0(-x0))
    & z2(-y0(-x0))
  );
}

/////   Hall:  -I 4
abstract::ptr  asu_087 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(y4(x4) & x2(-y0))
  );
}

/////   Hall:  -I 4ad
abstract::ptr  asu_088 ()
{
  return abstract_asu(
      x0
    & x4
    & y0(-x0(z2) | -x4(+z4))
    & y4(-x0(-z8 & z1*5/8))
    & z0
    & +z1
  );
}

/////   Hall:  P 4 2
abstract::ptr  asu_089 ()
{
  return abstract_asu(
      x0(p0)
    & x2
    & y0
    & y2(-x2)
    & z0(p0)
    & z2(p0)
  );
}

/////   Hall:  P 4ab 2ab
abstract::ptr  asu_090 ()
{
  return abstract_asu(
      x0
    & x2(-y0)
    & y0
    & y2(-x0)
    & z0(p0)
    & z2(p0)
  );
}

/////   Hall:  P 4w 2c
abstract::ptr  asu_091 ()
{
  return abstract_asu(
      x0(z8(-y0))
    & +x1
    & y0
    & +y1
    & z0(x2)
    & z8(m1)
  );
}

/////   Hall:  P 4abw 2nw
abstract::ptr  asu_092 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z8(-y2)
  );
}

/////   Hall:  P 4c 2
abstract::ptr  asu_093 ()
{
  return abstract_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(y2)
    & z4(-p0 & m1)
  );
}

/////   Hall:  P 4n 2n
abstract::ptr  asu_094 ()
{
  return abstract_asu(
      x0(-y0)
    & x2(z2(-y2))
    & y0(z2(-x0))
    & +y2
    & z0(p0)
    & z2(p0)
  );
}

/////   Hall:   P 4cw 2c
abstract::ptr  asu_095 ()
{
  return abstract_asu(
      x1(z8(-y0))
    & +x0
    & y0
    & +y1
    & z0(-x2)
    & z8(p0)
  );
}

/////   Hall:  P 4nw 2abw
abstract::ptr  asu_096 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z8(-x2)
  );
}

/////   Hall:  I 4 2
abstract::ptr  asu_097 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0(p0)
    & z4(m2)
  );
}

/////   Hall:  I 4bw 2bw
abstract::ptr  asu_098 ()
{
  return abstract_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(m1 & -p0)
    & z8(-y4 & y1*3/4)
  );
}

/////   Hall:  P 4 -2
abstract::ptr  asu_099 ()
{
  return abstract_asu(
      x0
    & y2
    & z0
    & +z1
    & -p0
  );
}

/////   Hall:  P 4 -2ab
abstract::ptr  asu_100 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z1
    & m2
  );
}

/////   Hall:  P 4c -2c
abstract::ptr  asu_101 ()
{
  return abstract_asu(
      x0(+z2)
    & y2(+z2)
    & z0
    & +z1
    & -p0
  );
}

/////   Hall:  P 4n -2n
abstract::ptr  asu_102 ()
{
  return abstract_asu(
      x0(+z2)
    & y2(+z2)
    & z0
    & +z1
    & -p0
  );
}

/////   Hall:  P 4 -2c
abstract::ptr  asu_103 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  );
}

/////   Hall:  P 4 -2n
abstract::ptr  asu_104 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  );
}

/////   Hall:  P 4c -2
abstract::ptr  asu_105 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  );
}

/////   Hall:  P 4c -2ab
abstract::ptr  asu_106 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & +y2
    & z0
    & +z2
  );
}

/////   Hall:  I 4 -2
abstract::ptr  asu_107 ()
{
  return abstract_asu(
      x0
    & y2
    & z0
    & +z2
    & -p0
  );
}

/////   Hall:  I 4 -2c
abstract::ptr  asu_108 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & m2
  );
}

/////   Hall:  I 4bw -2
abstract::ptr  asu_109 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z4
  );
}

/////   Hall:  I 4bw -2c
abstract::ptr  asu_110 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & +y2
    & z0
    & +z4
  );
}

/////   Hall:  P -4 2
abstract::ptr  asu_111 ()
{
  return abstract_asu(
      x0(z2)
    & y2(z2)
    & z0
    & +z1
    & -p0
  );
}

/////   Hall:  P -4 2c
abstract::ptr  asu_112 ()
{
  return abstract_asu(
      x0(z4 & z0(-y0))
    & x2(z4)
    & y0(z4)
    & y2(z4 & z0(-x2))
    & z0
    & +z2
  );
}

/////   Hall:  P -4 2ab
abstract::ptr  asu_113 ()
{
  return abstract_asu(
      x0(-y0(z2))
    & y0
    & z0
    & +z1
    & m2
  );
}

/////   Hall:  P -4 2n
abstract::ptr  asu_114 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2(-z0))
    & z0
    & +z2
  );
}

/////   Hall:  P -4 -2
abstract::ptr  asu_115 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0(p0)
    & z2(p0)
  );
}

/////   Hall:  P -4 -2c
abstract::ptr  asu_116 ()
{
  return abstract_asu(
      x0(y2)
    & x2(y2 & z0(-y2))
    & y0
    & +y1
    & z0(y2 & y0(-x0))
    & z4(m1 & -p0)
  );
}

/////   Hall:  P -4 -2ab
abstract::ptr  asu_117 ()
{
  return abstract_asu(
      x0(z0(-y0) & z2(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0(m2)
    & z2(m2)
  );
}

/////   Hall:  P -4 -2n
abstract::ptr  asu_118 ()
{
  return abstract_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(y2(-x2) & x0(-y0))
    & z4(~p2 & -m2)
  );
}

/////   Hall:  I -4 -2
abstract::ptr  asu_119 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0(p0)
    & z4(m2)
  );
}

/////   Hall:  I -4 -2c
abstract::ptr  asu_120 ()
{
  return abstract_asu(
      x0(z0(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0(m2)
    & z4(p0)
  );
}

/////   Hall:  I -4 2
abstract::ptr  asu_121 ()
{
  return abstract_asu(
      x0
    & y2(-x0(z4))
    & z0
    & z2(-x0)
    & -p0
  );
}

/////   Hall:  I -4 2bw
abstract::ptr  asu_122 ()
{
  return abstract_asu(
      x0(y2 & z0(-y0))
    & x2(y2)
    & y0
    & +y1
    & z0(y2(-x2))
    & z8(-y4 & y1*3/4)
  );
}

/////   Hall:  -P 4 2
abstract::ptr  asu_123 ()
{
  return abstract_asu(
      x0
    & y2
    & z0
    & z2
    & -p0
  );
}

/////   Hall:  -P 4 2c
abstract::ptr  asu_124 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(p0)
  );
}

/////   Hall:  -P 4a 2b
abstract::ptr  asu_125 ()
{
  return abstract_asu(
      ~x4(-~y4)
    & ~y4
    & z0(p0)
    & z2(p0)
    & -m0
  );
}

/////   Hall:  -P 4a 2bc
abstract::ptr  asu_126 ()
{
  return abstract_asu(
      ~x4(-~y4)
    & x4
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0) & x4(-~y4))
    & z4(p0)
  );
}

/////   Hall:  -P 4 2ab
abstract::ptr  asu_127 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & z2
    & m2
  );
}

/////   Hall:  -P 4 2n
abstract::ptr  asu_128 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(m2)
  );
}

/////   Hall:  -P 4a 2a
abstract::ptr  asu_129 ()
{
  return abstract_asu(
      ~x4
    & y4
    & z0(-m0)
    & z2(-m0)
    & -p0
  );
}

/////   Hall:  -P 4a 2ac
abstract::ptr  asu_130 ()
{
  return abstract_asu(
      ~x4(-~y4)
    & x4(z0(-~y4))
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0))
    & z4(-m0)
  );
}

/////   Hall:  -P 4c 2
abstract::ptr  asu_131 ()
{
  return abstract_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & z4(p0)
  );
}

/////   Hall:  -P 4c 2c
abstract::ptr  asu_132 ()
{
  return abstract_asu(
      x0(z4)
    & y2(z4)
    & z0
    & z2
    & -p0
  );
}

/////   Hall:  -P 4ac 2b
abstract::ptr  asu_133 ()
{
  return abstract_asu(
      ~x4
    & x4(-z0 | -~y4)
    & ~y4
    & +y4
    & z0(-y0(-x0))
    & z4(p0)
  );
}

/////   Hall:  -P 4ac 2bc
abstract::ptr  asu_134 ()
{
  return abstract_asu(
      ~x4(z4)
    & ~y4(z4)
    & z0(p0)
    & z2(p0)
    & -m0
  );
}

/////   Hall:  -P 4c 2ab
abstract::ptr  asu_135 ()
{
  return abstract_asu(
      x0(z4(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z4(m2)
  );
}

/////   Hall:  -P 4n 2n
abstract::ptr  asu_136 ()
{
  return abstract_asu(
      x0(z4)
    & y2(z4(-x0))
    & z0
    & z2
    & -p0
  );
}

/////   Hall:  -P 4ac 2a
abstract::ptr  asu_137 ()
{
  return abstract_asu(
      ~x4
    & x4
    & ~y4
    & y4
    & z0(-y0(-x0))
    & z4(-m0)
  );
}

/////   Hall:  -P 4ac 2ac
abstract::ptr  asu_138 ()
{
  return abstract_asu(
      ~x4
    & y4(-~x4(z4))
    & z0(-m0)
    & z2(-m0 & ~x4(-y4))
    & -p0
  );
}

/////   Hall:  -I 4 2
abstract::ptr  asu_139 ()
{
  return abstract_asu(
      x0
    & y2
    & z0
    & z4(m2)
    & -p0
  );
}

/////   Hall:  -I 4 2c
abstract::ptr  asu_140 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & z4(p0)
    & m2
  );
}

/////   Hall:  -I 4bd 2
abstract::ptr  asu_141 ()
{
  return abstract_asu(
      x0
    & x2
    & ~y4
    & y4
    & z0(-y0)
    & z8(-p4)
  );
}

/////   Hall:  -I 4bd 2c
abstract::ptr  asu_142 ()
{
  return abstract_asu(
      x0(z8(-~y4) & z0(-y0))
    & x2(-~y4)
    & ~y4
    & +y4
    & z0(x4)
    & z8(m1/4)
  );
}

/////   Hall:  P 3
abstract::ptr  asu_143 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z1
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 31
abstract::ptr  asu_144 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0
    & +z3
  );
}

/////   Hall:   P 32
abstract::ptr  asu_145 ()
{
  return abstract_asu(
      y0
    & +y1
    & x0
    & +x1
    & z0
    & +z3
  );
}

/////   Hall:  R 3
abstract::ptr  asu_146 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z3
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -P 3
abstract::ptr  asu_147 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z2(p0(-y0))
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -R 3
abstract::ptr  asu_148 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z6(-h0(x3) | -k0(-y0 | -m1))
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 3 2
abstract::ptr  asu_149 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(-h0 | -k0)
    & z2(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 3 2"
abstract::ptr  asu_150 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z2(p0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 31 2 (x,y,z+1/3)
abstract::ptr  asu_151 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(-h0 | -h1)
    & z6(-k0 | -k1)
  );
}

/////   Hall:  P 31 2"
abstract::ptr  asu_152 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(-p0)
    & z6(-p0)
  );
}

/////   Hall:  P 32 2 (x,y,z+1/6)
abstract::ptr  asu_153 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(-h0 | -h1)
    & z6(x0(-y0) & m1)
  );
}

/////   Hall:   P 32 2"
abstract::ptr  asu_154 ()
{
  return abstract_asu(
      y0
    & +y1
    & x0
    & +x1
    & z0(p0)
    & z6(p0)
  );
}

/////   Hall:  R 3 2"
abstract::ptr  asu_155 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z6(x3 & ~p3)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 3 -2"
abstract::ptr  asu_156 ()
{
  return abstract_asu(
      z0
    & +z1
    & h0
    & m1
    & k0
  );
}

/////   Hall:  P 3 -2
abstract::ptr  asu_157 ()
{
  return abstract_asu(
      y0
    & z0
    & +z1
    & k1
    & m1(y3)
    & p0
  );
}

/////   Hall:  P 3 -2"c
abstract::ptr  asu_158 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 3 -2c
abstract::ptr  asu_159 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  R 3 -2"
abstract::ptr  asu_160 ()
{
  return abstract_asu(
      z0
    & +z3
    & h0
    & m1
    & k0
  );
}

/////   Hall:  R 3 -2"c
abstract::ptr  asu_161 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z6
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -P 3 2
abstract::ptr  asu_162 ()
{
  return abstract_asu(
      y0
    & z0(-h0)
    & z2(-h0)
    & k1
    & m1(y3)
    & p0
  );
}

/////   Hall:  -P 3 2c
abstract::ptr  asu_163 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -P 3 2"
abstract::ptr  asu_164 ()
{
  return abstract_asu(
      y0(z2)
    & z0
    & +z1
    & k1
    & -h0
  );
}

/////   Hall:  -P 3 2"c
abstract::ptr  asu_165 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4(p0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -R 3 2"
abstract::ptr  asu_166 ()
{
  return abstract_asu(
      z0(p0)
    & z6(x3)
    & h0
    & m1
    & k0
  );
}

/////   Hall:  -R 3 2"c
abstract::ptr  asu_167 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z12(y3 & p3)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 6
abstract::ptr  asu_168 ()
{
  return abstract_asu(
      y0
    & z0
    & +z1
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

/////   Hall:  P 61
abstract::ptr  asu_169 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0
    & +z6
  );
}

/////   Hall:   P 65
abstract::ptr  asu_170 ()
{
  return abstract_asu(
      y0
    & +y1
    & x0
    & +x1
    & z0
    & +z6
  );
}

/////   Hall:  P 62
abstract::ptr  asu_171 ()
{
  return abstract_asu(
      x1(y2)
    & y0(x2)
    & z0
    & +z3
    & p0(y2)
  );
}

/////   Hall:   P 64
abstract::ptr  asu_172 ()
{
  return abstract_asu(
      y0(-x2)
    & x1(-y2)
    & z0
    & +z3
    & p0(-x2)
  );
}

/////   Hall:  P 6c
abstract::ptr  asu_173 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P -6
abstract::ptr  asu_174 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0
    & z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -P 6
abstract::ptr  asu_175 ()
{
  return abstract_asu(
      y0
    & z0
    & z2
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

/////   Hall:  -P 6c
abstract::ptr  asu_176 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 6 2
abstract::ptr  asu_177 ()
{
  return abstract_asu(
      y0
    & z0(-h0)
    & z2(-h0)
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

/////   Hall:  P 61 2 (x,y,z+5/12)
abstract::ptr  asu_178 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z12(-h0 | -h1)
  );
}

/////   Hall:  P 65 2 (x,y,z+1/12)
abstract::ptr  asu_179 ()
{
  return abstract_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z12(m1 & x0(-y0))
  );
}

/////   Hall:  P 62 2 (x,y,z+1/3)
abstract::ptr  asu_180 ()
{
  return abstract_asu(
      x1(y2)
    & y0(x2)
    & z0(k1)
    & z6(-h0)
    & p0(y2)
  );
}

/////   Hall:   P 64 2 (x,y,z+1/6)
abstract::ptr  asu_181 ()
{
  return abstract_asu(
      y0(p2)
    & p0(-y2)
    & z6(-m1)
    & z0(k1)
    & x1(p2)
  );
}

/////   Hall:  P 6c 2c
abstract::ptr  asu_182 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z4(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P 6 -2
abstract::ptr  asu_183 ()
{
  return abstract_asu(
      y0
    & z0
    & +z1
    & k1
    & -h0
  );
}

/////   Hall:  P 6 -2c
abstract::ptr  asu_184 ()
{
  return abstract_asu(
      y0
    & z0
    & +z2
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

/////   Hall:  P 6c -2
abstract::ptr  asu_185 ()
{
  return abstract_asu(
      y0
    & z0
    & +z2
    & k1
    & m1(y3)
    & p0
  );
}

/////   Hall:  P 6c -2c
abstract::ptr  asu_186 ()
{
  return abstract_asu(
      y0(+z2)
    & z0
    & +z1
    & k1
    & -h0
  );
}

/////   Hall:  P -6 2
abstract::ptr  asu_187 ()
{
  return abstract_asu(
      z0
    & z2
    & h0
    & m1
    & k0
  );
}

/////   Hall:  P -6c 2
abstract::ptr  asu_188 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(-h0 | -k0)
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  P -6 -2
abstract::ptr  asu_189 ()
{
  return abstract_asu(
      y0
    & z0
    & z2
    & k1
    & m1(y3)
    & p0
  );
}

/////   Hall:  P -6c -2c
abstract::ptr  asu_190 ()
{
  return abstract_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

/////   Hall:  -P 6 2
abstract::ptr  asu_191 ()
{
  return abstract_asu(
      y0
    & z0
    & z2
    & k1
    & -h0
  );
}

/////   Hall:  -P 6 2c
abstract::ptr  asu_192 ()
{
  return abstract_asu(
      y0
    & z0
    & z4(-h0)
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

/////   Hall:  -P 6c 2
abstract::ptr  asu_193 ()
{
  return abstract_asu(
      y0
    & z0(-h0)
    & z4
    & k1
    & m1(y3)
    & p0
  );
}

/////   Hall:  -P 6c 2c
abstract::ptr  asu_194 ()
{
  return abstract_asu(
      z0(p0)
    & z4
    & h0
    & m1
    & k0
  );
}

/////   Hall:  P 2 2 3
abstract::ptr  asu_195 ()
{
  return abstract_asu(
      z0(y2 & x2)
    & m1(-y2)
    & cut(vvv(0,1,-1),0)(cut(vvv(-1,0,1),0))
    & cut(vvv(1,0,-1),0)
  );
}

/////   Hall:  F 2 2 3
abstract::ptr  asu_196 ()
{
  return abstract_asu(
      p0(m2)
    & cut(vvv(-1,0,-1),r1/2)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,1),r1/2)(cut(vvv(0,-1,-1),0))
    & cut(vvv(0,1,1),0)
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  I 2 2 3
abstract::ptr  asu_197 ()
{
  return abstract_asu(
      z0(x2)
    & p0(cut(vvv(0,-1,1),0))
    & +m1
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  P 2ac 2ab 3
abstract::ptr  asu_198 ()
{
  return abstract_asu(
      x0(-y0)
    & x2
    & y2(+z0 & x2(+z2))
    & cut(vvv(-1,0,1),r1/2)(m2)
    & cut(vvv(1,0,-1),0)(p0)
    & cut(vvv(0,1,1),0)
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  I 2b 2c 3
abstract::ptr  asu_199 ()
{
  return abstract_asu(
      x2(-y4)
    & y2(-z4)
    & z0(x4)
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0)(+x2))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -P 2 2 3
abstract::ptr  asu_200 ()
{
  return abstract_asu(
      x2
    & y2
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -P 2ab 2bc 3
abstract::ptr  asu_201 ()
{
  return abstract_asu(
      ~z4(x4)
    & p0(cut(vvv(0,-1,1),0)(-x0))
    & m2(cut(vvv(0,-1,1),0)(x2))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -F 2 2 3
abstract::ptr  asu_202 ()
{
  return abstract_asu(
      z0
    & p0(x4)
    & cut(vvv(-1,0,-1),r1/2)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -F 2uv 2vw 3
abstract::ptr  asu_203 ()
{
  return abstract_asu(
      p0(cut(vvv(0,-1,1),0)(-x0))
    & (m1/4)(cut(vvv(0,-1,1),0) | cut(vvv(0,-1,-1),-r1/4)(~z4))
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),r1/4)
  );
}

/////   Hall:  -I 2 2 3
abstract::ptr  asu_204 ()
{
  return abstract_asu(
      x2
    & z0
    & p0(cut(vvv(0,-1,1),0)(x4))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -P 2ac 2ab 3
abstract::ptr  asu_205 ()
{
  return abstract_asu(
      x2(-z0(cut(vvv(0,-1,1),0)))
    & y2(cut(vvv(0,-1,1),0))
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -I 2b 2c 3
abstract::ptr  asu_206 ()
{
  return abstract_asu(
      z0(x4)
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,-1),r1/2)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,-1,-1),r1/2)(cut(vvv(-1,0,1),0))
  );
}

/////   Hall:  P 4 2 3
abstract::ptr  asu_207 ()
{
  return abstract_asu(
      z0(x2)
    & p0
    & m1(-p0)
    & cut(vvv(0,1,-1),0)(x2)
  );
}

/////   Hall:  P 4n 2 3
abstract::ptr  asu_208 ()
{
  return abstract_asu(
      cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(1,0,1),0)(cut(vvv(0,-1,-1),0))
    & cut(vvv(-1,0,1),r1/2)(y4)
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)
    & cut(vvv(0,-1,1),r1/2)(-x4)
    & cut(vvv(0,-1,-1),r1/2)(-x4)
  );
}

/////   Hall:  F 4 2 3
abstract::ptr  asu_209 ()
{
  return abstract_asu(
      p0(z0)
    & m2(z0)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)(cut(vvv(0,-1,1),0))
  );
}

/////   Hall:  F 4d 2 3
abstract::ptr  asu_210 ()
{
  return abstract_asu(
      y8(cut(vvv(1,0,1),-r1/4))
    & z8(m1/4)
    & p0(cut(vvv(-1,0,1),0))
    & m2(cut(vvv(1,0,1),-r1/2))
    & cut(vvv(0,1,1),0)(z0)
    & cut(vvv(1,0,-1),0)
    & cut(vvv(-1,0,-1),r1/2)
  );
}

/////   Hall:  I 4 2 3
abstract::ptr  asu_211 ()
{
  return abstract_asu(
      z0(p0)
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,-1,-1),r1/2)(-x4)
  );
}

/////   Hall:  P 4acd 2ab 3
abstract::ptr  asu_212 ()
{
  return abstract_asu(
      cut(vvv(-1,0,1),r1/2)
    & cut(vvv(0,1,1),0)(cut(vvv(1,0,-1),-r1/2))
    & cut(vvv(0,-1,-1),r1/2)(cut(vvv(-2,1,1),0))
    & cut(vvv(2,-1,-1),0)(x8)
    & cut(vvv(-1,2,-1),0)(y8)
    & cut(vvv(-2,1,-1),r1/2)(-x1*3/8)
  );
}

/////   Hall:   P 4bd 2ab 3
abstract::ptr  asu_213 ()
{
  return abstract_asu(
      cut(vvv(0,1,-1),0)
    & -p0(cut(vvv(0,-1,1),0))
    & ~p2(cut(vvv(-1,1,-2),0))
    & cut(vvv(1,-1,2),0)(z8)
    & cut(vvv(-2,-1,-1),r1*3/2)(-x1*3/8)
    & cut(vvv(-1,-1,-2),r1*3/2)(-z1*3/8)
  );
}

/////   Hall:  I 4bd 2c 3
abstract::ptr  asu_214 ()
{
  return abstract_asu(
      x8(cut(vvv(0,-1,-1),r1/4))
    & y8(cut(vvv(-1,0,-1),r1/4))
    & ~y8(cut(vvv(-1,0,1),-r1/4))
    & cut(vvv(-1,0,1),0)(cut(vvv(0,1,-1),0))
    & cut(vvv(0,-1,1),0)
    & cut(vvv(0,1,-1),r1/4)(-y0)
    & cut(vvv(1,-1,1),r1/8)(~p4)
  );
}

/////   Hall:  P -4 2 3
abstract::ptr  asu_215 ()
{
  return abstract_asu(
      z0(x2)
    & p0
    & m1
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  F -4 2 3
abstract::ptr  asu_216 ()
{
  return abstract_asu(
      p0
    & m2
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)
  );
}

/////   Hall:  I -4 2 3
abstract::ptr  asu_217 ()
{
  return abstract_asu(
      x2(-z0(y4))
    & z0
    & p0
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  P -4n 2 3
abstract::ptr  asu_218 ()
{
  return abstract_asu(
      x2(-z0(y4))
    & y2(-z0(x4))
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  F -4a 2 3
abstract::ptr  asu_219 ()
{
  return abstract_asu(
      p0
    & m2(-p0(z0))
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)(-p0 | cut(vvv(0,-1,1),0)(x4))
  );
}

/////   Hall:  I -4bd 2c 3
abstract::ptr  asu_220 ()
{
  return abstract_asu(
      -x4(-z0(-y1*3/8))
    & x2
    & -y4(-x2(-z8))
    & y2(-z4)
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -P 4 2 3
abstract::ptr  asu_221 ()
{
  return abstract_asu(
      x2
    & z0
    & p0
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -P 4a 2bc 3
abstract::ptr  asu_222 ()
{
  return abstract_asu(
      (x1*3/4)(z4(y2) | cut(vvv(0,-1,1),0))
    & -z4
    & p0(cut(vvv(0,-1,1),0)(z2))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -P 4n 2 3
abstract::ptr  asu_223 ()
{
  return abstract_asu(
      z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,-1,-1),r1/2)(x4)
  );
}

/////   Hall:  -P 4bc 2bc 3
abstract::ptr  asu_224 ()
{
  return abstract_asu(
      p0
    & cut(vvv(-1,0,-1),r1)(y2)
    & cut(vvv(-1,0,1),r1/2)(y2)
    & cut(vvv(0,1,1),-r1/2)
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -F 4 2 3
abstract::ptr  asu_225 ()
{
  return abstract_asu(
      z0
    & p0
    & m2
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -F 4a 2 3
abstract::ptr  asu_226 ()
{
  return abstract_asu(
      z0
    & p0
    & m2(-p0)
    & cut(vvv(0,1,-1),0)(x4)
  );
}

/////   Hall:  -F 4vw 2vw 3
abstract::ptr  asu_227 ()
{
  return abstract_asu(
      -y0(cut(vvv(1,0,1),0))
    & p0
    & m1/4
    & cut(vvv(0,1,1),r1/4)
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -F 4ud 2vw 3
abstract::ptr  asu_228 ()
{
  return abstract_asu(
      -y0(cut(vvv(-1,0,1),r1/4))
    & p0(cut(vvv(0,-1,1),0))
    & m1/4
    & cut(vvv(0,1,1),r1/4)(cut(vvv(0,-1,1),0)(x8))
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -I 4 2 3
abstract::ptr  asu_229 ()
{
  return abstract_asu(
      z0
    & p0
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
  );
}

/////   Hall:  -I 4bd 2c 3
abstract::ptr  asu_230 ()
{
  return abstract_asu(
      x8(cut(vvv(0,1,-1),r1/4) & cut(vvv(0,-1,-1),r1/4))
    & ~x8(y0(-z4))
    & y8(cut(vvv(1,0,1),-r1/4))
    & ~y8(cut(vvv(1,0,-1),r1/4))
    & z4(y0)
    & cut(vvv(-1,0,1),0)
    & cut(vvv(1,0,1),0)(-z0)
    & cut(vvv(0,-1,1),0)(cut(vvv(1,0,-1),0))
    & cut(vvv(0,1,1),0)
  );
}



asu_func asu_table[230] = {
  asu_001, asu_002, asu_003, asu_004, asu_005, asu_006, asu_007, asu_008, 
  asu_009, asu_010, asu_011, asu_012, asu_013, asu_014, asu_015, asu_016, 
  asu_017, asu_018, asu_019, asu_020, asu_021, asu_022, asu_023, asu_024, 
  asu_025, asu_026, asu_027, asu_028, asu_029, asu_030, asu_031, asu_032, 
  asu_033, asu_034, asu_035, asu_036, asu_037, asu_038, asu_039, asu_040, 
  asu_041, asu_042, asu_043, asu_044, asu_045, asu_046, asu_047, asu_048, 
  asu_049, asu_050, asu_051, asu_052, asu_053, asu_054, asu_055, asu_056, 
  asu_057, asu_058, asu_059, asu_060, asu_061, asu_062, asu_063, asu_064, 
  asu_065, asu_066, asu_067, asu_068, asu_069, asu_070, asu_071, asu_072, 
  asu_073, asu_074, asu_075, asu_076, asu_077, asu_078, asu_079, asu_080, 
  asu_081, asu_082, asu_083, asu_084, asu_085, asu_086, asu_087, asu_088, 
  asu_089, asu_090, asu_091, asu_092, asu_093, asu_094, asu_095, asu_096, 
  asu_097, asu_098, asu_099, asu_100, asu_101, asu_102, asu_103, asu_104, 
  asu_105, asu_106, asu_107, asu_108, asu_109, asu_110, asu_111, asu_112, 
  asu_113, asu_114, asu_115, asu_116, asu_117, asu_118, asu_119, asu_120, 
  asu_121, asu_122, asu_123, asu_124, asu_125, asu_126, asu_127, asu_128, 
  asu_129, asu_130, asu_131, asu_132, asu_133, asu_134, asu_135, asu_136, 
  asu_137, asu_138, asu_139, asu_140, asu_141, asu_142, asu_143, asu_144, 
  asu_145, asu_146, asu_147, asu_148, asu_149, asu_150, asu_151, asu_152, 
  asu_153, asu_154, asu_155, asu_156, asu_157, asu_158, asu_159, asu_160, 
  asu_161, asu_162, asu_163, asu_164, asu_165, asu_166, asu_167, asu_168, 
  asu_169, asu_170, asu_171, asu_172, asu_173, asu_174, asu_175, asu_176, 
  asu_177, asu_178, asu_179, asu_180, asu_181, asu_182, asu_183, asu_184, 
  asu_185, asu_186, asu_187, asu_188, asu_189, asu_190, asu_191, asu_192, 
  asu_193, asu_194, asu_195, asu_196, asu_197, asu_198, asu_199, asu_200, 
  asu_201, asu_202, asu_203, asu_204, asu_205, asu_206, asu_207, asu_208, 
  asu_209, asu_210, asu_211, asu_212, asu_213, asu_214, asu_215, asu_216, 
  asu_217, asu_218, asu_219, asu_220, asu_221, asu_222, asu_223, asu_224, 
  asu_225, asu_226, asu_227, asu_228, asu_229, asu_230, 
};

}}}

