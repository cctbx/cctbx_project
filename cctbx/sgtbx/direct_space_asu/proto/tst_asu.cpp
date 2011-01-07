#include <iostream>
#include <vector>
#include <set>
#include <memory>

#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>

//
//  Compilation flags:
// -I ~/Phenix/cctbx_sources  -I ~/Phenix/cctbx_build/include -I ~/Phenix/cctbx_sources/boost -L$HOME/Phenix/cctbx_build/lib -lcctbx_sgtbx_asu -lcctbx
//

using namespace cctbx::sgtbx::asu;
using namespace cctbx::sgtbx;


inline double get_double(rational_t r)
{
  return double(r.numerator())/r.denominator();
}

std::size_t loop_over_grid_points(const direct_space_asu &a, unsigned n_)
{
  register std::size_t result = 0;
  register const int n = n_;

  a.show_comprehensive_summary(std::cout);
  std::cout <<"\n";

  rvector3_t p, mn, mx;
  a.box_corners(mn, mx);
  std::cout << "low corner = " << mn << "   high corner = " << mx << std::endl;
  rvector3_t box = mx - mn;
  rvector3_t step = box/rational_t(n);
  std::cout << "step = " << step << std::endl;

  for(rational_t i=mn[0]; i<=mx[0]; i += step[0]) {
    p[0] = i;
  for(rational_t j=mn[1]; j<=mx[1]; j += step[1]) {
    p[1] = j;
  for(rational_t k=mn[2]; k<=mx[2]; k += step[2]) {
    p[2] = k;
    if (a.is_inside(p)) result += 1;
  }}}

  const char *asun = "n asu points:  ";
  const char *vol = "   volume: ";
  space_group grp(a.hall_symbol);
  std::cout << asun << result << vol << get_double(box[0]*box[1]*box[2]) * result/(double(n)*n*n)
    <<"   expected volume= "<< 1.0/grp.order_z() << std::endl;
  return result;
}


int main(int argc, const char* argv[])
{
  try {

    unsigned n = 100;
    if (argc > 1) {
      n = std::atoi(argv[1]);
    }
    std::string spgr;
    if( argc>2 )
      spgr = argv[2];

    if( spgr.empty()  )
      spgr = "P 21 21 21";

    cctbx::sgtbx::space_group_type grp( spgr );
    std::cout << "Space group: " << spgr << "  number: "<< grp.number() << "  hall: " << grp.hall_symbol() << std::endl;

    // direct_space_asu asu( grp );
    direct_space_asu asu( spgr );

    size_t ins = 0;

    std::cout << "n: " << n << std::endl << std::flush;

    ins = loop_over_grid_points(asu, n);

    asu.shape_only();
    std::cout << "\nAfter shape_only\n";

    ins = loop_over_grid_points(asu, n);
  }
  catch( const std::exception &err )
  {
    std::cerr << "\n===== ERROR =====\n" << err.what() << std::endl;
  }
  return 0;
}

