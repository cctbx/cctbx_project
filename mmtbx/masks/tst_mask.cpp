#ifdef NDEBUG
#undef NDEBUG
#endif

#include <iostream>
#include <vector>
#include <set>
#include <memory>
#include <algorithm>

#include <cassert>

#include <iotbx/pdb/input.h>
#include <cctbx/sgtbx/space_group.h>
#include <mmtbx/masks/atom_mask.h>
#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>
#include <scitbx/fftpack/real_to_complex_3d.h>
#include <cctbx/maptbx/copy.h>
#include <cctbx/maptbx/structure_factors.h>
#include <cctbx/miller.h>
#include <scitbx/array_family/flex_types.h>

//
//  Compilation flags:
// -I ~/Phenix/cctbx_sources  -I ~/Phenix/cctbx_build/include -I ~/Phenix/cctbx_sources/boost -L$HOME/Phenix/cctbx_build/lib -lcctbx_sgtbx_asu -lcctbx
//

using namespace cctbx;
using namespace scitbx;
using namespace scitbx::af;

int main(int argc, const char* argv[])
{
  using namespace mmtbx::masks;
  using namespace cctbx::sgtbx::asu;

  try {
    std::string spgr("P 21 21 21");
    double resolution = 2.0;
    int method = 2;
    if( argc>1 )
      spgr = argv[1];
    if( argc > 2 ) {
      resolution = std::atof(argv[2]);
    }
    if( argc>3 )
      method = std::atoi( argv[3] );

    cctbx::sgtbx::space_group_symbols symbol(spgr) ;
    cctbx::sgtbx::space_group grp ( symbol );
    std::cout << "Space group = " << spgr << "  Hall: " << grp.type().hall_symbol()
      << "  order= " << grp.order_z() << std::endl;
    double aparams[6] = {20.0, 20.0, 20.0, 90.0, 90.0, 90.0};
    scitbx::af::small<double,6> params(aparams, aparams+6); //;

    cctbx::uctbx::unit_cell cell(params);

    mmtbx::masks::coord_array_t coords;
    mmtbx::masks::double_array_t radii;

    mmtbx::masks::atom_mask msk(
      cell,
      grp,
      resolution,
      method);
    msk.compute(coords, radii);
    scitbx::af::shared< cctbx::miller::index<> > indices;
    for(size_t i=1; i<10; ++i)
    {
      cctbx::sg_vec3 h(i, i, i);
      indices.push_back( h );
    }
    std::cout << "number of indices = " << indices.size() << std::endl;
    scitbx::af::shared< std::complex<double> >  fc;
    fc = msk.structure_factors(indices.const_ref());
    assert( fc.size() == indices.size() );
    for( size_t i=0; i<fc.size() && i<10U; ++i)
      std::cout << fc[i] << std::endl;
    // iotbx::pdb::input inp("tst.pdb");
    // iotbx::pdb::hierarchy::root pdbh = inp.construct_hierarchy();
    // scitbx::af::shared< iotbx::pdb::hierarchy::atom > atoms = pdbh.atoms();

  }
  catch(const cctbx::error &err )
  {
    std::cerr << "\n" << err.what() << std::endl;
    assert( 0 );
    return 1;
  }
  catch( const std::exception &err )
  {
    std::cerr << "\n===== ERROR =====\n" << err.what() << std::endl;
    assert( 0 );
    return 1;
  }
  return 0;
}
