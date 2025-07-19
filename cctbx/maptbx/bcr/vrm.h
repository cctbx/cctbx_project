#ifndef VRM_H
#define VRM_H

#include <boost/python/list.hpp>
#include <cctbx/maptbx/bcr/bcr.h>

namespace cctbx { namespace maptbx {

using scitbx::constants::pi;
using scitbx::constants::four_pi_sq;
using scitbx::constants::four_pi;

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class OmegaMap {

public:
  FloatType acell, bcell, ccell, alpha, beta,  gamma;
  scitbx::mat3<double> OrthMatrix, DeortMatrix;
  int Mx, My, Mz;
  int Sx, Sy, Sz;
  int Nx, Ny, Nz;
  int Fx, Fy, Fz;
  FloatType StepX, StepY, StepZ;
  FloatType orthxx, orthxy, orthxz;
  FloatType orthyy, orthyz;
  FloatType orthzz;
  FloatType dortxx, dortxy, dortxz;
  FloatType dortyy, dortyz;
  FloatType dortzz;
  FloatType StepXX, StepXY, StepXZ;
  FloatType StepYY, StepYZ;
  FloatType StepZZ;
  FloatType StepXX2, StepXXS;
  FloatType StepYY2, StepYYS;
  FloatType StepZZ2, StepZZS;
  FloatType StepXXS2, StepYYS2, StepZZS2;
  FloatType cosalp, cosbet, cosgam;
  FloatType vol0;
  FloatType RprojX;
  FloatType RprojY;
  FloatType RprojZ;
  af::versa<FloatType, af::c_grid<3> > map;
  boost::python::list bcr_scatterers;

  OmegaMap(
    af::tiny<int, 3> const& Ncrs,
    af::tiny<int, 3> const& Scrs,
    af::tiny<int, 3> const& Nxyz,
    cctbx::uctbx::unit_cell const& unit_cell,
    const boost::python::list & bcr_scatterers_)
  :
  bcr_scatterers(bcr_scatterers_), // is this doing copy?
  map(af::c_grid<3>(Nxyz[0],Nxyz[1],Nxyz[2]))
  {
    scitbx::af::tiny<FloatType, 6> ucp = unit_cell.parameters();
    FloatType acell = ucp[0];
    FloatType bcell = ucp[1];
    FloatType ccell = ucp[2];
    FloatType alpha = ucp[3];
    FloatType beta  = ucp[4];
    FloatType gamma = ucp[5];
    OrthMatrix  = unit_cell.orthogonalization_matrix();
    DeortMatrix = unit_cell.fractionalization_matrix();
    Mx = Ncrs[0];
    My = Ncrs[1];
    Mz = Ncrs[2];
    Sx = Scrs[0];
    Sy = Scrs[1];
    Sz = Scrs[2];
    Nx = Nxyz[0];
    Ny = Nxyz[1];
    Nz = Nxyz[2];
    Fx = Sx + Mx;
    Fy = Sy + My;
    Fz = Sz + Mz;
    StepX = 1.0/static_cast<FloatType>(Nx);
    StepY = 1.0/static_cast<FloatType>(Ny);
    StepZ = 1.0/static_cast<FloatType>(Nz);
    orthxx = OrthMatrix[0];
    orthxy = OrthMatrix[3];
    orthxz = OrthMatrix[6];
    orthyy = OrthMatrix[4];
    orthyz = OrthMatrix[7];
    orthzz = OrthMatrix[8];
    dortxx = DeortMatrix[0];
    dortxy = DeortMatrix[3];
    dortxz = DeortMatrix[6];
    dortyy = DeortMatrix[4];
    dortyz = DeortMatrix[7];
    dortzz = DeortMatrix[8];
    StepXX = orthxx * StepX;
    StepXY = orthxy * StepY;
    StepXZ = orthxz * StepZ;
    StepYY = orthyy * StepY;
    StepYZ = orthyz * StepZ;
    StepZZ  = orthzz * StepZ;
    StepXX2 = StepXX * 2;
    StepXXS = StepXX * StepXX;
    StepYY2 = StepYY * 2;
    StepYYS = StepYY * StepYY;
    StepZZ2 = StepZZ * 2;
    StepZZS = StepZZ * StepZZ;
    StepXXS2 = StepXXS * 2;
    StepYYS2 = StepYYS * 2;
    StepZZS2 = StepZZS * 2;
    cosalp = std::cos(alpha);
    cosbet = std::cos(beta);
    cosgam = std::cos(gamma);
    vol0 = std::sqrt(
      1.-cosalp*cosalp-cosbet*cosbet-cosgam*cosgam+2.*cosalp*cosbet*cosgam);
    RprojX = std::sin(alpha) / (vol0 * acell);
    RprojY = std::sin(beta)  / (vol0 * bcell);
    RprojZ = std::sin(gamma) / (vol0 * ccell);
  }

  af::versa<FloatType, af::c_grid<3> > compute_map(bool arg_value) {
    for(std::size_t i=0; i<boost::python::len(bcr_scatterers); i++) {
      bcr_scatterer<FloatType> bcrs =
         boost::python::extract<bcr_scatterer<FloatType> >(bcr_scatterers[i])();
      FloatType RadAtom = bcrs.radius;
      FloatType RadAtom2  = RadAtom   * RadAtom;
      FloatType RadAtomX  = RadAtom   * RprojX;
      FloatType RadAtomY  = RadAtom   * RprojY;
      FloatType RadAtomZ  = RadAtom   * RprojZ;
      cctbx::cartesian<> r = bcrs.site_cart;
      FloatType xat = r[0];
      FloatType yat = r[1];
      FloatType zat = r[2];
      FloatType cat = bcrs.occ;
      FloatType bat = bcrs.u_iso;
      FloatType xfrac = dortxx * xat + dortxy * yat + dortxz * zat;
      FloatType yfrac =                dortyy * yat + dortyz * zat;
      FloatType zfrac =                               dortzz * zat;
      auto box = AtomBox(xfrac, yfrac, zfrac,
        RadAtomX,RadAtomY,RadAtomZ,StepX,StepY,StepZ,Sx,Sy,Sz,Fx,Fy,Fz);
      int Kx1       = std::get<0>(box);
      int Kx2       = std::get<1>(box);
      int Ky1       = std::get<2>(box);
      int Ky2       = std::get<3>(box);
      int Kz1       = std::get<4>(box);
      int Kz2       = std::get<5>(box);
      FloatType dxf = std::get<6>(box);
      FloatType dyf = std::get<7>(box);
      FloatType dzf = std::get<8>(box);
      FloatType dx0 = orthxx * dxf + orthxy * dyf + orthxz * dzf;
      FloatType dy0 =                orthyy * dyf + orthyz * dzf;
      FloatType dz0 =                               orthzz * dzf;
      struct AtomGridEntry {
        FloatType r2zyx;
        FloatType rzyx;
        FloatType rzyx2;
        int ix;
        int iy;
        int iz;
      };
      std::vector<AtomGridEntry> AtomGrids;

      af::tiny<int, 3> AtomGrid0(-1, -1, -1);
      int ix0 = -1;

      FloatType dz  = dz0;
      FloatType r2z = dz * dz;
      FloatType cz  = dz * StepZZ2 + StepZZS;

      for (int iz = Kz1; iz < Kz2; ++iz) {
        FloatType dy   = dy0;
        FloatType dt0  = dx0;
        FloatType r2zy = dy * dy + r2z;
        FloatType cy   = dy * StepYY2 + StepYYS;
        for (int iy = Ky1; iy < Ky2; ++iy) {
          FloatType dx    = dt0;
          FloatType r2zyx = dx * dx + r2zy;
          FloatType cx    = dx * StepXX2 + StepXXS;
          for (int ix = Kx1; ix < Kx2; ++ix) {
              if(r2zyx == 0.0) {
                 AtomGrid0 = af::tiny<int, 3>(ix,iy,iz);
                 ix0       = ix;
              }
              else if(r2zyx <= RadAtom2) {
                 FloatType rzyx  = std::sqrt(r2zyx);
                 FloatType rzyx2 = rzyx * 2;
                 AtomGrids.push_back({r2zyx, rzyx, rzyx2, ix, iy, iz});
              }
              r2zyx = r2zyx + cx;
              cx    = cx    + StepXXS2;
              dx    = dx    + StepXX;
          }
          r2zy = r2zy + cy;
          cy   = cy   + StepYYS2;
          dt0  = dt0  + StepXY;
          dy   = dy   + StepYY;
        }
        r2z = r2z + cz;
        cz  = cz  + StepZZS2;
        dx0 = dx0 + StepXZ;
        dy0 = dy0 + StepYZ;
        dz  = dz  + StepZZ;
      }
      dz  = dz0;
      r2z = dz * dz;
      cz  = dz * StepZZ2 + StepZZS;

      int Ngrids = AtomGrids.size();
      af::shared<FloatType> GridValues(Ngrids, 0.0);
      FloatType GridValue0 = 0.0;

      for (size_t iterm = 0; iterm < bcrs.mu.size(); ++iterm) {
        FloatType mu    = bcrs.mu[iterm];
        FloatType nu    = bcrs.nu[iterm];
        FloatType musq  = bcrs.musq[iterm];
        FloatType kappi = bcrs.kappi[iterm];
        FloatType nuatom  = nu + bat;
        FloatType nuatom2 = nuatom + nuatom;
        FloatType fact1   = kappi / std::pow(nuatom2, 1.5);
        FloatType munuat  = mu / nuatom;
        // contribution to the node coinciding with the atomic center
        if(ix0 >= 0) {
          FloatType argg  = musq / nuatom2;
          FloatType fact2 = std::exp(-argg);
          FloatType GridValue0 = GridValue0 + fact1 * fact2;
        }
        // contribution to nodes different from the atomic center
        if(mu == 0.0) {
          for (int ig = 0; ig < Ngrids; ++ig) {
            AtomGridEntry ag = AtomGrids[ig];
            FloatType argg  = ag.r2zyx / nuatom2;
            FloatType fact2 = std::exp(-argg);
            GridValues[ig] = GridValues[ig] + fact1 * fact2;
          }
        }
        else {
          for (int ig = 0; ig < Ngrids; ++ig) {
            AtomGridEntry ag = AtomGrids[ig];
            FloatType tterm   = ag.rzyx2 * munuat;
            FloatType argg  = std::pow(ag.rzyx-mu, 2) / nuatom2;
            FloatType fact2 = std::exp(-argg) * (1.0 - std::exp(-tterm)) / tterm;
            GridValues[ig] = GridValues[ig] + fact1 * fact2;
          }
        }
      }

      for (int ig = 0; ig < Ngrids; ++ig) {
        AtomGridEntry ag = AtomGrids[ig];
        map(ag.iz, ag.iy, ag.ix) += GridValues[ig] * cat;
      }

      ix0 = AtomGrid0[0];
      int iy0 = AtomGrid0[1];
      int iz0 = AtomGrid0[2];
      if(ix0 >= 0) {
        FloatType val = map(iz0, iy0, ix0);
        val += GridValue0 * cat;
        map(iz0, iy0, ix0) = val;
      }

    }

    return map;
  }

  inline std::tuple<int, int, int, int, int, int, double, double, double>
  AtomBox(double xfrac, double yfrac, double zfrac,
          double RadAtomX, double RadAtomY, double RadAtomZ,
          double StepX, double StepY, double StepZ,
          int Sx, int Sy, int Sz,
          int Fx, int Fy, int Fz)
  {
      double x1 = (xfrac - RadAtomX) / StepX;
      double y1 = (yfrac - RadAtomY) / StepY;
      double z1 = (zfrac - RadAtomZ) / StepZ;
      int Kx1 = static_cast<int>(x1);
      int Ky1 = static_cast<int>(y1);
      int Kz1 = static_cast<int>(z1);
      if (x1 >= 0.0) Kx1 += 1;
      if (y1 >= 0.0) Ky1 += 1;
      if (z1 >= 0.0) Kz1 += 1;
      double x2 = (xfrac + RadAtomX) / StepX;
      double y2 = (yfrac + RadAtomY) / StepY;
      double z2 = (zfrac + RadAtomZ) / StepZ;
      int Kx2 = static_cast<int>(x2) + 1;
      int Ky2 = static_cast<int>(y2) + 1;
      int Kz2 = static_cast<int>(z2) + 1;
      Kx1 = std::max(Kx1, Sx);
      Ky1 = std::max(Ky1, Sy);
      Kz1 = std::max(Kz1, Sz);
      Kx2 = std::min(Kx2, Fx);
      Ky2 = std::min(Ky2, Fy);
      Kz2 = std::min(Kz2, Fz);
      double dxf = Kx1 * StepX - xfrac;
      double dyf = Ky1 * StepY - yfrac;
      double dzf = Kz1 * StepZ - zfrac;
      Kx1 -= Sx;
      Ky1 -= Sy;
      Kz1 -= Sz;
      Kx2 -= Sx;
      Ky2 -= Sy;
      Kz2 -= Sz;
      return std::make_tuple(Kx1, Kx2, Ky1, Ky2, Kz1, Kz2, dxf, dyf, dzf);
  }

};


}} // namespace cctbx::maptbx

#endif
