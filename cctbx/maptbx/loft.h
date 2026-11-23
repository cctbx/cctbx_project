#ifndef CCTBX_MAPTBX_LOFT_H
#define CCTBX_MAPTBX_LOFT_H

#include <cctbx/miller.h>
#include <cctbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/constants.h>
#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>

namespace cctbx { namespace maptbx {

namespace loft_detail {
  struct less_by_key {
    std::vector<long> const& key;
    less_by_key(std::vector<long> const& key_) : key(key_) {}
    bool operator()(std::size_t a, std::size_t b) const {
      return key[a] < key[b];
    }
  };
}

/*

LoFT - Local Fourier Transform.

Authors: Pavel Afonine, Alexandre Urzhumtsev

First release: mid-July 2023. Refactored November 2026 into a general
library function (any unit cell, any set of Miller indices, region given as
a list of grid points).

Structure factors of a map that is one inside a region and zero elsewhere,

  F(h,k,l) = v * sum_p exp(2 pi i (h jx/nx + k jy/ny + l jz/nz)),

where the sum is over the grid points p=(jx,jy,jz) of the region and v is the
volume per grid point. The cost is proportional to the number of region
points times the number of reflections, so for a region that is a small
fraction of the unit cell this is much faster than FFT of the whole cell.
The sum is carried out one dimension at a time (rows along z, planes of
rows, planes), which makes the cost per reflection proportional to the
number of grid planes the region spans, not to its number of points.

No symmetry is used: the list of points must contain all symmetry copies of
the region (the P1 content of the unit cell), and values are computed for
exactly the Miller indices given. The values are exact (trigonometric
factors come from tables exp(2 pi i m/n), m=0..n-1, one per axis).

points: 1d indices into a grid of size n_real=(nx,ny,nz), last index fastest
  (jx*ny*nz + jy*nz + jz), as flex.grid does. Any order.
*/

inline
af::shared<std::complex<double> >
loft(
  af::const_ref<std::size_t>             const& points,
  af::tiny<int, 3>                       const& n_real,
  af::const_ref<cctbx::miller::index<> > const& miller_indices,
  double                                 const& volume_per_point)
{
  typedef std::complex<double> cmplx;
  int nx = n_real[0], ny = n_real[1], nz = n_real[2];
  CCTBX_ASSERT(nx > 0 && ny > 0 && nz > 0);
  std::size_t n_grid = static_cast<std::size_t>(nx) * ny * nz;
  std::size_t nyz = static_cast<std::size_t>(ny) * nz;
  std::size_t n_refl = miller_indices.size();
  af::shared<cmplx> result(n_refl, cmplx(0, 0));
  if(points.size() == 0 || n_refl == 0) return result;
  // Points sorted: planes of equal jx, made of rows of equal (jx,jy)
  std::vector<std::size_t> pts(points.begin(), points.end());
  std::sort(pts.begin(), pts.end());
  CCTBX_ASSERT(pts.back() < n_grid);
  std::vector<int> jz(pts.size());
  std::vector<int> row_y;                // jy of each row
  std::vector<std::size_t> row_start;    // first point of each row
  std::vector<int> plane_x;              // jx of each plane
  std::vector<std::size_t> plane_start;  // first row of each plane
  long prev_x = -1, prev_xy = -1;
  for(std::size_t i = 0; i < pts.size(); i++) {
    std::size_t p = pts[i];
    long x = static_cast<long>(p / nyz);
    long xy = static_cast<long>(p / nz);
    if(xy != prev_xy) {
      if(x != prev_x) {
        plane_x.push_back(static_cast<int>(x));
        plane_start.push_back(row_y.size());
        prev_x = x;
      }
      row_y.push_back(static_cast<int>((p / nz) % ny));
      row_start.push_back(i);
      prev_xy = xy;
    }
    jz[i] = static_cast<int>(p % nz);
  }
  row_start.push_back(pts.size());
  plane_start.push_back(row_y.size());
  std::size_t n_rows = row_y.size();
  std::size_t n_planes = plane_x.size();
  // Tables exp(2 pi i m/n)
  std::vector<cmplx> tx(nx), ty(ny), tz(nz);
  double two_pi = scitbx::constants::two_pi;
  for(int m = 0; m < nx; m++) tx[m] = std::polar(1., two_pi * m / nx);
  for(int m = 0; m < ny; m++) ty[m] = std::polar(1., two_pi * m / ny);
  for(int m = 0; m < nz; m++) tz[m] = std::polar(1., two_pi * m / nz);
  // Reflections ordered by (l,k), indices taken modulo the grid
  std::vector<int> hm(n_refl), km(n_refl), lm(n_refl);
  std::vector<std::size_t> order(n_refl);
  for(std::size_t i = 0; i < n_refl; i++) {
    cctbx::miller::index<> const& hkl = miller_indices[i];
    hm[i] = ((hkl[0] % nx) + nx) % nx;
    km[i] = ((hkl[1] % ny) + ny) % ny;
    lm[i] = ((hkl[2] % nz) + nz) % nz;
    order[i] = i;
  }
  std::vector<long> key(n_refl);
  for(std::size_t i = 0; i < n_refl; i++) {
    key[i] = static_cast<long>(lm[i]) * ny + km[i];
  }
  std::sort(order.begin(), order.end(),
    loft_detail::less_by_key(key));
  // The sum over points is done one dimension at a time:
  //   z_row(l)     = sum over the points of a row of exp(2 pi i l jz/nz)
  //   y_plane(k,l) = sum over the rows of a plane of
  //                  exp(2 pi i k jy/ny) * z_row(l)
  //   F(h,k,l)     = sum over planes of exp(2 pi i h jx/nx) * y_plane(k,l)
  std::vector<cmplx> ez(nz), ey(ny), z_row(n_rows), y_plane(n_planes);
  int l_cur = -1, k_cur = -1;
  for(std::size_t io = 0; io < n_refl; io++) {
    std::size_t ir = order[io];
    if(lm[ir] != l_cur) {
      l_cur = lm[ir];
      k_cur = -1;
      int m = 0;
      for(int j = 0; j < nz; j++) {
        ez[j] = tz[m];
        m += l_cur;
        if(m >= nz) m -= nz;
      }
      for(std::size_t r = 0; r < n_rows; r++) {
        double zr = 0, zi = 0;
        for(std::size_t i = row_start[r]; i < row_start[r+1]; i++) {
          cmplx const& e = ez[jz[i]];
          zr += e.real();
          zi += e.imag();
        }
        z_row[r] = cmplx(zr, zi);
      }
    }
    if(km[ir] != k_cur) {
      k_cur = km[ir];
      int m = 0;
      for(int j = 0; j < ny; j++) {
        ey[j] = ty[m];
        m += k_cur;
        if(m >= ny) m -= ny;
      }
      for(std::size_t p = 0; p < n_planes; p++) {
        cmplx s(0, 0);
        for(std::size_t r = plane_start[p]; r < plane_start[p+1]; r++) {
          s += ey[row_y[r]] * z_row[r];
        }
        y_plane[p] = s;
      }
    }
    cmplx f(0, 0);
    long h = hm[ir];
    for(std::size_t p = 0; p < n_planes; p++) {
      f += tx[(h * plane_x[p]) % nx] * y_plane[p];
    }
    result[ir] = f * volume_per_point;
  }
  return result;
}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_LOFT_H
