// $Id$

#include <cctbx/uctbx.h>
#include <cctbx/basic/matrixlite.h>
#include <cctbx/basic/define_range.h>

namespace { // Helper functions in anonymous namespace.

  using namespace cctbx;
  using uctbx::Vec3;
  using uctbx::Mx33;

  const double EpsPI = 1.e-6; // ARBITRARY

  inline double sinC(double arg) {
    if (constants::pi_2 - EpsPI <= arg && arg <= constants::pi_2 + EpsPI) {
      return 1.;
    }
    return std::sin(arg);
  }

  inline double cosC(double arg) {
    if (constants::pi_2 - EpsPI <= arg && arg <= constants::pi_2 + EpsPI) {
      return 0.;
    }
    return std::cos(arg);
  }

  inline double DotG(const Vec3& u, const Mx33& G, const Vec3& v) {
    return   u[0] * (G[0] * v[0] + G[1] * v[1] + G[2] * v[2])
           + u[1] * (G[3] * v[0] + G[4] * v[1] + G[5] * v[2])
           + u[2] * (G[6] * v[0] + G[7] * v[1] + G[8] * v[2]);
  }

  Vec3 CrossG(const double sqrtdetG, const Mx33& G,
              const Vec3& r, const Vec3& s) {
    Vec3 Gr, Gs, rxs;
    MatrixLite::multiply<double>(G.elems, r.elems, 3, 3, 1, Gr.elems);
    MatrixLite::multiply<double>(G.elems, s.elems, 3, 3, 1, Gs.elems);
    rxs[0] = sqrtdetG * (Gr[1] * Gs[2] - Gs[1] * Gr[2]);
    rxs[1] = sqrtdetG * (Gr[2] * Gs[0] - Gs[2] * Gr[0]);
    rxs[2] = sqrtdetG * (Gr[0] * Gs[1] - Gs[0] * Gr[1]);
    return rxs;
  }

  Mx33 ConstructMetricalMatrix(const Vec3& Len, const Vec3& cosAng)
  {
    Mx33 G;
    rangei(3) G[i * 4] = Len[i] * Len[i];
    G[1] = G[3] = Len[0] * Len[1] * cosAng[2];
    G[2] = G[6] = Len[0] * Len[2] * cosAng[1];
    G[5] = G[7] = Len[1] * Len[2] * cosAng[0];
    return G;
  }

  bool approx_equal(double a, double b, double tolerance = 1.e-6) {
    if (a == b) return true;
    if (a == -b) return false;
    double n = a - b;
    double d = a + b;
    if (n < 0.) n = -n;
    if (d < 0.) d = -d;
    if (n < d * tolerance) return true;
    return false;
  }
}

namespace uctbx {

  void UnitCell::SetVolume()
  {
    /* V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2
                              + 2 * cos(alpha) * cos(beta) * cos(gamma))
     */
    double D = 1.;
    rangei(3) D -= cosAng[i] * cosAng[i];
    D += 2. * cosAng[0] * cosAng[1] * cosAng[2];
    if (D < 0.) throw corrupt_unit_cell_parameters;

        Vol = Len[0] * Len[1] * Len[2] * sqrt(D);
    if (Vol <= 0.) throw corrupt_unit_cell_parameters;
  }

  void UnitCell::SetReciprocal()
  {
    // Transformation Lattice Constants -> Reciprocal Lattice Constants
    // after Kleber, W., 17. Aufl., Verlag Technik GmbH Berlin 1990, P.352

    int i;
    for(i=0;i<3;i++) R_Len[i] =   Len[(i + 1) % 3]
                                * Len[(i + 2) % 3]
                                * sinAng[i] / Vol;
    for(i=0;i<3;i++) R_cosAng[i] =   (  cosAng[(i + 1) % 3]
                                      * cosAng[(i + 2) % 3]
                                      - cosAng[i])
                                   / (  sinAng[(i + 1) % 3]
                                      * sinAng[(i + 2) % 3]);
    for(i=0;i<3;i++) R_Ang[i]    = std::acos(R_cosAng[i]);
    for(i=0;i<3;i++) R_sinAng[i] = sinC(R_Ang[i]);
    for(i=0;i<3;i++) R_cosAng[i] = cosC(R_Ang[i]);
  }

  void UnitCell::SetOrthAndFracMatrix()
  {
    // Crystallographic Basis: D = {a,b,c}
    // Cartesian Basis:        C = {i,j,k}
    //
    // PDB convention:
    //   i || a
    //   j is in (a,b) plane
    //   k = i x j

    double s1rca2 = sqrt(1. - R_cosAng[0] * R_cosAng[0]);
    if (s1rca2 == 0.) throw corrupt_unit_cell_parameters;

    // fractional to cartesian
    Orth[0] =  Len[0];
    Orth[1] =  cosAng[2] * Len[1];
    Orth[2] =  cosAng[1] * Len[2];
    Orth[3] =  0.;
    Orth[4] =  sinAng[2] * Len[1];
    Orth[5] = -sinAng[1] * R_cosAng[0] * Len[2];
    Orth[6] =  0.;
    Orth[7] =  0.;
    Orth[8] =  sinAng[1] * Len[2] * s1rca2;

    // cartesian to fractional
    Frac[0] =  1. / Len[0];
    Frac[1] = -cosAng[2] / (sinAng[2] * Len[0]);
    Frac[2] = -(  cosAng[2] * sinAng[1] * R_cosAng[0]
                + cosAng[1] * sinAng[2])
              / (sinAng[1] * s1rca2 * sinAng[2] * Len[0]);
    Frac[3] =  0.;
    Frac[4] =  1. / (sinAng[2] * Len[1]);
    Frac[5] =  R_cosAng[0] / (s1rca2 * sinAng[2] * Len[1]);
    Frac[6] =  0.;
    Frac[7] =  0.;
    Frac[8] =  1. / (sinAng[1] * s1rca2 * Len[2]);
  }

  void UnitCell::SetMetricalMatrices()
  {
    G = ConstructMetricalMatrix(Len, cosAng);
    R_G = ConstructMetricalMatrix(R_Len, R_cosAng);
  }

  void UnitCell::Initialize()
  {
    int i;
    for(i=0;i<3;i++) if (Len[i] <= 0.) throw corrupt_unit_cell_parameters;
    for(i=0;i<3;i++) if (Ang[i] <= 0. || Ang[i] >= constants::pi)
      throw corrupt_unit_cell_parameters;

    for(i=0;i<3;i++) sinAng[i] = std::sin(Ang[i]);
    for(i=0;i<3;i++) cosAng[i] = std::cos(Ang[i]);

    for(i=0;i<3;i++) if (sinAng[i] == 0.) throw corrupt_unit_cell_parameters;

    SetVolume();
    SetReciprocal();
    SetMetricalMatrices();
    SetOrthAndFracMatrix();
  }

  UnitCell::UnitCell(const uc_params& ucp) {
    int i;
    for(i=0;i<3;i++) Len[i] = ucp.Len(i);
    for(i=0;i<3;i++) Ang[i] = deg_as_rad(ucp.Ang(i));
    Initialize();
  }

  UnitCell::UnitCell(const Mx33& MetricalMatrix)
  {
    for (int i = 0; i < 9; i += 4) {
      if (MetricalMatrix[i] <= 0.) throw corrupt_metrical_matrix;
    }
    if (   ! approx_equal(MetricalMatrix[1], MetricalMatrix[3])
        || ! approx_equal(MetricalMatrix[2], MetricalMatrix[6])
        || ! approx_equal(MetricalMatrix[5], MetricalMatrix[7]))
      throw corrupt_metrical_matrix;
    Len[0] = sqrt(MetricalMatrix[0]);
    Len[1] = sqrt(MetricalMatrix[4]);
    Len[2] = sqrt(MetricalMatrix[8]);
    Ang[0] = acos(MetricalMatrix[5] / Len[1] / Len[2]);
    Ang[1] = acos(MetricalMatrix[2] / Len[2] / Len[0]);
    Ang[2] = acos(MetricalMatrix[1] / Len[0] / Len[1]);
    try {
      Initialize();
    }
    catch (error) {
      throw corrupt_metrical_matrix;
    }
  }

  uc_params UnitCell::getParameters(bool reciprocal) const
  {
    uc_params ucp;
    if (reciprocal == false) ucp = uc_params(Len, Ang);
    else                     ucp = uc_params(R_Len, R_Ang);
    rangei(3) ucp.Ang()[i] = rad_as_deg(ucp.Ang()[i]);
    return ucp;
  }

  Miller::Index UnitCell::MaxMillerIndices(double dmin) const
  {
    Miller::Index MaxMIx;
    int i, j;
    for(i=0;i<3;i++) {
      Vec3 u, v, uxv;
      for(j=0;j<3;j++) u[j] = 0.; u[(i + 1) % 3] = 1.;
      for(j=0;j<3;j++) v[j] = 0.; v[(i + 2) % 3] = 1.;
      uxv = CrossG(1., R_G, u, v); // Since length of uxv is not used
                                   //   sqrt(det(G)) is set to 1
      double uxv2 = DotG(uxv, R_G, uxv);
      MaxMIx[i] = (int) (uxv[i] / std::sqrt(uxv2) / dmin + 1.e-4);
                                                        // ARBITRARY
    }
    return MaxMIx;
  }

  UnitCell UnitCell::ChangeBasis(const Mx33& InvCBMxR, double RBF) const
  {
    Mx33 R = InvCBMxR;
    if (RBF != 0.) rangei(9) R[i] /= RBF;
    Mx33 RtGR = getRtGR(G, R);
    return UnitCell(RtGR);
  }

  UnitCell UnitCell::ChangeBasis(const sgtbx::RotMx& InvCBMxR) const
  {
    return ChangeBasis(InvCBMxR.as_array(static_cast<double>(0)), 1.);
  }

  std::ostream& UnitCell::print(std::ostream& os) const {
    os << Len[0] << " " << Len[1] << " " << Len[2] << " "
       << rad_as_deg(Ang[0]) << " "
       << rad_as_deg(Ang[1]) << " "
       << rad_as_deg(Ang[2]);
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const UnitCell& uc) {
    return uc.print(os);
  }

} // namespace uctbx
