// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/uctbx.h>

namespace { // Helper functions in anonymous namespace.

  using namespace cctbx;

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

  inline
  double
  DotG(const af::double3& u, const af::double9& G, const af::double3& v) {
    return   u[0] * (G[0] * v[0] + G[1] * v[1] + G[2] * v[2])
           + u[1] * (G[3] * v[0] + G[4] * v[1] + G[5] * v[2])
           + u[2] * (G[6] * v[0] + G[7] * v[1] + G[8] * v[2]);
  }

  af::double3 CrossG(const double& sqrtdetG, const af::double9& G,
                     const af::double3& r, const af::double3& s) {
    af::double3 Gr, Gs;
    MatrixLite::multiply<double>(G.begin(), r.begin(), 3, 3, 1, Gr.begin());
    MatrixLite::multiply<double>(G.begin(), s.begin(), 3, 3, 1, Gs.begin());
    return sqrtdetG * MatrixLite::cross_product(Gr, Gs);
  }

  af::double9
  ConstructMetricalMatrix(const af::double3& Len, const af::double3& cosAng)
  {
    af::double9 G;
    for(std::size_t i=0;i<3;i++) G[i * 4] = Len[i] * Len[i];
    G[1] = G[3] = Len[0] * Len[1] * cosAng[2];
    G[2] = G[6] = Len[0] * Len[2] * cosAng[1];
    G[5] = G[7] = Len[1] * Len[2] * cosAng[0];
    return G;
  }

  bool isSymmetric(const af::double9& M, double tolerance = 1.e-6)
  {
    double maxelem = M[0];
    for(int i=1;i<9;i++) if (maxelem < M[i]) maxelem = M[i];
    return    af::approx_equal_scaled(M[1], M[3], maxelem * tolerance)
           && af::approx_equal_scaled(M[2], M[6], maxelem * tolerance)
           && af::approx_equal_scaled(M[5], M[7], maxelem * tolerance);
  }
}

namespace cctbx { namespace uctbx {

  void UnitCell::SetVolume()
  {
    /* V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2
                              + 2 * cos(alpha) * cos(beta) * cos(gamma))
     */
    double D = 1.;
    for(std::size_t i=0;i<3;i++) D -= cosAng[i] * cosAng[i];
    D += 2. * cosAng[0] * cosAng[1] * cosAng[2];
    if (D < 0.) throw corrupt_unit_cell_parameters;

        Vol = Len[0] * Len[1] * Len[2] * std::sqrt(D);
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

    double s1rca2 = std::sqrt(1. - R_cosAng[0] * R_cosAng[0]);
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

  void UnitCell::SetLongestVector2()
  {
    LongestVector2 = 0.;
    int Corner[3];
    for (Corner[0] = 0; Corner[0] <= 1; Corner[0]++)
    for (Corner[1] = 0; Corner[1] <= 1; Corner[1]++)
    for (Corner[2] = 0; Corner[2] <= 1; Corner[2]++) {
      fractional<double> Frac;
      for(std::size_t i=0;i<3;i++) Frac[i] = Corner[i];
      double Cart2 = orthogonalize(Frac).Length2();
      if (LongestVector2 < Cart2) LongestVector2 = Cart2;
    }
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
    SetLongestVector2();
  }

  UnitCell::UnitCell() {
    int i;
    for(i=0;i<3;i++) Len[i] = 1.;
    for(i=0;i<3;i++) Ang[i] = deg_as_rad(90.);
    Initialize();
  }

  UnitCell::UnitCell(const uc_params& ucp) {
    int i;
    for(i=0;i<3;i++) Len[i] = ucp.Len(i);
    for(i=0;i<3;i++) Ang[i] = deg_as_rad(ucp.Ang(i));
    Initialize();
  }

  UnitCell::UnitCell(const af::double9& MetricalMatrix)
  {
    for (int i = 0; i < 9; i += 4) {
      if (MetricalMatrix[i] <= 0.) throw corrupt_metrical_matrix;
    }
    if (!isSymmetric(MetricalMatrix)) throw corrupt_metrical_matrix;
    Len[0] = std::sqrt(MetricalMatrix[0]);
    Len[1] = std::sqrt(MetricalMatrix[4]);
    Len[2] = std::sqrt(MetricalMatrix[8]);
    Ang[0] = std::acos(MetricalMatrix[5] / Len[1] / Len[2]);
    Ang[1] = std::acos(MetricalMatrix[2] / Len[2] / Len[0]);
    Ang[2] = std::acos(MetricalMatrix[1] / Len[0] / Len[1]);
    try {
      Initialize();
    }
    catch (const error&) {
      throw corrupt_metrical_matrix;
    }
  }

  uc_params UnitCell::getParameters(bool reciprocal) const
  {
    uc_params ucp;
    if (reciprocal == false) ucp = uc_params(Len, Ang);
    else                     ucp = uc_params(R_Len, R_Ang);
    for(std::size_t i=0;i<3;i++) ucp.Ang()[i] = rad_as_deg(ucp.Ang()[i]);
    return ucp;
  }

  // XXX design as isSimilarTo(other, len_tolerance, ang_tolerance)
  bool
  UnitCell::isEqual(const UnitCell& other, double tolerance) const
  {
    for(int i=0;i<3;i++) {
      if (!af::approx_equal_unscaled(Len[i],other.Len[i], tolerance))
        return false;
      if (!af::approx_equal_unscaled(Ang[i],other.Ang[i], tolerance))
        return false;
    }
    return true;
  }

  miller::Index UnitCell::MaxMillerIndices(double dmin) const
  {
    miller::Index MaxMIx;
    int i, j;
    for(i=0;i<3;i++) {
      af::double3 u, v, uxv;
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

  af::shared<double>
  UnitCell::Q(const af::shared<miller::Index>& MIx) const
  {
    af::shared<double> result;
    result.reserve(MIx.size());
    for(std::size_t i=0;i<MIx.size();i++) {
      result.push_back(Q(MIx[i]));
    }
    return result;
  }

  double
  UnitCell::max_Q(const af::shared<miller::Index>& MIx) const
  {
    double result = 0;
    for(std::size_t i=0;i<MIx.size();i++) {
      math::update_max(result, Q(MIx[i]));
    }
    return result;
  }

  af::tiny<double, 2>
  UnitCell::min_max_Q(const af::shared<miller::Index>& MIx) const
  {
    af::tiny<double, 2> result(0, 0);
    if (MIx.size()) {
      result.fill(Q(MIx[0]));
      for(std::size_t i=1;i<MIx.size();i++) {
        double q = Q(MIx[i]);
        math::update_min(result[0], q);
        math::update_max(result[1], q);
      }
    }
    return result;
  }

  af::shared<double>
  UnitCell::stol2(const af::shared<miller::Index>& MIx) const
  {
    af::shared<double> result;
    result.reserve(MIx.size());
    for(std::size_t i=0;i<MIx.size();i++) {
      result.push_back(stol2(MIx[i]));
    }
    return result;
  }

  af::shared<double>
  UnitCell::two_stol(const af::shared<miller::Index>& MIx) const
  {
    af::shared<double> result;
    result.reserve(MIx.size());
    for(std::size_t i=0;i<MIx.size();i++) {
      result.push_back(two_stol(MIx[i]));
    }
    return result;
  }

  af::shared<double>
  UnitCell::stol(const af::shared<miller::Index>& MIx) const
  {
    af::shared<double> result;
    result.reserve(MIx.size());
    for(std::size_t i=0;i<MIx.size();i++) {
      result.push_back(stol(MIx[i]));
    }
    return result;
  }

  af::shared<double>
  UnitCell::d(const af::shared<miller::Index>& MIx) const
  {
    af::shared<double> result;
    result.reserve(MIx.size());
    for(std::size_t i=0;i<MIx.size();i++) {
      result.push_back(d(MIx[i]));
    }
    return result;
  }

  af::shared<double>
  UnitCell::two_theta(
    af::shared<miller::Index> MIx, double wavelength, bool deg) const
  {
    af::shared<double> result;
    result.reserve(MIx.size());
    for(std::size_t i=0;i<MIx.size();i++) {
      result.push_back(two_theta(MIx[i], wavelength, deg));
    }
    return result;
  }

  UnitCell UnitCell::ChangeBasis(const af::double9& InvCBMxR, double RBF) const
  {
    af::double9 R = InvCBMxR;
    if (RBF != 0.) for(std::size_t i=0;i<9;i++) R[i] /= RBF;
    af::double9 RtGR = getRtGR(G, R);
    return UnitCell(RtGR);
  }

  UnitCell UnitCell::ChangeBasis(const sgtbx::RotMx& InvCBMxR) const
  {
    return ChangeBasis(InvCBMxR.as_array(double()), 1.);
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

}} // namespace cctbx::uctbx
