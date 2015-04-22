#ifndef CCTBX_GEOMETRY_RESTRAINTS_PARALLELITY_H
#define CCTBX_GEOMETRY_RESTRAINTS_PARALLELITY_H

#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/geometry_restraints/utils.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/misc_functions.h>

namespace cctbx { namespace geometry_restraints {

  //! Grouping of indices into array of sites (i_seqs) and weights.
  struct parallelity_proxy
  {
    //! Support for shared_proxy_select.
    typedef af::shared<std::size_t> i_seqs_type;

    //! Default constructor. Some data members are not initialized!
    parallelity_proxy() {}

    //! Constructor.
    parallelity_proxy(
      i_seqs_type const& i_seqs_,
      i_seqs_type const& j_seqs_,
      double weight_,
      double target_angle_deg_ = 0,
      double slack_ = 0,
      double limit_ = 1,
      bool top_out_ = false,
      unsigned char origin_id_ = 0)
    :
      i_seqs(i_seqs_),
      j_seqs(j_seqs_),
      weight(weight_),
      target_angle_deg(target_angle_deg_),
      slack(slack_),
      limit(limit_),
      top_out(top_out_),
      origin_id(origin_id_)
    {
      CCTBX_ASSERT(i_seqs.size() > 2);
      CCTBX_ASSERT(j_seqs.size() > 2);
      CCTBX_ASSERT(weight > 0);
      CCTBX_ASSERT(slack >= 0);
      CCTBX_ASSERT(slack <= 90);
      CCTBX_ASSERT(limit >= 1);
    }

    //! Constructor.
    parallelity_proxy(
      i_seqs_type const& i_seqs_,
      i_seqs_type const& j_seqs_,
      optional_container<af::shared<sgtbx::rt_mx> > const& sym_ops_,
      double weight_,
      double target_angle_deg_ = 0,
      double slack_ = 0,
      double limit_ = 1,
      bool top_out_ = false,
      unsigned char origin_id_ = 0)
    :
      i_seqs(i_seqs_),
      j_seqs(j_seqs_),
      sym_ops(sym_ops_),
      weight(weight_),
      target_angle_deg(target_angle_deg_),
      slack(slack_),
      limit(limit_),
      top_out(top_out_),
      origin_id(origin_id_)
    {
      CCTBX_ASSERT(i_seqs.size() > 2);
      CCTBX_ASSERT(j_seqs.size() > 2);
      CCTBX_ASSERT(weight > 0);
      CCTBX_ASSERT(slack >= 0);
      CCTBX_ASSERT(slack <= 90);
      CCTBX_ASSERT(limit >= 1);
      if ( sym_ops.get() != 0 ) {
        CCTBX_ASSERT(sym_ops.get()->size() == i_seqs.size());
      }
    }

    //! Support for proxy_select (and similar operations).
    parallelity_proxy(
      i_seqs_type const& i_seqs_,
      i_seqs_type const& j_seqs_,
      parallelity_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      j_seqs(j_seqs_),
      sym_ops(proxy.sym_ops),
      weight(proxy.weight),
      target_angle_deg(proxy.target_angle_deg),
      slack(proxy.slack),
      limit(proxy.limit),
      top_out(proxy.top_out),
      origin_id(proxy.origin_id)
    {
      CCTBX_ASSERT(i_seqs.size() > 2);
      CCTBX_ASSERT(j_seqs.size() > 2);
      CCTBX_ASSERT(weight > 0);
      CCTBX_ASSERT(slack >= 0);
      CCTBX_ASSERT(slack <= 90);
      CCTBX_ASSERT(limit >= 1);
      if ( sym_ops.get() != 0 ) {
        CCTBX_ASSERT(sym_ops.get()->size() == i_seqs.size());
      }
    }


    parallelity_proxy
    scale_weight(double factor) const
    {
      return parallelity_proxy(
        i_seqs, j_seqs, sym_ops, weight*factor, target_angle_deg,
        slack, limit, top_out, origin_id);
    }

    //! Sorts i_seqs and j_seqs such that, e.g. i_seq[0] < i_seq[2].
    parallelity_proxy
    sort_ij_seqs() const
    {
      af::const_ref<std::size_t> i_seqs_cr = i_seqs.const_ref();
      af::const_ref<std::size_t> j_seqs_cr = j_seqs.const_ref();
      i_seqs_type i_seqs_result, j_seqs_result;
      i_seqs_result.reserve(i_seqs_cr.size());
      j_seqs_result.reserve(j_seqs_cr.size());
      i_seqs_type i_perm = af::sort_permutation(i_seqs_cr);
      i_seqs_type j_perm = af::sort_permutation(j_seqs_cr);
      af::const_ref<std::size_t> i_perm_cr = i_perm.const_ref();
      af::const_ref<std::size_t> j_perm_cr = j_perm.const_ref();
      for(std::size_t i=0;i<i_seqs_cr.size();i++) {
        i_seqs_result.push_back(i_seqs_cr[i_perm_cr[i]]);
      }
      for(std::size_t j=0;j<j_seqs_cr.size();j++) {
        j_seqs_result.push_back(j_seqs_cr[j_perm_cr[j]]);
      }

      // Need to figure out what's going on here with sym_ops
      // They are not used yet.
      if ( sym_ops.get() != 0 ) {
        af::const_ref<sgtbx::rt_mx> sym_ops_cr = sym_ops.get()->const_ref();
        af::shared<sgtbx::rt_mx> sym_ops_result;
        sym_ops_result.reserve(sym_ops_cr.size());
        for(std::size_t i=0;i<i_seqs_cr.size();i++) {
          sym_ops_result.push_back(sym_ops_cr[i]);
        }
        return parallelity_proxy(
          i_seqs_result,
          j_seqs_result,
          optional_container<af::shared<sgtbx::rt_mx> >(sym_ops_result),
          weight,
          target_angle_deg,
          slack,
          limit,
          top_out,
          origin_id);
      }
      else {
        return parallelity_proxy(
          i_seqs_result,j_seqs_result, weight, target_angle_deg, slack,
          limit, top_out, origin_id);
      }
    }

    //! Indices into array of sites.
    i_seqs_type i_seqs;
    i_seqs_type j_seqs;
    //! Array of symmetry operations.
    optional_container<af::shared<sgtbx::rt_mx> > sym_ops;
    //! Weight.
    double weight;
    //! Target angle in degrees
    double target_angle_deg;
    double slack;
    //! Parameter l in top-out potential
    double limit;
    //! Use top-out potential instead of harmonic
    bool top_out;
    unsigned char origin_id;
  };


  class parallelity
  {
    protected:
      scitbx::vec3<double> i_n, j_n;
      af::shared<scitbx::vec3<double> > i_dF__doriginal, j_dF__doriginal;

      // Function to calculate commonly used values for residual and
      // gradient calculations
      scitbx::vec3<double>
      calculate_C(af::const_ref<scitbx::vec3<double> > sites)
      {
        scitbx::vec3<double> result = 0;
        for(std::size_t i_site=0;i_site<sites.size();i_site++) {
          result += sites[i_site];
        }
        result /= double(sites.size());
        return result;
      }

      scitbx::mat3<double>
      calculate_S(
        af::const_ref<scitbx::vec3<double> > sites,
        scitbx::vec3<double> const& center_of_mass)
      {
        scitbx::mat3<double> result(0,0,0,0,0,0,0,0,0);
        for(std::size_t i_site=0;i_site<sites.size();i_site++) {
          scitbx::vec3<double> x = sites[i_site] - center_of_mass;
          result(0,0) += x[0]*x[0];
          result(1,1) += x[1]*x[1];
          result(2,2) += x[2]*x[2];
          result(0,1) += x[0]*x[1];
          result(0,2) += x[0]*x[2];
          result(1,2) += x[1]*x[2];
        }
        result(1,0) = result(0,1);
        result(2,0) = result(0,2);
        result(2,1) = result(1,2);
        return result;
      }

      double
      calculate_a_S(scitbx::mat3<double> const& S)
      {  return -S(0,0)-S(1,1)-S(2,2);}

      double
      calculate_b_S(scitbx::mat3<double> const& S)
      {
        return S(0,0)*S(1,1)-scitbx::fn::pow2(S(0,1)) +
               S(1,1)*S(2,2)-scitbx::fn::pow2(S(1,2)) +
               S(2,2)*S(0,0)-scitbx::fn::pow2(S(0,2));
      }

      double
      calculate_c_S(scitbx::mat3<double> const& S)
      {
        return -S.determinant();
      }

      double
      calculate_g_S(double a_S, double b_S)
      {
        return scitbx::fn::pow2(a_S)/3.0-b_S;
      }

      double
      calculate_h_S(double a_S, double b_S, double c_S)
      {
        return 2.0*scitbx::fn::pow3(a_S)/27.0-a_S*b_S/3.0+c_S;
      }

      double calculate_tau_S(double g_S, double h_S)
      {
        double e_S = -0.5*h_S*pow(1.0/3.0*g_S,-1.5);
        return cos(acos(e_S)/3.0+2.0/3.0*scitbx::constants::pi);
      }

      double
      calculate_lambda_plus(double a_S, double g_S, double tau_S)
      {
        return 2*sqrt(g_S/3.0)*tau_S-a_S/3.0;
      }

      int
      calculate_N(
        scitbx::mat3<double> const& S,
        double lambda_plus,
        scitbx::vec3<double> & t_1,
        scitbx::vec3<double> & t_2,
        scitbx::vec3<double> & N)
      {
        scitbx::vec3<double> N1, N2, N3;
        scitbx::vec3<double> t1(S(0,0)-lambda_plus, S(0,1), S(0,2)),
                             t2(S(1,0), S(1,1)-lambda_plus, S(1,2)),
                             t3(S(2,0), S(2,1), S(2,2)-lambda_plus);
        int result;
        N1 = t1.cross(t2);
        double nN1 = N1.length();
        N2 = t2.cross(t3);
        double nN2 = N2.length();
        double nN12;
        if (nN2 > nN1) {
          result = 2;
          nN12 = nN2;
          N = N2;
          t_1 = t2;
          t_2 = t3;
        }
        else {
          result = 1;
          N = N1;
          nN12 = nN1;
          t_1 = t1;
          t_2 = t2;
        }
        N3 = t3.cross(t1);
        double nN3 = N3.length();
        if (nN3 > nN12) {
          result = 3;
          N = N3;
          t_1 = t3;
          t_2 = t1;
        }
        return result;
      }

      void
      derive_dFparallelity__dn(
        scitbx::vec3<double> const& i_n,
        scitbx::vec3<double> const& j_n,
        double dFparallelity__dCtheta,
        double dFparallelity__dStheta,
        scitbx::vec3<double> & dFparallelity__dn_1,
        scitbx::vec3<double> & dFparallelity__dn_2)
      {
        scitbx::vec3<double> n = i_n.cross(j_n);
        double Stheta = n.length();
        if (fabs(Stheta) > 1.e-100) {

          dFparallelity__dn_1[0] = j_n[0]*dFparallelity__dCtheta +
            (j_n[1]*n[2]/Stheta - j_n[2]*n[1]/Stheta) * dFparallelity__dStheta;

          dFparallelity__dn_1[1] = j_n[1]*dFparallelity__dCtheta +
            (j_n[2]*n[0]/Stheta - j_n[0]*n[2]/Stheta) * dFparallelity__dStheta;

          dFparallelity__dn_1[2] = j_n[2]*dFparallelity__dCtheta +
            (j_n[0]*n[1]/Stheta - j_n[1]*n[0]/Stheta) * dFparallelity__dStheta;

          dFparallelity__dn_2[0] = i_n[0]*dFparallelity__dCtheta -
            (i_n[1]*n[2]/Stheta - i_n[2]*n[1]/Stheta) * dFparallelity__dStheta;

          dFparallelity__dn_2[1] = i_n[1]*dFparallelity__dCtheta -
            (i_n[2]*n[0]/Stheta - i_n[0]*n[2]/Stheta) * dFparallelity__dStheta;

          dFparallelity__dn_2[2] = i_n[2]*dFparallelity__dCtheta -
            (i_n[0]*n[1]/Stheta - i_n[1]*n[0]/Stheta) * dFparallelity__dStheta;
        }
      }

      scitbx::vec3<double>
      derive_dFparallelity__dN(
        scitbx::vec3<double> const& dFparallelity__dn,
        scitbx::vec3<double> const& N)
      {
        scitbx::vec3<double> result;
        double before_brackets = 1.0/pow((N*N),1.5);
        result[0] = before_brackets*((N[1]*N[1]+N[2]*N[2])*dFparallelity__dn[0]-
                                      N[0]*N[1]*dFparallelity__dn[1]-
                                      N[0]*N[2]*dFparallelity__dn[2]);
        result[1] = before_brackets*((N[2]*N[2]+N[0]*N[0])*dFparallelity__dn[1]-
                                      N[1]*N[2]*dFparallelity__dn[2]-
                                      N[1]*N[0]*dFparallelity__dn[0]);
        result[2] = before_brackets*((N[0]*N[0]+N[1]*N[1])*dFparallelity__dn[2]-
                                      N[2]*N[0]*dFparallelity__dn[0]-
                                      N[2]*N[1]*dFparallelity__dn[1]);
        return result;
      }

      void
      derive_dFparallelity__dt(
        scitbx::vec3<double> & dF__dt_1,
        scitbx::vec3<double> & dF__dt_2,
        scitbx::vec3<double> const& dF__dN,
        scitbx::vec3<double> const& t_1,
        scitbx::vec3<double> const& t_2)
      {
        dF__dt_1[0] = -t_2[2]*dF__dN[1] + t_2[1]*dF__dN[2];
        dF__dt_1[1] = -t_2[0]*dF__dN[2] + t_2[2]*dF__dN[0];
        dF__dt_1[2] = -t_2[1]*dF__dN[0] + t_2[0]*dF__dN[1];
        dF__dt_2[0] = t_1[2]*dF__dN[1] - t_1[1]*dF__dN[2];
        dF__dt_2[1] = t_1[0]*dF__dN[2] - t_1[2]*dF__dN[0];
        dF__dt_2[2] = t_1[1]*dF__dN[0] - t_1[0]*dF__dN[1];
      }

      double
      derive_dF__dS__dlambda_plus(
        scitbx::mat3<double> & dF__dS,
        scitbx::vec3<double> const& dF__dt_1,
        scitbx::vec3<double> const& dF__dt_2,
        int variant)
      {
        if (variant == 1) {
          dF__dS(0,0) = dF__dt_1[0];
          dF__dS(0,1) = dF__dt_1[1];
          dF__dS(0,2) = dF__dt_1[2];
          dF__dS(1,0) = dF__dt_2[0];
          dF__dS(1,1) = dF__dt_2[1];
          dF__dS(1,2) = dF__dt_2[2];
          dF__dS(2,0) = dF__dS(2,1) = dF__dS(2,2) = 0.0;
          return -dF__dt_1[0] - dF__dt_2[1];
        }
        if (variant == 2) {
          dF__dS(1,0) = dF__dt_1[0];
          dF__dS(1,1) = dF__dt_1[1];
          dF__dS(1,2) = dF__dt_1[2];
          dF__dS(2,0) = dF__dt_2[0];
          dF__dS(2,1) = dF__dt_2[1];
          dF__dS(2,2) = dF__dt_2[2];
          dF__dS(0,0) = 0.0; dF__dS(0,1) = 0.0; dF__dS(0,2) = 0.0;
          return -dF__dt_1[1] - dF__dt_2[2];
        }
        if (variant == 3) {
          dF__dS(2,0) = dF__dt_1[0];
          dF__dS(2,1) = dF__dt_1[1];
          dF__dS(2,2) = dF__dt_1[2];
          dF__dS(0,0) = dF__dt_2[0];
          dF__dS(0,1) = dF__dt_2[1];
          dF__dS(0,2) = dF__dt_2[2];
          dF__dS(1,0) = 0.0; dF__dS(1,1) = 0.0; dF__dS(1,2) = 0.0;
          return -dF__dt_1[2] - dF__dt_2[0];
        }
        std::cout << "Variant number:" << variant << "\n";
        CCTBX_ASSERT(1 == 2); // Should never get here
      }

      double det(double a00, double a01, double a10, double a11)
      {
        scitbx::mat2<double> m(a00, a01, a10, a11);
        return m.determinant();
      }
      scitbx::mat3<double> derive_da__dS()
      {
        scitbx::mat3<double> result;
        result.fill(0.0);
        result(0,0) = result(1,1) = result(2,2) = -1;
        return result;
      }

      scitbx::mat3<double> derive_db__dS(
        scitbx::mat3<double> const& S)
      {
        scitbx::mat3<double> result;
        result(0,0) = S(1,1)+S(2,2);
        result(1,1) = S(2,2)+S(0,0);
        result(2,2) = S(0,0)+S(1,1);
        result(0,1) = -S(1,0);
        result(1,0) = -S(0,1);
        result(1,2) = -S(2,1);
        result(2,1) = -S(1,2);
        result(2,0) = -S(0,2);
        result(0,2) = -S(2,0);
        return result;
      }

      scitbx::mat3<double> derive_dc__dS(
        scitbx::mat3<double> const& S)
      {
        scitbx::mat3<double> result;
        result(0,0) = -det(S(1,1), S(1,2), S(2,1), S(2,2));
        result(0,1) = -det(S(1,2), S(1,0), S(2,2), S(2,0));
        result(0,2) = -det(S(1,0), S(1,1), S(2,0), S(2,1));
        result(1,0) = -det(S(2,1), S(2,2), S(0,1), S(0,2));
        result(1,1) = -det(S(2,2), S(2,0), S(0,2), S(0,0));
        result(1,2) = -det(S(2,0), S(2,1), S(0,0), S(0,1));
        result(2,0) = -det(S(0,1), S(0,2), S(1,1), S(1,2));
        result(2,1) = -det(S(0,2), S(0,0), S(1,2), S(1,0));
        result(2,2) = -det(S(0,0), S(0,1), S(1,0), S(1,1));
        return result;
      }


      scitbx::mat3<double>
      derive_DF__DS(
        double DF__Da_S,
        double DF__Db_S,
        double DF__Dc_S,
        scitbx::mat3<double> const& dF__dS,
        scitbx::mat3<double> const& S)
      {
        scitbx::mat3<double> result, da__dS, db__dS, dc__dS;
        da__dS = derive_da__dS();
        db__dS = derive_db__dS(S);
        dc__dS = derive_dc__dS(S);
        for (int i=0; i<3; i++) {
          for (int j=0; j<3; j++) {
            result(i,j) = da__dS(i,j)*DF__Da_S +
                          db__dS(i,j)*DF__Db_S +
                          dc__dS(i,j)*DF__Dc_S +
                          dF__dS(i,j);
          }
        }
        return result;
      }

      void
      derive_dF__dcentered(
        af::shared<scitbx::vec3<double> > &       dF__dcentered,
        scitbx::mat3<double> const&               DF__DS,
        af::const_ref<scitbx::vec3<double> >      sites,
        scitbx::vec3<double> const&               center_of_mass)
      {
        for(std::size_t i_site=0;i_site<sites.size();i_site++) {
          scitbx::vec3<double> x = sites[i_site] - center_of_mass;
          scitbx::vec3<double> dF__dX_K;
          dF__dX_K[0] = 2.0*x[0]* DF__DS(0,0) +
                            x[1]*(DF__DS(0,1)+DF__DS(1,0)) +
                            x[2]*(DF__DS(0,2)+DF__DS(2,0));
          dF__dX_K[1] = 2.0*x[1]* DF__DS(1,1) +
                            x[2]*(DF__DS(1,2)+DF__DS(2,1)) +
                            x[0]*(DF__DS(1,0)+DF__DS(0,1));
          dF__dX_K[2] = 2.0*x[2]* DF__DS(2,2) +
                            x[0]*(DF__DS(2,0)+DF__DS(0,2)) +
                            x[1]*(DF__DS(2,1)+DF__DS(1,2));
          dF__dcentered.push_back(dF__dX_K);
        }
      }

      void
      derive_dF__doriginal(
        af::shared<scitbx::vec3<double> > &       dF__doriginal,
        af::shared<scitbx::vec3<double> > const & dF__dcentered)
      {
        scitbx::vec3<double> sum_dF__dY_divK = 0;
        for (std::size_t i_site=0;i_site<dF__dcentered.size();i_site++) {
          sum_dF__dY_divK += dF__dcentered[i_site];
        }
        sum_dF__dY_divK /= double(dF__dcentered.size());
        for(std::size_t i_site=0;i_site<dF__dcentered.size();i_site++) {
          scitbx::vec3<double> dF__dx;
          dF__dx =  dF__dcentered[i_site] - sum_dF__dY_divK;
          dF__doriginal.push_back(dF__dx);
        }
      }

      void
      init_deltas()
      {
        af::const_ref<scitbx::vec3<double> > i_sites_ref = i_sites.const_ref();
        af::const_ref<scitbx::vec3<double> > j_sites_ref = j_sites.const_ref();
        scitbx::vec3<double> i_center_of_mass_ = calculate_C(i_sites_ref);
        scitbx::vec3<double> j_center_of_mass_ = calculate_C(j_sites_ref);
        scitbx::mat3<double> i_S_ = calculate_S(i_sites_ref, i_center_of_mass_);
        scitbx::mat3<double> j_S_ = calculate_S(j_sites_ref, j_center_of_mass_);
        double i_a_S = calculate_a_S(i_S_);
        double i_b_S = calculate_b_S(i_S_);
        double i_c_S = calculate_c_S(i_S_);
        double j_a_S = calculate_a_S(j_S_);
        double j_b_S = calculate_b_S(j_S_);
        double j_c_S = calculate_c_S(j_S_);
        double i_g_S = calculate_g_S(i_a_S, i_b_S);
        double i_h_S = calculate_h_S(i_a_S, i_b_S, i_c_S);
        double j_g_S = calculate_g_S(j_a_S, j_b_S);
        double j_h_S = calculate_h_S(j_a_S, j_b_S, j_c_S);
        double i_tau_S = calculate_tau_S(i_g_S, i_h_S);
        double j_tau_S = calculate_tau_S(j_g_S, j_h_S);
        double i_lambda_plus = calculate_lambda_plus(i_a_S, i_g_S, i_tau_S);
        double j_lambda_plus = calculate_lambda_plus(j_a_S, j_g_S, j_tau_S);
        scitbx::vec3<double> i_t_1, i_t_2, i_N;
        int i_variant = calculate_N(i_S_, i_lambda_plus, i_t_1, i_t_2, i_N);
        scitbx::vec3<double> j_t_1, j_t_2, j_N;
        int j_variant = calculate_N(j_S_, j_lambda_plus, j_t_1, j_t_2, j_N);
        if (i_N * j_N < 0) {
          j_N = -j_N;
        }
        i_n = i_N.normalize();
        j_n = j_N.normalize();

        double n1n2 = i_n * j_n;
        double theta = acos(n1n2);
        theta_deg = scitbx::rad_as_deg(theta);
        delta = theta_deg - target_angle_deg;
        theta1=target_angle_deg;

        if (fabs(delta) <= slack) {
          delta_slack = 0;
        }
        else if (delta > slack) {
          delta_slack = delta - slack;
          theta1 = target_angle_deg+slack;
        }
        else if (delta < slack) {
          delta_slack = delta+slack;
          theta1 = target_angle_deg-slack;
        }
        i_dF__doriginal.reserve(i_sites_ref.size());
        j_dF__doriginal.reserve(j_sites_ref.size());
        if (fabs(delta_slack) < 1.e-100) {
          // No need to calculate gradients
          scitbx::vec3<double> zero_vector(0,0,0);
          for(std::size_t i_site=0; i_site<i_sites_ref.size();i_site++) {
            i_dF__doriginal.push_back(zero_vector);
          }
          for(std::size_t i_site=0; i_site<j_sites_ref.size();i_site++) {
            j_dF__doriginal.push_back(zero_vector);
          }
        }
        else {

          //==================================
          // gradients
          //==================================
          double dF__dCtheta, dF__dStheta;
          scitbx::vec3<double> dFparallelity__dn_1(0,0,0),
                               dFparallelity__dn_2(0,0,0);

          if (top_out) {
            double l2 = limit * limit;
            double dF__dCDtheta = -weight*std::exp(
                (cos(scitbx::deg_as_rad(theta_deg - theta1))-1)/l2);
            dF__dCtheta = dF__dCDtheta*cos(scitbx::deg_as_rad(theta1));
            dF__dStheta = dF__dCDtheta*sin(scitbx::deg_as_rad(theta1));
          }
          else {
            dF__dCtheta = -weight*cos(scitbx::deg_as_rad(theta1));
            dF__dStheta = -weight*sin(scitbx::deg_as_rad(theta1));
          }
          derive_dFparallelity__dn(i_n, j_n, dF__dCtheta, dF__dStheta,
              dFparallelity__dn_1, dFparallelity__dn_2);

          scitbx::vec3<double> dF__dN_1, dF__dN_2;
          dF__dN_1 = derive_dFparallelity__dN(dFparallelity__dn_1, i_N);
          dF__dN_2 = derive_dFparallelity__dN(dFparallelity__dn_2, j_N);
          scitbx::vec3<double> i_dF__dt_1, i_dF__dt_2, j_dF__dt_1, j_dF__dt_2;
          derive_dFparallelity__dt(
              i_dF__dt_1, i_dF__dt_2, dF__dN_1, i_t_1, i_t_2);
          derive_dFparallelity__dt(
              j_dF__dt_1, j_dF__dt_2, dF__dN_2, j_t_1, j_t_2);
          scitbx::mat3<double> i_dF__dS, j_dF__dS;
          double i_dF__dlambda_plus, j_dF__dlambda_plus;
          i_dF__dlambda_plus = derive_dF__dS__dlambda_plus(
                                    i_dF__dS, i_dF__dt_1, i_dF__dt_2, i_variant);
          j_dF__dlambda_plus = derive_dF__dS__dlambda_plus(
                                    j_dF__dS, j_dF__dt_1, j_dF__dt_2, j_variant);
          double i_dF__da_S = -i_dF__dlambda_plus/3.0;
          double j_dF__da_S = -j_dF__dlambda_plus/3.0;
          double i_dF__dg_S = i_dF__dlambda_plus*pow(3*i_g_S,-0.5)*i_tau_S;
          double j_dF__dg_S = j_dF__dlambda_plus*pow(3*j_g_S,-0.5)*j_tau_S;
          double i_dF__dtau_S = i_dF__dlambda_plus*2*pow(i_g_S/3.0,0.5);
          double j_dF__dtau_S = j_dF__dlambda_plus*2*pow(j_g_S/3.0,0.5);
          double i_dF_dksi_S = i_dF__dtau_S/(12.0*i_tau_S*i_tau_S-3.0);
          double j_dF_dksi_S = j_dF__dtau_S/(12.0*j_tau_S*j_tau_S-3.0);
          //p.21
          double i_DF__Dg_S = i_dF_dksi_S*i_h_S/4.0*pow(i_g_S/3.0,-2.5)+
                              i_dF__dg_S;
          double j_DF__Dg_S = j_dF_dksi_S*j_h_S/4.0*pow(j_g_S/3.0,-2.5)+
                              j_dF__dg_S;
          double i_dF__dh_S = -i_dF_dksi_S/2.0*pow(i_g_S/3.0, -1.5);
          double j_dF__dh_S = -j_dF_dksi_S/2.0*pow(j_g_S/3.0, -1.5);
          //p.22
          double i_DF__Dc_S =  i_dF__dh_S;
          double j_DF__Dc_S =  j_dF__dh_S;
          double i_DF__Db_S = -i_dF__dh_S/3.0*i_a_S - i_DF__Dg_S;
          double j_DF__Db_S = -j_dF__dh_S/3.0*j_a_S - j_DF__Dg_S;
          double i_DF__Da_S = i_dF__dh_S*(2.0/9.0*i_a_S*i_a_S-i_b_S/3.0) +
                                  2.0/3.0*i_a_S*i_DF__Dg_S + i_dF__da_S;
          double j_DF__Da_S = j_dF__dh_S*(2.0/9.0*j_a_S*j_a_S-j_b_S/3.0) +
                                  2.0/3.0*j_a_S*j_DF__Dg_S + j_dF__da_S;
          //p.23
          scitbx::mat3<double> i_DF__DS, j_DF__DS;
          i_DF__DS = derive_DF__DS(
              i_DF__Da_S, i_DF__Db_S, i_DF__Dc_S, i_dF__dS, i_S_);
          j_DF__DS = derive_DF__DS(
              j_DF__Da_S, j_DF__Db_S, j_DF__Dc_S, j_dF__dS, j_S_);
          af::shared<scitbx::vec3<double> > i_dF__dcentered, j_dF__dcentered;
          i_dF__dcentered.reserve(i_sites_ref.size());
          j_dF__dcentered.reserve(j_sites_ref.size());
          derive_dF__dcentered(
              i_dF__dcentered, i_DF__DS, i_sites_ref, i_center_of_mass_);
          derive_dF__dcentered(
              j_dF__dcentered, j_DF__DS, j_sites_ref, j_center_of_mass_);
          derive_dF__doriginal(i_dF__doriginal, i_dF__dcentered);
          derive_dF__doriginal(j_dF__doriginal, j_dF__dcentered);
        }
      }

    public:
      //! Cartesian coordinates of the sites defining the first (i_sites) and
      //! second (j_sites) plane.
      af::shared<scitbx::vec3<double> > i_sites, j_sites;
      //! Weight.
      double weight;
      //! Target angle in degrees
      double target_angle_deg;
      //! Difference target_angle - model_angle
      double delta;
      //! sign(delta) * max(0, (abs(delta) - slack))
      double delta_slack;
      double theta1, theta_deg;
      double slack;
      double limit;
      bool top_out;

      //! Default constructor. Some data members are not initialized!
      parallelity() {}

      //! Constructor.
      parallelity(
        af::shared<scitbx::vec3<double> > const& i_sites_,
        af::shared<scitbx::vec3<double> > const& j_sites_,
        double weight_,
        double target_angle_deg_,
        double slack_ = 0,
        double limit_ = 1,
        bool top_out_ = false)
      :
        i_sites(i_sites_),
        j_sites(j_sites_),
        weight(weight_),
        target_angle_deg(target_angle_deg_),
        slack(slack_),
        limit(limit_),
        top_out(top_out_)
      {
        CCTBX_ASSERT(i_sites.size() > 2);
        CCTBX_ASSERT(j_sites.size() > 2);
        CCTBX_ASSERT(limit >= 1);
        init_deltas();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs and proxy.j_seqs, weights are copied from proxy.
       */
      parallelity(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        parallelity_proxy const& proxy)
      :
        weight(proxy.weight),
        target_angle_deg(proxy.target_angle_deg),
        slack(proxy.slack),
        limit(proxy.limit),
        top_out(proxy.top_out)
      {
        af::const_ref<std::size_t> i_seqs_ref = proxy.i_seqs.const_ref();
        af::const_ref<std::size_t> j_seqs_ref = proxy.j_seqs.const_ref();
        i_sites.reserve(i_seqs_ref.size());
        j_sites.reserve(j_seqs_ref.size());
        for(std::size_t i=0;i<i_seqs_ref.size();i++) {
          std::size_t i_seq = i_seqs_ref[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          i_sites.push_back(sites_cart[i_seq]);
        }
        for(std::size_t i=0;i<j_seqs_ref.size();i++) {
          std::size_t j_seq = j_seqs_ref[i];
          CCTBX_ASSERT(j_seq < sites_cart.size());
          j_sites.push_back(sites_cart[j_seq]);
        }
        init_deltas();
      }

      /*! \brief Coordinates are obtained from sites_cart according
          to proxy.i_seqs by applying proxy.sym_ops and unit_cell,
          weights are copied from proxy.
       */
      parallelity(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        parallelity_proxy const& proxy)
      :
        weight(proxy.weight),
        target_angle_deg(proxy.target_angle_deg),
        slack(proxy.slack),
        limit(proxy.limit),
        top_out(proxy.top_out)
      {
        af::const_ref<std::size_t> i_seqs_ref = proxy.i_seqs.const_ref();
        af::const_ref<std::size_t> j_seqs_ref = proxy.j_seqs.const_ref();
        i_sites.reserve(i_seqs_ref.size());
        j_sites.reserve(j_seqs_ref.size());
        for(std::size_t i=0;i<i_seqs_ref.size();i++) {
          std::size_t i_seq = i_seqs_ref[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          i_sites.push_back(sites_cart[i_seq]);
          if ( proxy.sym_ops.get() != 0 ) {
            sgtbx::rt_mx rt_mx = proxy.sym_ops[i];
            if ( !rt_mx.is_unit_mx() ) {
              i_sites[i] = unit_cell.orthogonalize(
                rt_mx * unit_cell.fractionalize(i_sites[i]));
            }
          }
        }
        for(std::size_t i=0;i<j_seqs_ref.size();i++) {
          std::size_t j_seq = j_seqs_ref[i];
          CCTBX_ASSERT(j_seq < sites_cart.size());
          j_sites.push_back(sites_cart[j_seq]);
          if ( proxy.sym_ops.get() != 0 ) {
            sgtbx::rt_mx rt_mx = proxy.sym_ops[i];
            if ( !rt_mx.is_unit_mx() ) {
              j_sites[i+i_seqs_ref.size()] = unit_cell.orthogonalize(
                rt_mx * unit_cell.fractionalize(j_sites[i+i_seqs_ref.size()]));
            }
          }
        }
        init_deltas();
      }

      //! Sum of weight * delta**2 over all sites.
      double
      residual() const
      {
        if (fabs(delta_slack) < 1.e-100) {return 0;}
        else {
          if (top_out) {
            double l2 = limit * limit;
            return weight*l2*(1-std::exp(
                (cos(scitbx::deg_as_rad(theta_deg - theta1))-1)/l2));
          }
          else
            return weight*(1-cos(scitbx::deg_as_rad(theta_deg - theta1)));
        }
      }

      //! Gradients with respect to the sites.
      af::shared<scitbx::vec3<double> >
      gradients() const
      {
        af::shared<scitbx::vec3<double> > i_result, j_result, result;
        result.reserve(i_sites.size()+j_sites.size());
        for(std::size_t i_site=0;i_site<i_sites.size();i_site++) {
          result.push_back(i_dF__doriginal[i_site]);
        }
        for(std::size_t j_site=0;j_site<j_sites.size();j_site++) {
          result.push_back(j_dF__doriginal[j_site]);
        }
        return result;
      }

      //! Support for planarity_residual_sum ????
      /*! Not available in Python.
       */
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        parallelity_proxy::i_seqs_type const& i_seqs,
        parallelity_proxy::i_seqs_type const& j_seqs) const
      {
        std::size_t i_size = i_seqs.size();
        af::const_ref<std::size_t> i_seqs_ref = i_seqs.const_ref();
        af::const_ref<std::size_t> j_seqs_ref = j_seqs.const_ref();
        af::shared<scitbx::vec3<double> > grads = gradients();
        af::const_ref<scitbx::vec3<double> > grads_ref = grads.const_ref();
        for(std::size_t i=0;i<i_size;i++) {
          gradient_array[i_seqs_ref[i]] += grads_ref[i];
        }
        for (std::size_t j=i_size; j<j_seqs.size(); j++) {
          gradient_array[j_seqs_ref[j-i_size]] += grads_ref[j];
        }
      }

      //! Support for planarity_residual_sum. ????
      /*! Not available in Python.

          Inefficient implementation, r_inv_cart is not cached.
          TODO: use asu_mappings to take advantage of caching of r_inv_cart.
       */
      void
      add_gradients(
        uctbx::unit_cell const& unit_cell,
        af::ref<scitbx::vec3<double> > const& gradient_array,
        parallelity_proxy const& proxy) const
      {
        CCTBX_ASSERT(1 == 2);
      }
  };

  inline
  af::shared<double>
  parallelity_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<parallelity_proxy> const& proxies)
  {
    return detail::generic_deltas<parallelity_proxy, parallelity>::get(
      sites_cart, proxies);
  }

  /*! \brief Fast computation of parallelity::residual() given an array
      of parallelity proxies, ignoring proxy.sym_ops.
   */
  inline
  af::shared<double>
  parallelity_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<parallelity_proxy> const& proxies)
  {
    return detail::generic_residuals<parallelity_proxy, parallelity>::get(
      sites_cart, proxies);
  }

  inline
  double
  parallelity_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<parallelity_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    // due to presence j_seqs we cannot directly use this functions, so...
    //return detail::generic_residual_sum<parallelity_proxy, parallelity>::get(
    //  sites_cart, proxies, gradient_array);
    CCTBX_ASSERT(   gradient_array.size() == 0
                 || gradient_array.size() == sites_cart.size());
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      parallelity_proxy const& proxy = proxies[i];
      parallelity restraint(sites_cart, proxy);
      result += restraint.residual();
      if (gradient_array.size() != 0) {
        restraint.add_gradients(gradient_array, proxy.i_seqs, proxy.j_seqs);
      }
    }
    return result;
  }

  /*! \brief Fast computation of parallelity::rms_deltas() given an array
      of parallelity sym proxies, taking account of proxy.sym_ops.
   */
  inline
  af::shared<double>
  parallelity_deltas(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<parallelity_proxy> const& proxies)
  {
    CCTBX_ASSERT(1 == 2);
    return detail::generic_deltas<parallelity_proxy, parallelity>::get(
      sites_cart, proxies);
  }

  /*! \brief Fast computation of parallelity::residual() given an array
      of parallelity sym proxies, taking account of proxy.sym_ops.
   */
  inline
  af::shared<double>
  parallelity_residuals(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<parallelity_proxy> const& proxies)
  {
    return detail::generic_residuals<parallelity_proxy, parallelity>::get(
      unit_cell, sites_cart, proxies);
  }

  inline
  double
  parallelity_residual_sum(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<parallelity_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<parallelity_proxy, parallelity>::get(
      unit_cell, sites_cart, proxies, gradient_array);
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_PARALLELITY_H
