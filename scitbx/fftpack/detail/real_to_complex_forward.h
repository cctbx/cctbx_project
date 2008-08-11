#ifndef SCITBX_FFTPACK_DETAIL_REAL_TO_COMPLEX_FORWARD_H
#define SCITBX_FFTPACK_DETAIL_REAL_TO_COMPLEX_FORWARD_H

namespace scitbx { namespace fftpack {

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType, ComplexType>::forward_compressed(
    real_type* c, /* seq */
    real_type* ch /* scratch */)
  {
    if (n_ < 2) return;
    const real_type* wa = &(*(wa_.begin()));
    std::size_t idl1;
    std::size_t ido;
    std::size_t ip;
    std::size_t iw;
    std::size_t ix2;
    std::size_t ix3;
    std::size_t ix4;
    std::size_t kh;
    std::size_t l1;
    std::size_t l2;
    std::size_t na;
    na = 1;
    l2 = n_;
    iw = n_;
    for (std::size_t k1 = 1; k1 <= factors_.size(); k1++) {
      kh = factors_.size()-k1;
      ip = factors_[kh+3-2-1];
      l1 = l2/ip;
      ido = n_/l2;
      idl1 = ido*l1;
      iw = iw-(ip-1)*ido;
      na = 1-na;
      if (ip == 4) {
        ix2 = iw+ido;
        ix3 = ix2+ido;
        if (na == 0) {
          passf4(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1);
        }
        else {
          passf4(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1);
        }
      }
      else if (ip == 2) {
        if (na == 0) {
          passf2(ido,l1,c,ch,wa+iw-1);
        }
        else {
          passf2(ido,l1,ch,c,wa+iw-1);
        }
      }
      else if (ip == 3) {
        ix2 = iw+ido;
        if (na == 0) {
          passf3(ido,l1,c,ch,wa+iw-1,wa+ix2-1);
        }
        else {
          passf3(ido,l1,ch,c,wa+iw-1,wa+ix2-1);
        }
      }
      else if (ip == 5) {
        ix2 = iw+ido;
        ix3 = ix2+ido;
        ix4 = ix3+ido;
        if (na == 0) {
          passf5(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1,wa+ix4-1);
        }
        else {
          passf5(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1,wa+ix4-1);
        }
      }
      else {
        if (ido == 1) na = 1-na;
        if (na == 0) {
          passfg(ido,ip,l1,idl1,c,c,c,ch,ch,wa+iw-1);
          na = 1;
        }
        else {
          passfg(ido,ip,l1,idl1,ch,ch,ch,c,c,wa+iw-1);
          na = 0;
        }
      }
      l2 = l1;
    }
    if (na == 1) return;
    for (std::size_t i = 1; i <= n_; i++) {
      c[i-1] = ch[i-1];
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passf2(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1)
  {
    dim3 cc(cc_begin, ido, l1, 2);
    dim3 ch(ch_begin, ido, 2, l1);
    std::size_t ic;
    real_type ti2;
    real_type tr2;
    std::size_t k;
    for (k = 0; k < l1; k++) {
      ch(0,0,k) = cc(0,k,0)+cc(0,k,1);
      ch(ido-1,1,k) = cc(0,k,0)-cc(0,k,1);
    }
    if (ido < 2) return;
    if (ido > 2) {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
          std::size_t i1 = i0 + 1;
          ic = ido-i1;
          tr2 = wa1[i0-1]*cc(i0,k,1)+wa1[i0]*cc(i1,k,1);
          ti2 = wa1[i0-1]*cc(i1,k,1)-wa1[i0]*cc(i0,k,1);
          ch(i1,0,k) = cc(i1,k,0)+ti2;
          ch(ic,1,k) = ti2-cc(i1,k,0);
          ch(i0,0,k) = cc(i0,k,0)+tr2;
          ch(ic-1,1,k) = cc(i0,k,0)-tr2;
        }
      }
      if (ido % 2 != 0) return;
    }
    for (k = 0; k < l1; k++) {
      ch(0,1,k) = -cc(ido-1,k,1);
      ch(ido-1,0,k) = cc(ido-1,k,0);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passf3(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1,
                                       const real_type* wa2)
  {
    dim3 cc(cc_begin, ido, l1, 3);
    dim3 ch(ch_begin, ido, 3, l1);
    real_type ci2;
    real_type cr2;
    real_type di2;
    real_type di3;
    real_type dr2;
    real_type dr3;
    std::size_t ic;
    real_type ti2;
    real_type ti3;
    real_type tr2;
    real_type tr3;
    const real_type taur = -.5;
    const real_type taui = -taur * std::sqrt(real_type(3));
    std::size_t k;
    for (k = 0; k < l1; k++) {
      cr2 = cc(0,k,1)+cc(0,k,2);
      ch(0,0,k) = cc(0,k,0)+cr2;
      ch(0,2,k) = taui*(cc(0,k,2)-cc(0,k,1));
      ch(ido-1,1,k) = cc(0,k,0)+taur*cr2;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; k++) {
      for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
        std::size_t i1 = i0 + 1;
        ic = ido-i1;
        dr2 = wa1[i0-1]*cc(i0,k,1)+wa1[i0]*cc(i1,k,1);
        di2 = wa1[i0-1]*cc(i1,k,1)-wa1[i0]*cc(i0,k,1);
        dr3 = wa2[i0-1]*cc(i0,k,2)+wa2[i0]*cc(i1,k,2);
        di3 = wa2[i0-1]*cc(i1,k,2)-wa2[i0]*cc(i0,k,2);
        cr2 = dr2+dr3;
        ci2 = di2+di3;
        ch(i0,0,k) = cc(i0,k,0)+cr2;
        ch(i1,0,k) = cc(i1,k,0)+ci2;
        tr2 = cc(i0,k,0)+taur*cr2;
        ti2 = cc(i1,k,0)+taur*ci2;
        tr3 = taui*(di2-di3);
        ti3 = taui*(dr3-dr2);
        ch(i0,2,k) = tr2+tr3;
        ch(ic-1,1,k) = tr2-tr3;
        ch(i1,2,k) = ti2+ti3;
        ch(ic,1,k) = ti3-ti2;
      }
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passf4(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1,
                                       const real_type* wa2,
                                       const real_type* wa3)
  {
    dim3 cc(cc_begin, ido, l1, 4);
    dim3 ch(ch_begin, ido, 4, l1);
    real_type ci2;
    real_type ci3;
    real_type ci4;
    real_type cr2;
    real_type cr3;
    real_type cr4;
    std::size_t ic;
    real_type ti1;
    real_type ti2;
    real_type ti3;
    real_type ti4;
    real_type tr1;
    real_type tr2;
    real_type tr3;
    real_type tr4;
    const real_type hsqt2 = real_type(.5) * std::sqrt(real_type(2));
    std::size_t k;
    for (k = 0; k < l1; k++) {
      tr1 = cc(0,k,1)+cc(0,k,3);
      tr2 = cc(0,k,0)+cc(0,k,2);
      ch(0,0,k) = tr1+tr2;
      ch(ido-1,3,k) = tr2-tr1;
      ch(ido-1,1,k) = cc(0,k,0)-cc(0,k,2);
      ch(0,2,k) = cc(0,k,3)-cc(0,k,1);
    }
    if (ido < 2) return;
    if (ido > 2) {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
          std::size_t i1 = i0 + 1;
          ic = ido-i1;
          cr2 = wa1[i0-1]*cc(i0,k,1)+wa1[i0]*cc(i1,k,1);
          ci2 = wa1[i0-1]*cc(i1,k,1)-wa1[i0]*cc(i0,k,1);
          cr3 = wa2[i0-1]*cc(i0,k,2)+wa2[i0]*cc(i1,k,2);
          ci3 = wa2[i0-1]*cc(i1,k,2)-wa2[i0]*cc(i0,k,2);
          cr4 = wa3[i0-1]*cc(i0,k,3)+wa3[i0]*cc(i1,k,3);
          ci4 = wa3[i0-1]*cc(i1,k,3)-wa3[i0]*cc(i0,k,3);
          tr1 = cr2+cr4;
          tr4 = cr4-cr2;
          ti1 = ci2+ci4;
          ti4 = ci2-ci4;
          ti2 = cc(i1,k,0)+ci3;
          ti3 = cc(i1,k,0)-ci3;
          tr2 = cc(i0,k,0)+cr3;
          tr3 = cc(i0,k,0)-cr3;
          ch(i0,0,k) = tr1+tr2;
          ch(ic-1,3,k) = tr2-tr1;
          ch(i1,0,k) = ti1+ti2;
          ch(ic,3,k) = ti1-ti2;
          ch(i0,2,k) = ti4+tr3;
          ch(ic-1,1,k) = tr3-ti4;
          ch(i1,2,k) = tr4+ti3;
          ch(ic,1,k) = tr4-ti3;
        }
      }
      if (ido % 2 != 0) return;
    }
    for (k = 0; k < l1; k++) {
      ti1 = -hsqt2*(cc(ido-1,k,1)+cc(ido-1,k,3));
      tr1 = hsqt2*(cc(ido-1,k,1)-cc(ido-1,k,3));
      ch(ido-1,0,k) = tr1+cc(ido-1,k,0);
      ch(ido-1,2,k) = cc(ido-1,k,0)-tr1;
      ch(0,1,k) = ti1-cc(ido-1,k,2);
      ch(0,3,k) = ti1+cc(ido-1,k,2);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passf5(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1,
                                       const real_type* wa2,
                                       const real_type* wa3,
                                       const real_type* wa4)
  {
    dim3 cc(cc_begin, ido, l1, 5);
    dim3 ch(ch_begin, ido, 5, l1);
    real_type ci2;
    real_type ci3;
    real_type ci4;
    real_type ci5;
    real_type cr2;
    real_type cr3;
    real_type cr4;
    real_type cr5;
    real_type di2;
    real_type di3;
    real_type di4;
    real_type di5;
    real_type dr2;
    real_type dr3;
    real_type dr4;
    real_type dr5;
    std::size_t ic;
    real_type ti2;
    real_type ti3;
    real_type ti4;
    real_type ti5;
    real_type tr2;
    real_type tr3;
    real_type tr4;
    real_type tr5;
    // sin(18 deg)
    const real_type tr11 =  .309016994374947424102293417182819058860154589903;
    // cos(18 deg)
    const real_type ti11 =  .951056516295153572116439333379382143405698634126;
    // -cos(36 deg)
    const real_type tr12 = -.809016994374947424102293417182819058860154589903;
    // sin(36 deg)
    const real_type ti12 =  .587785252292473129168705954639072768597652437643;
    std::size_t k;
    for (k = 0; k < l1; k++) {
      cr2 = cc(0,k,4)+cc(0,k,1);
      ci5 = cc(0,k,4)-cc(0,k,1);
      cr3 = cc(0,k,3)+cc(0,k,2);
      ci4 = cc(0,k,3)-cc(0,k,2);
      ch(0,0,k) = cc(0,k,0)+cr2+cr3;
      ch(ido-1,1,k) = cc(0,k,0)+tr11*cr2+tr12*cr3;
      ch(0,2,k) = ti11*ci5+ti12*ci4;
      ch(ido-1,3,k) = cc(0,k,0)+tr12*cr2+tr11*cr3;
      ch(0,4,k) = ti12*ci5-ti11*ci4;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; k++) {
      for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
        std::size_t i1 = i0 + 1;
        ic = ido-i1;
        dr2 = wa1[i0-1]*cc(i0,k,1)+wa1[i0]*cc(i1,k,1);
        di2 = wa1[i0-1]*cc(i1,k,1)-wa1[i0]*cc(i0,k,1);
        dr3 = wa2[i0-1]*cc(i0,k,2)+wa2[i0]*cc(i1,k,2);
        di3 = wa2[i0-1]*cc(i1,k,2)-wa2[i0]*cc(i0,k,2);
        dr4 = wa3[i0-1]*cc(i0,k,3)+wa3[i0]*cc(i1,k,3);
        di4 = wa3[i0-1]*cc(i1,k,3)-wa3[i0]*cc(i0,k,3);
        dr5 = wa4[i0-1]*cc(i0,k,4)+wa4[i0]*cc(i1,k,4);
        di5 = wa4[i0-1]*cc(i1,k,4)-wa4[i0]*cc(i0,k,4);
        cr2 = dr2+dr5;
        ci5 = dr5-dr2;
        cr5 = di2-di5;
        ci2 = di2+di5;
        cr3 = dr3+dr4;
        ci4 = dr4-dr3;
        cr4 = di3-di4;
        ci3 = di3+di4;
        ch(i0,0,k) = cc(i0,k,0)+cr2+cr3;
        ch(i1,0,k) = cc(i1,k,0)+ci2+ci3;
        tr2 = cc(i0,k,0)+tr11*cr2+tr12*cr3;
        ti2 = cc(i1,k,0)+tr11*ci2+tr12*ci3;
        tr3 = cc(i0,k,0)+tr12*cr2+tr11*cr3;
        ti3 = cc(i1,k,0)+tr12*ci2+tr11*ci3;
        tr5 = ti11*cr5+ti12*cr4;
        ti5 = ti11*ci5+ti12*ci4;
        tr4 = ti12*cr5-ti11*cr4;
        ti4 = ti12*ci5-ti11*ci4;
        ch(i0,2,k) = tr2+tr5;
        ch(ic-1,1,k) = tr2-tr5;
        ch(i1,2,k) = ti2+ti5;
        ch(ic,1,k) = ti5-ti2;
        ch(i0,4,k) = tr3+tr4;
        ch(ic-1,3,k) = tr3-tr4;
        ch(i1,4,k) = ti3+ti4;
        ch(ic,3,k) = ti4-ti3;
      }
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passfg(std::size_t ido,
                                       std::size_t ip,
                                       std::size_t l1,
                                       std::size_t idl1,
                                       real_type* cc_begin,
                                       real_type* c1_begin,
                                       real_type* c2_begin,
                                       real_type* ch_begin,
                                       real_type* ch2_begin,
                                       const real_type* wa)
  {
    dim3 cc(cc_begin, ido, ip, l1);
    dim3 c1(c1_begin, ido, l1, ip);
    dim2 c2(c2_begin, idl1, ip);
    dim3 ch(ch_begin, ido, l1, ip);
    dim2 ch2(ch2_begin, idl1, ip);
    real_type ai1;
    real_type ai2;
    real_type ar1;
    real_type ar1h;
    real_type ar2;
    real_type ar2h;
    real_type dc2;
    real_type dcp;
    real_type ds2;
    real_type dsp;
    std::size_t ic;
    std::size_t idij;
    std::size_t ipph;
    std::size_t is;
    std::size_t j2;
    std::size_t jc;
    std::size_t lc;
    std::size_t nbd;
    const real_type tpi = real_type(8) * std::atan(real_type(1));
    real_type arg = tpi / real_type(ip);
    dcp = std::cos(arg);
    dsp = std::sin(arg);
    ipph = (ip+1)/2;
    nbd = (ido-1)/2;
    if (ido == 1) {
      for (std::size_t ik = 0; ik < idl1; ik++) {
        c2(ik,0) = ch2(ik,0);
      }
    }
    else {
      for (std::size_t ik = 0; ik < idl1; ik++) {
        ch2(ik,0) = c2(ik,0);
      }
      for (std::size_t j = 1; j < ip; j++) {
        for (std::size_t k = 0; k < l1; k++) {
          ch(0,k,j) = c1(0,k,j);
        }
      }
      if (nbd <= l1) {
        is = 0;
        for (std::size_t j = 1; j < ip; j++) {
          idij = is;
          for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
            std::size_t i1 = i0 + 1;
            for (std::size_t k = 0; k < l1; k++) {
              ch(i0,k,j) = wa[idij]*c1(i0,k,j)+wa[idij+1]*c1(i1,k,j);
              ch(i1,k,j) = wa[idij]*c1(i1,k,j)-wa[idij+1]*c1(i0,k,j);
            }
            idij += 2;
          }
          is += ido;
        }
      }
      else {
        is = 0;
        for (std::size_t j = 1; j < ip; j++) {
          for (std::size_t k = 0; k < l1; k++) {
            idij = is;
            for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
              std::size_t i1 = i0 + 1;
              ch(i0,k,j) = wa[idij]*c1(i0,k,j)+wa[idij+1]*c1(i1,k,j);
              ch(i1,k,j) = wa[idij]*c1(i1,k,j)-wa[idij+1]*c1(i0,k,j);
              idij = idij+2;
            }
          }
          is += ido;
        }
      }
      if (nbd >= l1) {
        for (std::size_t j = 1; j < ipph; j++) {
          jc = ip-j;
          for (std::size_t k = 0; k < l1; k++) {
            for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
              std::size_t i1 = i0 + 1;
              c1(i0,k,j) = ch(i0,k,j)+ch(i0,k,jc);
              c1(i0,k,jc) = ch(i1,k,j)-ch(i1,k,jc);
              c1(i1,k,j) = ch(i1,k,j)+ch(i1,k,jc);
              c1(i1,k,jc) = ch(i0,k,jc)-ch(i0,k,j);
            }
          }
        }
      }
      else {
        for (std::size_t j = 1; j < ipph; j++) {
          jc = ip-j;
          for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
            std::size_t i1 = i0 + 1;
            for (std::size_t k = 0; k < l1; k++) {
              c1(i0,k,j) = ch(i0,k,j)+ch(i0,k,jc);
              c1(i0,k,jc) = ch(i1,k,j)-ch(i1,k,jc);
              c1(i1,k,j) = ch(i1,k,j)+ch(i1,k,jc);
              c1(i1,k,jc) = ch(i0,k,jc)-ch(i0,k,j);
            }
          }
        }
      }
    }
    std::size_t j;
    for (j = 1; j < ipph; j++) {
      jc = ip-j;
      for (std::size_t k = 0; k < l1; k++) {
        c1(0,k,j) = ch(0,k,j)+ch(0,k,jc);
        c1(0,k,jc) = ch(0,k,jc)-ch(0,k,j);
      }
    }
    ar1 = 1.;
    ai1 = 0.;
    for (std::size_t l = 1; l < ipph; l++) {
      lc = ip-l;
      ar1h = dcp*ar1-dsp*ai1;
      ai1 = dcp*ai1+dsp*ar1;
      ar1 = ar1h;
      for (std::size_t ik = 0; ik < idl1; ik++) {
        ch2(ik,l) = c2(ik,0)+ar1*c2(ik,1);
        ch2(ik,lc) = ai1*c2(ik,ip-1);
      }
      dc2 = ar1;
      ds2 = ai1;
      ar2 = ar1;
      ai2 = ai1;
      for (std::size_t j = 2; j < ipph; j++) {
        jc = ip-j;
        ar2h = dc2*ar2-ds2*ai2;
        ai2 = dc2*ai2+ds2*ar2;
        ar2 = ar2h;
        for (std::size_t ik = 0; ik < idl1; ik++) {
          ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j);
          ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc);
        }
      }
    }
    for (j = 1; j < ipph; j++) {
      for (std::size_t ik = 0; ik < idl1; ik++) {
        ch2(ik,0) = ch2(ik,0)+c2(ik,j);
      }
    }
    if (ido >= l1) {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i = 0; i < ido; i++) {
          cc(i,0,k) = ch(i,k,0);
        }
      }
    }
    else {
      for (std::size_t i = 0; i < ido; i++) {
        for (std::size_t k = 0; k < l1; k++) {
          cc(i,0,k) = ch(i,k,0);
        }
      }
    }
    for (j = 1; j < ipph; j++) {
      jc = ip-j;
      j2 = j+j;
      for (std::size_t k = 0; k < l1; k++) {
        cc(ido-1,j2-1,k) = ch(0,k,j);
        cc(0,j2,k) = ch(0,k,jc);
      }
    }
    if (ido == 1) return;
    if (nbd >= l1) {
      for (std::size_t j = 1; j < ipph; j++) {
        jc = ip-j;
        j2 = j+j;
        for (std::size_t k = 0; k < l1; k++) {
          for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
            std::size_t i1 = i0 + 1;
            ic = ido-i1;
            cc(i0,j2,k) = ch(i0,k,j)+ch(i0,k,jc);
            cc(ic-1,j2-1,k) = ch(i0,k,j)-ch(i0,k,jc);
            cc(i1,j2,k) = ch(i1,k,j)+ch(i1,k,jc);
            cc(ic,j2-1,k) = ch(i1,k,jc)-ch(i1,k,j);
          }
        }
      }
    }
    else {
      for (std::size_t j = 1; j < ipph; j++) {
        jc = ip-j;
        j2 = j+j;
        for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
          std::size_t i1 = i0 + 1;
          ic = ido-i1;
          for (std::size_t k = 0; k < l1; k++) {
            cc(i0,j2,k) = ch(i0,k,j)+ch(i0,k,jc);
            cc(ic-1,j2-1,k) = ch(i0,k,j)-ch(i0,k,jc);
            cc(i1,j2,k) = ch(i1,k,j)+ch(i1,k,jc);
            cc(ic,j2-1,k) = ch(i1,k,jc)-ch(i1,k,j);
          }
        }
      }
    }
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_DETAIL_REAL_TO_COMPLEX_FORWARD_H
