#ifndef SCITBX_FFTPACK_DETAIL_REAL_TO_COMPLEX_BACKWARD_H
#define SCITBX_FFTPACK_DETAIL_REAL_TO_COMPLEX_BACKWARD_H

namespace scitbx { namespace fftpack {

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType, ComplexType>::backward_compressed(
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
    std::size_t l1;
    std::size_t l2;
    std::size_t na;
    na = 0;
    l1 = 1;
    iw = 1;
    for (std::size_t k1 = 0; k1 < factors_.size(); k1++) {
      ip = factors_[k1];
      l2 = ip*l1;
      ido = n_/l2;
      idl1 = ido*l1;
      if (ip == 4) {
        ix2 = iw+ido;
        ix3 = ix2+ido;
        if (na == 0) {
          passb4(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1);
        }
        else {
          passb4(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1);
        }
        na = 1-na;
      }
      else if (ip == 2) {
        if (na == 0) {
          passb2(ido,l1,c,ch,wa+iw-1);
        }
        else {
          passb2(ido,l1,ch,c,wa+iw-1);
        }
        na = 1-na;
      }
      else if (ip == 3) {
        ix2 = iw+ido;
        if (na == 0) {
          passb3(ido,l1,c,ch,wa+iw-1,wa+ix2-1);
        }
        else {
          passb3(ido,l1,ch,c,wa+iw-1,wa+ix2-1);
        }
        na = 1-na;
      }
      else if (ip == 5) {
        ix2 = iw+ido;
        ix3 = ix2+ido;
        ix4 = ix3+ido;
        if (na == 0) {
          passb5(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1,wa+ix4-1);
        }
        else {
          passb5(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1,wa+ix4-1);
        }
        na = 1-na;
      }
      else {
        if (na == 0) {
          passbg(ido,ip,l1,idl1,c,c,c,ch,ch,wa+iw-1);
        }
        else {
          passbg(ido,ip,l1,idl1,ch,ch,ch,c,c,wa+iw-1);
        }
        if (ido == 1) na = 1-na;
      }
      l1 = l2;
      iw = iw+(ip-1)*ido;
    }
    if (na == 0) return;
    for (std::size_t i = 0; i < n_; i++) {
      c[i] = ch[i];
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb2(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1)
  {
    dim3 cc(cc_begin, ido, 2, l1);
    dim3 ch(ch_begin, ido, l1, 2);
    std::size_t ic;
    real_type ti2;
    real_type tr2;
    std::size_t k;
    for (k = 0; k < l1; k++) {
      ch(0,k,0) = cc(0,0,k)+cc(ido-1,1,k);
      ch(0,k,1) = cc(0,0,k)-cc(ido-1,1,k);
    }
    if (ido < 2) return;
    if (ido > 2) {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
          std::size_t i1 = i0 + 1;
          ic = ido-i1;
          ch(i0,k,0) = cc(i0,0,k)+cc(ic-1,1,k);
          tr2 = cc(i0,0,k)-cc(ic-1,1,k);
          ch(i1,k,0) = cc(i1,0,k)-cc(ic,1,k);
          ti2 = cc(i1,0,k)+cc(ic,1,k);
          ch(i0,k,1) = wa1[i0-1]*tr2-wa1[i0]*ti2;
          ch(i1,k,1) = wa1[i0-1]*ti2+wa1[i0]*tr2;
        }
      }
      if (ido % 2 != 0) return;
    }
    for (k = 0; k < l1; k++) {
      ch(ido-1,k,0) = cc(ido-1,0,k)+cc(ido-1,0,k);
      ch(ido-1,k,1) = -(cc(0,1,k)+cc(0,1,k));
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb3(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1,
                                       const real_type* wa2)
  {
    dim3 cc(cc_begin, ido, 3, l1);
    dim3 ch(ch_begin, ido, l1, 3);
    real_type ci2;
    real_type ci3;
    real_type cr2;
    real_type cr3;
    real_type di2;
    real_type di3;
    real_type dr2;
    real_type dr3;
    std::size_t ic;
    real_type ti2;
    real_type tr2;
    const real_type taur = -.5;
    const real_type taui = -taur * std::sqrt(real_type(3));
    std::size_t k;
    for (k = 0; k < l1; k++) {
      tr2 = cc(ido-1,1,k)+cc(ido-1,1,k);
      cr2 = cc(0,0,k)+taur*tr2;
      ch(0,k,0) = cc(0,0,k)+tr2;
      ci3 = taui*(cc(0,2,k)+cc(0,2,k));
      ch(0,k,1) = cr2-ci3;
      ch(0,k,2) = cr2+ci3;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; k++) {
      for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
        std::size_t i1 = i0 + 1;
        ic = ido-i1;
        tr2 = cc(i0,2,k)+cc(ic-1,1,k);
        cr2 = cc(i0,0,k)+taur*tr2;
        ch(i0,k,0) = cc(i0,0,k)+tr2;
        ti2 = cc(i1,2,k)-cc(ic,1,k);
        ci2 = cc(i1,0,k)+taur*ti2;
        ch(i1,k,0) = cc(i1,0,k)+ti2;
        cr3 = taui*(cc(i0,2,k)-cc(ic-1,1,k));
        ci3 = taui*(cc(i1,2,k)+cc(ic,1,k));
        dr2 = cr2-ci3;
        dr3 = cr2+ci3;
        di2 = ci2+cr3;
        di3 = ci2-cr3;
        ch(i0,k,1) = wa1[i0-1]*dr2-wa1[i0]*di2;
        ch(i1,k,1) = wa1[i0-1]*di2+wa1[i0]*dr2;
        ch(i0,k,2) = wa2[i0-1]*dr3-wa2[i0]*di3;
        ch(i1,k,2) = wa2[i0-1]*di3+wa2[i0]*dr3;
      }
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb4(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1,
                                       const real_type* wa2,
                                       const real_type* wa3)
  {
    dim3 cc(cc_begin, ido, 4, l1);
    dim3 ch(ch_begin, ido, l1, 4);
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
    real_type sqrt2 = std::sqrt(real_type(2));
    std::size_t k;
    for (k = 0; k < l1; k++) {
      tr1 = cc(0,0,k)-cc(ido-1,3,k);
      tr2 = cc(0,0,k)+cc(ido-1,3,k);
      tr3 = cc(ido-1,1,k)+cc(ido-1,1,k);
      tr4 = cc(0,2,k)+cc(0,2,k);
      ch(0,k,0) = tr2+tr3;
      ch(0,k,1) = tr1-tr4;
      ch(0,k,2) = tr2-tr3;
      ch(0,k,3) = tr1+tr4;
    }
    if (ido < 2) return;
    if (ido > 2) {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
          std::size_t i1 = i0 + 1;
          ic = ido-i1;
          ti1 = cc(i1,0,k)+cc(ic,3,k);
          ti2 = cc(i1,0,k)-cc(ic,3,k);
          ti3 = cc(i1,2,k)-cc(ic,1,k);
          tr4 = cc(i1,2,k)+cc(ic,1,k);
          tr1 = cc(i0,0,k)-cc(ic-1,3,k);
          tr2 = cc(i0,0,k)+cc(ic-1,3,k);
          ti4 = cc(i0,2,k)-cc(ic-1,1,k);
          tr3 = cc(i0,2,k)+cc(ic-1,1,k);
          ch(i0,k,0) = tr2+tr3;
          cr3 = tr2-tr3;
          ch(i1,k,0) = ti2+ti3;
          ci3 = ti2-ti3;
          cr2 = tr1-tr4;
          cr4 = tr1+tr4;
          ci2 = ti1+ti4;
          ci4 = ti1-ti4;
          ch(i0,k,1) = wa1[i0-1]*cr2-wa1[i0]*ci2;
          ch(i1,k,1) = wa1[i0-1]*ci2+wa1[i0]*cr2;
          ch(i0,k,2) = wa2[i0-1]*cr3-wa2[i0]*ci3;
          ch(i1,k,2) = wa2[i0-1]*ci3+wa2[i0]*cr3;
          ch(i0,k,3) = wa3[i0-1]*cr4-wa3[i0]*ci4;
          ch(i1,k,3) = wa3[i0-1]*ci4+wa3[i0]*cr4;
        }
      }
      if (ido % 2 != 0) return;
    }
    for (k = 0; k < l1; k++) {
      ti1 = cc(0,1,k)+cc(0,3,k);
      ti2 = cc(0,3,k)-cc(0,1,k);
      tr1 = cc(ido-1,0,k)-cc(ido-1,2,k);
      tr2 = cc(ido-1,0,k)+cc(ido-1,2,k);
      ch(ido-1,k,0) = tr2+tr2;
      ch(ido-1,k,1) = sqrt2*(tr1-ti1);
      ch(ido-1,k,2) = ti2+ti2;
      ch(ido-1,k,3) = -sqrt2*(tr1+ti1);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb5(std::size_t ido,
                                       std::size_t l1,
                                       real_type* cc_begin,
                                       real_type* ch_begin,
                                       const real_type* wa1,
                                       const real_type* wa2,
                                       const real_type* wa3,
                                       const real_type* wa4)
  {
    dim3 cc(cc_begin, ido, 5, l1);
    dim3 ch(ch_begin, ido, l1, 5);
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
      ti5 = cc(0,2,k)+cc(0,2,k);
      ti4 = cc(0,4,k)+cc(0,4,k);
      tr2 = cc(ido-1,1,k)+cc(ido-1,1,k);
      tr3 = cc(ido-1,3,k)+cc(ido-1,3,k);
      ch(0,k,0) = cc(0,0,k)+tr2+tr3;
      cr2 = cc(0,0,k)+tr11*tr2+tr12*tr3;
      cr3 = cc(0,0,k)+tr12*tr2+tr11*tr3;
      ci5 = ti11*ti5+ti12*ti4;
      ci4 = ti12*ti5-ti11*ti4;
      ch(0,k,1) = cr2-ci5;
      ch(0,k,2) = cr3-ci4;
      ch(0,k,3) = cr3+ci4;
      ch(0,k,4) = cr2+ci5;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; k++) {
      for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
        std::size_t i1 = i0 + 1;
        ic = ido-i1;
        ti5 = cc(i1,2,k)+cc(ic,1,k);
        ti2 = cc(i1,2,k)-cc(ic,1,k);
        ti4 = cc(i1,4,k)+cc(ic,3,k);
        ti3 = cc(i1,4,k)-cc(ic,3,k);
        tr5 = cc(i0,2,k)-cc(ic-1,1,k);
        tr2 = cc(i0,2,k)+cc(ic-1,1,k);
        tr4 = cc(i0,4,k)-cc(ic-1,3,k);
        tr3 = cc(i0,4,k)+cc(ic-1,3,k);
        ch(i0,k,0) = cc(i0,0,k)+tr2+tr3;
        ch(i1,k,0) = cc(i1,0,k)+ti2+ti3;
        cr2 = cc(i0,0,k)+tr11*tr2+tr12*tr3;
        ci2 = cc(i1,0,k)+tr11*ti2+tr12*ti3;
        cr3 = cc(i0,0,k)+tr12*tr2+tr11*tr3;
        ci3 = cc(i1,0,k)+tr12*ti2+tr11*ti3;
        cr5 = ti11*tr5+ti12*tr4;
        ci5 = ti11*ti5+ti12*ti4;
        cr4 = ti12*tr5-ti11*tr4;
        ci4 = ti12*ti5-ti11*ti4;
        dr3 = cr3-ci4;
        dr4 = cr3+ci4;
        di3 = ci3+cr4;
        di4 = ci3-cr4;
        dr5 = cr2+ci5;
        dr2 = cr2-ci5;
        di5 = ci2-cr5;
        di2 = ci2+cr5;
        ch(i0,k,1) = wa1[i0-1]*dr2-wa1[i0]*di2;
        ch(i1,k,1) = wa1[i0-1]*di2+wa1[i0]*dr2;
        ch(i0,k,2) = wa2[i0-1]*dr3-wa2[i0]*di3;
        ch(i1,k,2) = wa2[i0-1]*di3+wa2[i0]*dr3;
        ch(i0,k,3) = wa3[i0-1]*dr4-wa3[i0]*di4;
        ch(i1,k,3) = wa3[i0-1]*di4+wa3[i0]*dr4;
        ch(i0,k,4) = wa4[i0-1]*dr5-wa4[i0]*di5;
        ch(i1,k,4) = wa4[i0-1]*di5+wa4[i0]*dr5;
      }
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passbg(std::size_t ido,
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
    nbd = (ido-1)/2;
    ipph = (ip+1)/2;
    if (ido >= l1) {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i = 0; i < ido; i++) {
          ch(i,k,0) = cc(i,0,k);
        }
      }
    }
    else {
      for (std::size_t i = 0; i < ido; i++) {
        for (std::size_t k = 0; k < l1; k++) {
          ch(i,k,0) = cc(i,0,k);
        }
      }
    }
    std::size_t j;
    for (j = 1; j < ipph; j++) {
       jc = ip-j;
       j2 = j+j;
       for (std::size_t k = 0; k < l1; k++) {
         ch(0,k,j) = cc(ido-1,j2-1,k)+cc(ido-1,j2-1,k);
         ch(0,k,jc) = cc(0,j2,k)+cc(0,j2,k);
      }
    }
    if (ido != 1) {
      if (nbd >= l1) {
        for (std::size_t j = 1; j < ipph; j++) {
          jc = ip-j;
          j2 = j+j;
          for (std::size_t k = 0; k < l1; k++) {
            for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
              std::size_t i1 = i0 + 1;
              ic = ido-i1;
              ch(i0,k,j) = cc(i0,j2,k)+cc(ic-1,j2-1,k);
              ch(i0,k,jc) = cc(i0,j2,k)-cc(ic-1,j2-1,k);
              ch(i1,k,j) = cc(i1,j2,k)-cc(ic,j2-1,k);
              ch(i1,k,jc) = cc(i1,j2,k)+cc(ic,j2-1,k);
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
              ch(i0,k,j) = cc(i0,j2,k)+cc(ic-1,j2-1,k);
              ch(i0,k,jc) = cc(i0,j2,k)-cc(ic-1,j2-1,k);
              ch(i1,k,j) = cc(i1,j2,k)-cc(ic,j2-1,k);
              ch(i1,k,jc) = cc(i1,j2,k)+cc(ic,j2-1,k);
            }
          }
        }
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
        c2(ik,l) = ch2(ik,0)+ar1*ch2(ik,1);
        c2(ik,lc) = ai1*ch2(ik,ip-1);
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
          c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j);
          c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc);
        }
      }
    }
    for (j = 1; j < ipph; j++) {
      for (std::size_t ik = 0; ik < idl1; ik++) {
        ch2(ik,0) += ch2(ik,j);
      }
    }
    for (j = 1; j < ipph; j++) {
      jc = ip-j;
      for (std::size_t k = 0; k < l1; k++) {
        ch(0,k,j) = c1(0,k,j)-c1(0,k,jc);
        ch(0,k,jc) = c1(0,k,j)+c1(0,k,jc);
      }
    }
    if (ido != 1) {
      if (nbd >= l1) {
        for (std::size_t j = 1; j < ipph; j++) {
          jc = ip-j;
          for (std::size_t k = 0; k < l1; k++) {
            for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
              std::size_t i1 = i0 + 1;
              ch(i0,k,j) = c1(i0,k,j)-c1(i1,k,jc);
              ch(i0,k,jc) = c1(i0,k,j)+c1(i1,k,jc);
              ch(i1,k,j) = c1(i1,k,j)+c1(i0,k,jc);
              ch(i1,k,jc) = c1(i1,k,j)-c1(i0,k,jc);
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
              ch(i0,k,j) = c1(i0,k,j)-c1(i1,k,jc);
              ch(i0,k,jc) = c1(i0,k,j)+c1(i1,k,jc);
              ch(i1,k,j) = c1(i1,k,j)+c1(i0,k,jc);
              ch(i1,k,jc) = c1(i1,k,j)-c1(i0,k,jc);
            }
          }
        }
      }
    }
    if (ido == 1) return;
    for (std::size_t ik = 0; ik < idl1; ik++) {
      c2(ik,0) = ch2(ik,0);
    }
    for (j = 1; j < ip; j++) {
      for (std::size_t k = 0; k < l1; k++) {
        c1(0,k,j) = ch(0,k,j);
      }
    }
    if (nbd <= l1) {
      is = 0;
      for (std::size_t j = 1; j < ip; j++) {
        idij = is;
        for (std::size_t i0 = 1; i0 < ido-1; i0 += 2) {
          std::size_t i1 = i0 + 1;
          for (std::size_t k = 0; k < l1; k++) {
            c1(i0,k,j) = wa[idij]*ch(i0,k,j)-wa[idij+1]*ch(i1,k,j);
            c1(i1,k,j) = wa[idij]*ch(i1,k,j)+wa[idij+1]*ch(i0,k,j);
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
            c1(i0,k,j) = wa[idij]*ch(i0,k,j)-wa[idij+1]*ch(i1,k,j);
            c1(i1,k,j) = wa[idij]*ch(i1,k,j)+wa[idij+1]*ch(i0,k,j);
            idij += 2;
          }
        }
        is += ido;
      }
    }
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_DETAIL_REAL_TO_COMPLEX_BACKWARD_H
