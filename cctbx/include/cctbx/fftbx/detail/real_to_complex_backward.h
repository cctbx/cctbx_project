// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Dec 21: iterator-based interface (rwgk)
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_DETAIL_REAL_TO_COMPLEX_BACKWARD_H
#define CCTBX_FFTBX_DETAIL_REAL_TO_COMPLEX_BACKWARD_H

namespace cctbx { namespace fftbx {

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::backward_compressed(real_type* Seq_begin)
  {
    if (m_N < 2) return;
    real_type* C = Seq_begin;
    real_type* CH = &(*(m_CH.begin()));
    const real_type* WA = &(*(m_WA.begin()));
    std::size_t IDL1;
    std::size_t IDO;
    std::size_t IP;
    std::size_t IW;
    std::size_t IX2;
    std::size_t IX3;
    std::size_t IX4;
    std::size_t L1;
    std::size_t L2;
    std::size_t NA;
    NA = 0;
    L1 = 1;
    IW = 1;
    for (std::size_t K1 = 0; K1 < m_Factors.size(); K1++) {
      IP = m_Factors[K1];
      L2 = IP*L1;
      IDO = m_N/L2;
      IDL1 = IDO*L1;
      if (IP == 4) {
        IX2 = IW+IDO;
        IX3 = IX2+IDO;
        if (NA == 0) {
          passb4(IDO,L1,C,CH,WA+IW-1,WA+IX2-1,WA+IX3-1);
        }
        else {
          passb4(IDO,L1,CH,C,WA+IW-1,WA+IX2-1,WA+IX3-1);
        }
        NA = 1-NA;
      }
      else if (IP == 2) {
        if (NA == 0) {
          passb2(IDO,L1,C,CH,WA+IW-1);
        }
        else {
          passb2(IDO,L1,CH,C,WA+IW-1);
        }
        NA = 1-NA;
      }
      else if (IP == 3) {
        IX2 = IW+IDO;
        if (NA == 0) {
          passb3(IDO,L1,C,CH,WA+IW-1,WA+IX2-1);
        }
        else {
          passb3(IDO,L1,CH,C,WA+IW-1,WA+IX2-1);
        }
        NA = 1-NA;
      }
      else if (IP == 5) {
        IX2 = IW+IDO;
        IX3 = IX2+IDO;
        IX4 = IX3+IDO;
        if (NA == 0) {
          passb5(IDO,L1,C,CH,WA+IW-1,WA+IX2-1,WA+IX3-1,WA+IX4-1);
        }
        else {
          passb5(IDO,L1,CH,C,WA+IW-1,WA+IX2-1,WA+IX3-1,WA+IX4-1);
        }
        NA = 1-NA;
      }
      else {
        if (NA == 0) {
          passbg(IDO,IP,L1,IDL1,C,C,C,CH,CH,WA+IW-1);
        }
        else {
          passbg(IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA+IW-1);
        }
        if (IDO == 1) NA = 1-NA;
      }
      L1 = L2;
      IW = IW+(IP-1)*IDO;
    }
    if (NA == 0) return;
    for (std::size_t I = 0; I < m_N; I++) {
      C[I] = CH[I];
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb2(std::size_t IDO,
                                       std::size_t L1,
                                       real_type* CC_begin,
                                       real_type* CH_begin,
                                       const real_type* WA1)
  {
    dim3 CC(CC_begin, IDO, 2, L1);
    dim3 CH(CH_begin, IDO, L1, 2);
    std::size_t IC;
    real_type TI2;
    real_type TR2;
    std::size_t K;
    for (K = 0; K < L1; K++) {
      CH(0,K,0) = CC(0,0,K)+CC(IDO-1,1,K);
      CH(0,K,1) = CC(0,0,K)-CC(IDO-1,1,K);
    }
    if (IDO < 2) return;
    if (IDO > 2) {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
          std::size_t I1 = I0 + 1;
          IC = IDO-I1;
          CH(I0,K,0) = CC(I0,0,K)+CC(IC-1,1,K);
          TR2 = CC(I0,0,K)-CC(IC-1,1,K);
          CH(I1,K,0) = CC(I1,0,K)-CC(IC,1,K);
          TI2 = CC(I1,0,K)+CC(IC,1,K);
          CH(I0,K,1) = WA1[I0-1]*TR2-WA1[I0]*TI2;
          CH(I1,K,1) = WA1[I0-1]*TI2+WA1[I0]*TR2;
        }
      }
      if (IDO % 2 != 0) return;
    }
    for (K = 0; K < L1; K++) {
      CH(IDO-1,K,0) = CC(IDO-1,0,K)+CC(IDO-1,0,K);
      CH(IDO-1,K,1) = -(CC(0,1,K)+CC(0,1,K));
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb3(std::size_t IDO,
                                       std::size_t L1,
                                       real_type* CC_begin,
                                       real_type* CH_begin,
                                       const real_type* WA1,
                                       const real_type* WA2)
  {
    dim3 CC(CC_begin, IDO, 3, L1);
    dim3 CH(CH_begin, IDO, L1, 3);
    real_type CI2;
    real_type CI3;
    real_type CR2;
    real_type CR3;
    real_type DI2;
    real_type DI3;
    real_type DR2;
    real_type DR3;
    std::size_t IC;
    real_type TI2;
    real_type TR2;
    const real_type TAUR = -.5;
    const real_type TAUI = -TAUR * std::sqrt(real_type(3));
    std::size_t K;
    for (K = 0; K < L1; K++) {
      TR2 = CC(IDO-1,1,K)+CC(IDO-1,1,K);
      CR2 = CC(0,0,K)+TAUR*TR2;
      CH(0,K,0) = CC(0,0,K)+TR2;
      CI3 = TAUI*(CC(0,2,K)+CC(0,2,K));
      CH(0,K,1) = CR2-CI3;
      CH(0,K,2) = CR2+CI3;
    }
    if (IDO == 1) return;
    for (K = 0; K < L1; K++) {
      for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
        std::size_t I1 = I0 + 1;
        IC = IDO-I1;
        TR2 = CC(I0,2,K)+CC(IC-1,1,K);
        CR2 = CC(I0,0,K)+TAUR*TR2;
        CH(I0,K,0) = CC(I0,0,K)+TR2;
        TI2 = CC(I1,2,K)-CC(IC,1,K);
        CI2 = CC(I1,0,K)+TAUR*TI2;
        CH(I1,K,0) = CC(I1,0,K)+TI2;
        CR3 = TAUI*(CC(I0,2,K)-CC(IC-1,1,K));
        CI3 = TAUI*(CC(I1,2,K)+CC(IC,1,K));
        DR2 = CR2-CI3;
        DR3 = CR2+CI3;
        DI2 = CI2+CR3;
        DI3 = CI2-CR3;
        CH(I0,K,1) = WA1[I0-1]*DR2-WA1[I0]*DI2;
        CH(I1,K,1) = WA1[I0-1]*DI2+WA1[I0]*DR2;
        CH(I0,K,2) = WA2[I0-1]*DR3-WA2[I0]*DI3;
        CH(I1,K,2) = WA2[I0-1]*DI3+WA2[I0]*DR3;
      }
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb4(std::size_t IDO,
                                       std::size_t L1,
                                       real_type* CC_begin,
                                       real_type* CH_begin,
                                       const real_type* WA1,
                                       const real_type* WA2,
                                       const real_type* WA3)
  {
    dim3 CC(CC_begin, IDO, 4, L1);
    dim3 CH(CH_begin, IDO, L1, 4);
    real_type CI2;
    real_type CI3;
    real_type CI4;
    real_type CR2;
    real_type CR3;
    real_type CR4;
    std::size_t IC;
    real_type TI1;
    real_type TI2;
    real_type TI3;
    real_type TI4;
    real_type TR1;
    real_type TR2;
    real_type TR3;
    real_type TR4;
    real_type SQRT2 = std::sqrt(real_type(2));
    std::size_t K;
    for (K = 0; K < L1; K++) {
      TR1 = CC(0,0,K)-CC(IDO-1,3,K);
      TR2 = CC(0,0,K)+CC(IDO-1,3,K);
      TR3 = CC(IDO-1,1,K)+CC(IDO-1,1,K);
      TR4 = CC(0,2,K)+CC(0,2,K);
      CH(0,K,0) = TR2+TR3;
      CH(0,K,1) = TR1-TR4;
      CH(0,K,2) = TR2-TR3;
      CH(0,K,3) = TR1+TR4;
    }
    if (IDO < 2) return;
    if (IDO > 2) {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
          std::size_t I1 = I0 + 1;
          IC = IDO-I1;
          TI1 = CC(I1,0,K)+CC(IC,3,K);
          TI2 = CC(I1,0,K)-CC(IC,3,K);
          TI3 = CC(I1,2,K)-CC(IC,1,K);
          TR4 = CC(I1,2,K)+CC(IC,1,K);
          TR1 = CC(I0,0,K)-CC(IC-1,3,K);
          TR2 = CC(I0,0,K)+CC(IC-1,3,K);
          TI4 = CC(I0,2,K)-CC(IC-1,1,K);
          TR3 = CC(I0,2,K)+CC(IC-1,1,K);
          CH(I0,K,0) = TR2+TR3;
          CR3 = TR2-TR3;
          CH(I1,K,0) = TI2+TI3;
          CI3 = TI2-TI3;
          CR2 = TR1-TR4;
          CR4 = TR1+TR4;
          CI2 = TI1+TI4;
          CI4 = TI1-TI4;
          CH(I0,K,1) = WA1[I0-1]*CR2-WA1[I0]*CI2;
          CH(I1,K,1) = WA1[I0-1]*CI2+WA1[I0]*CR2;
          CH(I0,K,2) = WA2[I0-1]*CR3-WA2[I0]*CI3;
          CH(I1,K,2) = WA2[I0-1]*CI3+WA2[I0]*CR3;
          CH(I0,K,3) = WA3[I0-1]*CR4-WA3[I0]*CI4;
          CH(I1,K,3) = WA3[I0-1]*CI4+WA3[I0]*CR4;
        }
      }
      if (IDO % 2 != 0) return;
    }
    for (K = 0; K < L1; K++) {
      TI1 = CC(0,1,K)+CC(0,3,K);
      TI2 = CC(0,3,K)-CC(0,1,K);
      TR1 = CC(IDO-1,0,K)-CC(IDO-1,2,K);
      TR2 = CC(IDO-1,0,K)+CC(IDO-1,2,K);
      CH(IDO-1,K,0) = TR2+TR2;
      CH(IDO-1,K,1) = SQRT2*(TR1-TI1);
      CH(IDO-1,K,2) = TI2+TI2;
      CH(IDO-1,K,3) = -SQRT2*(TR1+TI1);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passb5(std::size_t IDO,
                                       std::size_t L1,
                                       real_type* CC_begin,
                                       real_type* CH_begin,
                                       const real_type* WA1,
                                       const real_type* WA2,
                                       const real_type* WA3,
                                       const real_type* WA4)
  {
    dim3 CC(CC_begin, IDO, 5, L1);
    dim3 CH(CH_begin, IDO, L1, 5);
    real_type CI2;
    real_type CI3;
    real_type CI4;
    real_type CI5;
    real_type CR2;
    real_type CR3;
    real_type CR4;
    real_type CR5;
    real_type DI2;
    real_type DI3;
    real_type DI4;
    real_type DI5;
    real_type DR2;
    real_type DR3;
    real_type DR4;
    real_type DR5;
    std::size_t IC;
    real_type TI2;
    real_type TI3;
    real_type TI4;
    real_type TI5;
    real_type TR2;
    real_type TR3;
    real_type TR4;
    real_type TR5;
    // sin(18 deg)
    const real_type TR11 =  .309016994374947424102293417182819058860154589903;
    // cos(18 deg)
    const real_type TI11 =  .951056516295153572116439333379382143405698634126;
    // -cos(36 deg)
    const real_type TR12 = -.809016994374947424102293417182819058860154589903;
    // sin(36 deg)
    const real_type TI12 =  .587785252292473129168705954639072768597652437643;
    std::size_t K;
    for (K = 0; K < L1; K++) {
      TI5 = CC(0,2,K)+CC(0,2,K);
      TI4 = CC(0,4,K)+CC(0,4,K);
      TR2 = CC(IDO-1,1,K)+CC(IDO-1,1,K);
      TR3 = CC(IDO-1,3,K)+CC(IDO-1,3,K);
      CH(0,K,0) = CC(0,0,K)+TR2+TR3;
      CR2 = CC(0,0,K)+TR11*TR2+TR12*TR3;
      CR3 = CC(0,0,K)+TR12*TR2+TR11*TR3;
      CI5 = TI11*TI5+TI12*TI4;
      CI4 = TI12*TI5-TI11*TI4;
      CH(0,K,1) = CR2-CI5;
      CH(0,K,2) = CR3-CI4;
      CH(0,K,3) = CR3+CI4;
      CH(0,K,4) = CR2+CI5;
    }
    if (IDO == 1) return;
    for (K = 0; K < L1; K++) {
      for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
        std::size_t I1 = I0 + 1;
        IC = IDO-I1;
        TI5 = CC(I1,2,K)+CC(IC,1,K);
        TI2 = CC(I1,2,K)-CC(IC,1,K);
        TI4 = CC(I1,4,K)+CC(IC,3,K);
        TI3 = CC(I1,4,K)-CC(IC,3,K);
        TR5 = CC(I0,2,K)-CC(IC-1,1,K);
        TR2 = CC(I0,2,K)+CC(IC-1,1,K);
        TR4 = CC(I0,4,K)-CC(IC-1,3,K);
        TR3 = CC(I0,4,K)+CC(IC-1,3,K);
        CH(I0,K,0) = CC(I0,0,K)+TR2+TR3;
        CH(I1,K,0) = CC(I1,0,K)+TI2+TI3;
        CR2 = CC(I0,0,K)+TR11*TR2+TR12*TR3;
        CI2 = CC(I1,0,K)+TR11*TI2+TR12*TI3;
        CR3 = CC(I0,0,K)+TR12*TR2+TR11*TR3;
        CI3 = CC(I1,0,K)+TR12*TI2+TR11*TI3;
        CR5 = TI11*TR5+TI12*TR4;
        CI5 = TI11*TI5+TI12*TI4;
        CR4 = TI12*TR5-TI11*TR4;
        CI4 = TI12*TI5-TI11*TI4;
        DR3 = CR3-CI4;
        DR4 = CR3+CI4;
        DI3 = CI3+CR4;
        DI4 = CI3-CR4;
        DR5 = CR2+CI5;
        DR2 = CR2-CI5;
        DI5 = CI2-CR5;
        DI2 = CI2+CR5;
        CH(I0,K,1) = WA1[I0-1]*DR2-WA1[I0]*DI2;
        CH(I1,K,1) = WA1[I0-1]*DI2+WA1[I0]*DR2;
        CH(I0,K,2) = WA2[I0-1]*DR3-WA2[I0]*DI3;
        CH(I1,K,2) = WA2[I0-1]*DI3+WA2[I0]*DR3;
        CH(I0,K,3) = WA3[I0-1]*DR4-WA3[I0]*DI4;
        CH(I1,K,3) = WA3[I0-1]*DI4+WA3[I0]*DR4;
        CH(I0,K,4) = WA4[I0-1]*DR5-WA4[I0]*DI5;
        CH(I1,K,4) = WA4[I0-1]*DI5+WA4[I0]*DR5;
      }
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::passbg(std::size_t IDO,
                                       std::size_t IP,
                                       std::size_t L1,
                                       std::size_t IDL1,
                                       real_type* CC_begin,
                                       real_type* C1_begin,
                                       real_type* C2_begin,
                                       real_type* CH_begin,
                                       real_type* CH2_begin,
                                       const real_type* WA)
  {
    dim3 CC(CC_begin, IDO, IP, L1);
    dim3 C1(C1_begin, IDO, L1, IP);
    dim2 C2(C2_begin, IDL1, IP);
    dim3 CH(CH_begin, IDO, L1, IP);
    dim2 CH2(CH2_begin, IDL1, IP);
    real_type AI1;
    real_type AI2;
    real_type AR1;
    real_type AR1H;
    real_type AR2;
    real_type AR2H;
    real_type DC2;
    real_type DCP;
    real_type DS2;
    real_type DSP;
    std::size_t IC;
    std::size_t IDIJ;
    std::size_t IPPH;
    std::size_t IS;
    std::size_t J2;
    std::size_t JC;
    std::size_t LC;
    std::size_t NBD;
    const real_type TPI = real_type(8) * std::atan(real_type(1));
    real_type ARG = TPI / real_type(IP);
    DCP = std::cos(ARG);
    DSP = std::sin(ARG);
    NBD = (IDO-1)/2;
    IPPH = (IP+1)/2;
    if (IDO >= L1) {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I = 0; I < IDO; I++) {
          CH(I,K,0) = CC(I,0,K);
        }
      }
    }
    else {
      for (std::size_t I = 0; I < IDO; I++) {
        for (std::size_t K = 0; K < L1; K++) {
          CH(I,K,0) = CC(I,0,K);
        }
      }
    }
    std::size_t J;
    for (J = 1; J < IPPH; J++) {
       JC = IP-J;
       J2 = J+J;
       for (std::size_t K = 0; K < L1; K++) {
         CH(0,K,J) = CC(IDO-1,J2-1,K)+CC(IDO-1,J2-1,K);
         CH(0,K,JC) = CC(0,J2,K)+CC(0,J2,K);
      }
    }
    if (IDO != 1) {
      if (NBD >= L1) {
        for (std::size_t J = 1; J < IPPH; J++) {
          JC = IP-J;
          J2 = J+J;
          for (std::size_t K = 0; K < L1; K++) {
            for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
              std::size_t I1 = I0 + 1;
              IC = IDO-I1;
              CH(I0,K,J) = CC(I0,J2,K)+CC(IC-1,J2-1,K);
              CH(I0,K,JC) = CC(I0,J2,K)-CC(IC-1,J2-1,K);
              CH(I1,K,J) = CC(I1,J2,K)-CC(IC,J2-1,K);
              CH(I1,K,JC) = CC(I1,J2,K)+CC(IC,J2-1,K);
            }
          }
        }
      }
      else {
        for (std::size_t J = 1; J < IPPH; J++) {
          JC = IP-J;
          J2 = J+J;
          for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
            std::size_t I1 = I0 + 1;
            IC = IDO-I1;
            for (std::size_t K = 0; K < L1; K++) {
              CH(I0,K,J) = CC(I0,J2,K)+CC(IC-1,J2-1,K);
              CH(I0,K,JC) = CC(I0,J2,K)-CC(IC-1,J2-1,K);
              CH(I1,K,J) = CC(I1,J2,K)-CC(IC,J2-1,K);
              CH(I1,K,JC) = CC(I1,J2,K)+CC(IC,J2-1,K);
            }
          }
        }
      }
    }
    AR1 = 1.;
    AI1 = 0.;
    for (std::size_t L = 1; L < IPPH; L++) {
      LC = IP-L;
      AR1H = DCP*AR1-DSP*AI1;
      AI1 = DCP*AI1+DSP*AR1;
      AR1 = AR1H;
      for (std::size_t IK = 0; IK < IDL1; IK++) {
        C2(IK,L) = CH2(IK,0)+AR1*CH2(IK,1);
        C2(IK,LC) = AI1*CH2(IK,IP-1);
      }
      DC2 = AR1;
      DS2 = AI1;
      AR2 = AR1;
      AI2 = AI1;
      for (std::size_t J = 2; J < IPPH; J++) {
        JC = IP-J;
        AR2H = DC2*AR2-DS2*AI2;
        AI2 = DC2*AI2+DS2*AR2;
        AR2 = AR2H;
        for (std::size_t IK = 0; IK < IDL1; IK++) {
          C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J);
          C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC);
        }
      }
    }
    for (J = 1; J < IPPH; J++) {
      for (std::size_t IK = 0; IK < IDL1; IK++) {
        CH2(IK,0) += CH2(IK,J);
      }
    }
    for (J = 1; J < IPPH; J++) {
      JC = IP-J;
      for (std::size_t K = 0; K < L1; K++) {
        CH(0,K,J) = C1(0,K,J)-C1(0,K,JC);
        CH(0,K,JC) = C1(0,K,J)+C1(0,K,JC);
      }
    }
    if (IDO != 1) {
      if (NBD >= L1) {
        for (std::size_t J = 1; J < IPPH; J++) {
          JC = IP-J;
          for (std::size_t K = 0; K < L1; K++) {
            for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
              std::size_t I1 = I0 + 1;
              CH(I0,K,J) = C1(I0,K,J)-C1(I1,K,JC);
              CH(I0,K,JC) = C1(I0,K,J)+C1(I1,K,JC);
              CH(I1,K,J) = C1(I1,K,J)+C1(I0,K,JC);
              CH(I1,K,JC) = C1(I1,K,J)-C1(I0,K,JC);
            }
          }
        }
      }
      else {
        for (std::size_t J = 1; J < IPPH; J++) {
          JC = IP-J;
          for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
            std::size_t I1 = I0 + 1;
            for (std::size_t K = 0; K < L1; K++) {
              CH(I0,K,J) = C1(I0,K,J)-C1(I1,K,JC);
              CH(I0,K,JC) = C1(I0,K,J)+C1(I1,K,JC);
              CH(I1,K,J) = C1(I1,K,J)+C1(I0,K,JC);
              CH(I1,K,JC) = C1(I1,K,J)-C1(I0,K,JC);
            }
          }
        }
      }
    }
    if (IDO == 1) return;
    for (std::size_t IK = 0; IK < IDL1; IK++) {
      C2(IK,0) = CH2(IK,0);
    }
    for (J = 1; J < IP; J++) {
      for (std::size_t K = 0; K < L1; K++) {
        C1(0,K,J) = CH(0,K,J);
      }
    }
    if (NBD <= L1) {
      IS = 0;
      for (std::size_t J = 1; J < IP; J++) {
        IDIJ = IS;
        for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
          std::size_t I1 = I0 + 1;
          for (std::size_t K = 0; K < L1; K++) {
            C1(I0,K,J) = WA[IDIJ]*CH(I0,K,J)-WA[IDIJ+1]*CH(I1,K,J);
            C1(I1,K,J) = WA[IDIJ]*CH(I1,K,J)+WA[IDIJ+1]*CH(I0,K,J);
          }
          IDIJ += 2;
        }
        IS += IDO;
      }
    }
    else {
      IS = 0;
      for (std::size_t J = 1; J < IP; J++) {
        for (std::size_t K = 0; K < L1; K++) {
          IDIJ = IS;
          for (std::size_t I0 = 1; I0 < IDO-1; I0 += 2) {
            std::size_t I1 = I0 + 1;
            C1(I0,K,J) = WA[IDIJ]*CH(I0,K,J)-WA[IDIJ+1]*CH(I1,K,J);
            C1(I1,K,J) = WA[IDIJ]*CH(I1,K,J)+WA[IDIJ+1]*CH(I0,K,J);
            IDIJ += 2;
          }
        }
        IS += IDO;
      }
    }
  }

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_DETAIL_REAL_TO_COMPLEX_BACKWARD_H
