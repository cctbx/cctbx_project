// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: split of miller_asu.cpp (R.W. Grosse-Kunstleve)
     2001 Oct 31: Redesign: AsymIndex (rwgk)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#include <cctbx/sgtbx/miller_ref_asu.h>

namespace cctbx { namespace sgtbx {

  namespace detail {

    class ReferenceReciprocalSpaceASU_1b
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_1b;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[2]>0 || (H[2]==0 && (H[0]>0 || (H[0]==0 && H[1]>=0))));
        }
        virtual const char* representation() const {
          return "l>0 or (l==0 and (h>0 or (h==0 and k>=0)))";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(-1, -1, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_2_m
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_2_m;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[1]>=0 && (H[2]>0 || (H[2]==0 && H[0]>=0)));
        }
        virtual const char* representation() const {
          return "k>=0 and (l>0 or (l=0 and h>=0))";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(-1, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_mmm
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_mmm;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=0 && H[1]>=0 && H[2]>=0);
        }
        virtual const char* representation() const {
          return "h>=0 and k>=0 and l>=0";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_4_m
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_4_m;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[2]>=0 && ((H[0]>=0 && H[1]>0) || (H[0]==0 && H[1]==0)));
        }
        virtual const char* representation() const {
          return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_4_mmm
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_4_mmm;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=H[1] && H[1]>=0 && H[2]>=0);
        }
        virtual const char* representation() const {
          return "h>=k and k>=0 and l>=0";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_3b
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_3b;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=0 && H[1]>0) || (H[0]==0 && H[1]==0 && H[2]>=0);
        }
        virtual const char* representation() const {
          return "(h>=0 and k>0) or (h=0 and k=0 and l>=0)";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, -1);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_3b1m
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_3b1m;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=H[1] && H[1]>=0 && (H[1]>0 || H[2]>=0));
        }
        virtual const char* representation() const {
          return "h>=k and k>=0 and (k>0 or l>=0)";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, -1);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_3bm1
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_3bm1;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=H[1] && H[1]>=0 && (H[0]>H[1] || H[2]>=0));
        }
        virtual const char* representation() const {
          return "h>=k and k>=0 and (h>k or l>=0)";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, -1);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_6_m
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_6_m;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[2]>=0 && ((H[0]>=0 && H[1]>0) || (H[0]==0 && H[1]==0)));
        }
        virtual const char* representation() const {
          return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_6_mmm
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_6_mmm;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=H[1] && H[1]>=0 && H[2]>=0);
        }
        virtual const char* representation() const {
          return "h>=k and k>=0 and l>=0";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_m3b
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_m3b;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[0]>=0 && (   (H[2]>=H[0] && H[1]>H[0])
                              || (H[2]==H[0] && H[1]==H[0])));
        }
        virtual const char* representation() const {
          return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };
    class ReferenceReciprocalSpaceASU_m3bm
      : public ReferenceReciprocalSpaceASU {
      public:
        virtual tables::MatrixGroup::Code LaueGroupCode() const {
          return tables::MatrixGroup::MGC_m3bm;
        }
        virtual bool isInASU(const Miller::Index& H) const {
          return (H[1]>=H[2] && H[2]>=H[0] && H[0]>=0);
        }
        virtual const char* representation() const {
          return "k>=l and l>=h and h>=0";
        }
        virtual const af::int3& getCutParameters() const {
          static const af::int3 result(0, 0, 0);
          return result;
        }
    };

    static const ReferenceReciprocalSpaceASU_1b    ReferenceASU_1b =
                 ReferenceReciprocalSpaceASU_1b();
    static const ReferenceReciprocalSpaceASU_2_m   ReferenceASU_2_m =
                 ReferenceReciprocalSpaceASU_2_m();
    static const ReferenceReciprocalSpaceASU_mmm   ReferenceASU_mmm =
                 ReferenceReciprocalSpaceASU_mmm();
    static const ReferenceReciprocalSpaceASU_4_m   ReferenceASU_4_m =
                 ReferenceReciprocalSpaceASU_4_m();
    static const ReferenceReciprocalSpaceASU_4_mmm ReferenceASU_4_mmm =
                 ReferenceReciprocalSpaceASU_4_mmm();
    static const ReferenceReciprocalSpaceASU_3b    ReferenceASU_3b =
                 ReferenceReciprocalSpaceASU_3b();
    static const ReferenceReciprocalSpaceASU_3b1m  ReferenceASU_3b1m =
                 ReferenceReciprocalSpaceASU_3b1m();
    static const ReferenceReciprocalSpaceASU_3bm1  ReferenceASU_3bm1 =
                 ReferenceReciprocalSpaceASU_3bm1();
    static const ReferenceReciprocalSpaceASU_6_m   ReferenceASU_6_m =
                 ReferenceReciprocalSpaceASU_6_m();
    static const ReferenceReciprocalSpaceASU_6_mmm ReferenceASU_6_mmm =
                 ReferenceReciprocalSpaceASU_6_mmm();
    static const ReferenceReciprocalSpaceASU_m3b   ReferenceASU_m3b =
                 ReferenceReciprocalSpaceASU_m3b();
    static const ReferenceReciprocalSpaceASU_m3bm  ReferenceASU_m3bm =
                 ReferenceReciprocalSpaceASU_m3bm();

    static const ReferenceReciprocalSpaceASU*
    TableReferenceReciprocalSpaceASU[] = {
      &ReferenceASU_1b,
      &ReferenceASU_2_m,
      &ReferenceASU_mmm,
      &ReferenceASU_4_m, &ReferenceASU_4_mmm,
      &ReferenceASU_3b, &ReferenceASU_3b1m, &ReferenceASU_3bm1,
      &ReferenceASU_6_m, &ReferenceASU_6_mmm,
      &ReferenceASU_m3b, &ReferenceASU_m3bm,
      0,
    };

  } // namespace detail

  bool isInReferenceReciprocalSpaceASU_1b(Miller::Index const& h)
  {
    return detail::ReferenceASU_1b.isInASU(h);
  }

  const ReferenceReciprocalSpaceASU*
  LookupReferenceReciprocalSpaceASU(tables::MatrixGroup::Code group_code)
  {
    using namespace tables::MatrixGroup;
    group_code = group_code.LaueGroupType();
    if (group_code == MGC_3bm) {
      if (   group_code == MGC_312
          || group_code == MGC_31m
          || group_code == MGC_3b1m) group_code = MGC_3b1m;
      else                           group_code = MGC_3bm1;
    }
    for(std::size_t i=0;; i++) {
      const ReferenceReciprocalSpaceASU*
      ref_asu = detail::TableReferenceReciprocalSpaceASU[i];
      if (ref_asu == 0) throw cctbx_internal_error();
      if (ref_asu->LaueGroupCode() == group_code) return ref_asu;
    }
  }

}} // namespace cctbx::sgtbx
