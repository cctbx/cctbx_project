/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Renamed miller_ref_asu.cpp -> reciprocal_space_ref_asu.cpp
     2002 Jul: split of miller_asu.cpp (R.W. Grosse-Kunstleve)
     2001 Oct: Redesign: AsymIndex (rwgk)
     2001 Sep: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#include <cctbx/sgtbx/reciprocal_space_reference_asu.h>

namespace cctbx { namespace sgtbx { namespace reciprocal_space {

  namespace reference_asu_definitions {
  namespace {

    struct asu_1b : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_1b;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[2]>0 || (h[2]==0 && (h[0]>0 || (h[0]==0 && h[1]>=0))));
      }
      virtual const char* as_string() const {
        return "l>0 or (l==0 and (h>0 or (h==0 and k>=0)))";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(-1, -1, 0);
        return result;
      }
    };

    struct asu_2_m : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_2_m;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[1]>=0 && (h[2]>0 || (h[2]==0 && h[0]>=0)));
      }
      virtual const char* as_string() const {
        return "k>=0 and (l>0 or (l=0 and h>=0))";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(-1, 0, 0);
        return result;
      }
    };

    struct asu_mmm : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_mmm;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=0 && h[1]>=0 && h[2]>=0);
      }
      virtual const char* as_string() const {
        return "h>=0 and k>=0 and l>=0";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    struct asu_4_m : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_4_m;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[2]>=0 && ((h[0]>=0 && h[1]>0) || (h[0]==0 && h[1]==0)));
      }
      virtual const char* as_string() const {
        return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    struct asu_4_mmm : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_4_mmm;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=h[1] && h[1]>=0 && h[2]>=0);
      }
      virtual const char* as_string() const {
        return "h>=k and k>=0 and l>=0";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    struct asu_3b : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_3b;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=0 && h[1]>0) || (h[0]==0 && h[1]==0 && h[2]>=0);
      }
      virtual const char* as_string() const {
        return "(h>=0 and k>0) or (h=0 and k=0 and l>=0)";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, -1);
        return result;
      }
    };

    struct asu_3b1m : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_3b1m;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=h[1] && h[1]>=0 && (h[1]>0 || h[2]>=0));
      }
      virtual const char* as_string() const {
        return "h>=k and k>=0 and (k>0 or l>=0)";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, -1);
        return result;
      }
    };

    struct asu_3bm1 : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_3bm1;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=h[1] && h[1]>=0 && (h[0]>h[1] || h[2]>=0));
      }
      virtual const char* as_string() const {
        return "h>=k and k>=0 and (h>k or l>=0)";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, -1);
        return result;
      }
    };

    struct asu_6_m : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_6_m;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[2]>=0 && ((h[0]>=0 && h[1]>0) || (h[0]==0 && h[1]==0)));
      }
      virtual const char* as_string() const {
        return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    struct asu_6_mmm : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_6_mmm;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=h[1] && h[1]>=0 && h[2]>=0);
      }
      virtual const char* as_string() const {
        return "h>=k and k>=0 and l>=0";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    struct asu_m3b : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_m3b;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[0]>=0 && (   (h[2]>=h[0] && h[1]>h[0])
                            || (h[2]==h[0] && h[1]==h[0])));
      }
      virtual const char* as_string() const {
        return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    struct asu_m3bm : reference_asu
    {
      virtual matrix_group::code laue_group() const {
        return matrix_group::code_m3bm;
      }
      virtual bool is_inside(miller::index<> const& h) const {
        return (h[1]>=h[2] && h[2]>=h[0] && h[0]>=0);
      }
      virtual const char* as_string() const {
        return "k>=l and l>=h and h>=0";
      }
      virtual af::int3 const& cut_parameters() const {
        static const af::int3 result(0, 0, 0);
        return result;
      }
    };

    static const asu_1b    asu_instance_1b    = asu_1b();
    static const asu_2_m   asu_instance_2_m   = asu_2_m();
    static const asu_mmm   asu_instance_mmm   = asu_mmm();
    static const asu_4_m   asu_instance_4_m   = asu_4_m();
    static const asu_4_mmm asu_instance_4_mmm = asu_4_mmm();
    static const asu_3b    asu_instance_3b    = asu_3b();
    static const asu_3b1m  asu_instance_3b1m  = asu_3b1m();
    static const asu_3bm1  asu_instance_3bm1  = asu_3bm1();
    static const asu_6_m   asu_instance_6_m   = asu_6_m();
    static const asu_6_mmm asu_instance_6_mmm = asu_6_mmm();
    static const asu_m3b   asu_instance_m3b   = asu_m3b();
    static const asu_m3bm  asu_instance_m3bm  = asu_m3bm();

    static const reference_asu* table[] = {
      &asu_instance_1b,
      &asu_instance_2_m,
      &asu_instance_mmm,
      &asu_instance_4_m, &asu_instance_4_mmm,
      &asu_instance_3b, &asu_instance_3b1m, &asu_instance_3bm1,
      &asu_instance_6_m, &asu_instance_6_mmm,
      &asu_instance_m3b, &asu_instance_m3bm,
      0,
    };

  } // namespace <anonymous>
  } // namespace reference_asu_definitions

  bool is_in_reference_asu_1b(miller::index<> const& h)
  {
    return reference_asu_definitions::asu_instance_1b.is_inside(h);
  }

  const reference_asu*
  lookup_reference_asu(matrix_group::code const& group)
  {
    using namespace matrix_group;
    code laue_group = group.laue_group_type();
    if (laue_group == code_3bm) {
      if (   group == code_312
          || group == code_31m
          || group == code_3b1m) laue_group = code_3b1m;
      else                       laue_group = code_3bm1;
    }
    for(std::size_t i=0;; i++) {
      const reference_asu* ref_asu = reference_asu_definitions::table[i];
      if (ref_asu == 0) throw CCTBX_INTERNAL_ERROR();
      if (ref_asu->laue_group() == laue_group) return ref_asu;
    }
  }

}}} // namespace cctbx::sgtbx::reciprocal_space
