// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created (R.W. Grosse-Kunstleve)
 */

/* Reference:
     C.M. Weeks, P.D. Adams, J. Berendzen, A.T. Brunger, E.J. Dodson,
     R.W. Grosse-Kunstleve, T.R. Schneider, G.M. Sheldrick,
     T.C. Terwilliger, M. Turkenburg, I. Uson
     Automatic solution of heavy-atom substructures.
     Methods in Enzymology, in press.
 */

#ifndef CCTBX_TRIPLETS_H
#define CCTBX_TRIPLETS_H

#include <map>
#include <cctbx/array_family/shared.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_asu.h>

#include <cctbx/array_family/simple_io.h> // XXX
#include <cctbx/error.h> // XXX

namespace cctbx { namespace dmtbx {

  typedef af::tiny<std::size_t, 3> triplet_indices;

  template <typename FloatType>
  class triplet_invariants
  {
    public:
      triplet_invariants() {}

      triplet_invariants(sgtbx::SpaceGroupInfo const& SgInfo,
                         af::shared<Miller::Index> miller_indices,
                         af::shared<FloatType> e_values)
      {
        cctbx_assert(miller_indices.size() == e_values.size());
        Miller::index_span miller_index_span(miller_indices);
        typedef std::map<std::size_t, std::size_t> lookup_dict_type;
        lookup_dict_type lookup_dict;
        std::size_t i;
        for(i=0;i<miller_indices.size();i++) {
          lookup_dict[miller_index_span.pack(miller_indices[i])] = i;
        }
        cctbx_assert(lookup_dict.size() == miller_indices.size());
        sgtbx::ReciprocalSpaceASU asu(SgInfo);
        for(i=0;i<miller_indices.size();i++) {
          Miller::Index mh = -miller_indices[i];
          for(std::size_t j=0;j<miller_indices.size();j++) {
            Miller::Index mhmk = mh - miller_indices[j];
            Miller::Index asym_mhmk = Miller::AsymIndex(
              SgInfo.SgOps(), asu, mhmk).HermitianLayout().H();
            if (miller_index_span.is_in_domain(asym_mhmk)) {
              typename lookup_dict_type::const_iterator
              ld_pos = lookup_dict.find(miller_index_span.pack(asym_mhmk));
              if (ld_pos != lookup_dict.end()) {
                cctbx_assert(miller_indices[ld_pos->second] == asym_mhmk);
                std::cout << miller_indices[i].ref()
                   << " " << miller_indices[j].ref()
                   << " " << mhmk.ref()
                   << " " << e_values[i]
                   << " " << e_values[j]
                   << " " << e_values[ld_pos->second]
                   << std::endl;
              }
            }
          }
        }
      }
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_TRIPLETS_H
