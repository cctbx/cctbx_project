#ifndef SCITBX_LINE_SEARCH_MORE_THUENTE_1994_H
#define SCITBX_LINE_SEARCH_MORE_THUENTE_1994_H

#include <scitbx/line_search/more_thuente_1994_raw.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace line_search {

  template <typename FloatType=double>
  class more_thuente_1994
  {
    protected:
      mcsrch<> raw;

    public:
      FloatType xtol;
      FloatType ftol;
      FloatType gtol;
      FloatType stpmin;
      FloatType stpmax;
      unsigned maxfev;
      int info_code;
      const char* info_meaning;
      unsigned nfev;
      af::shared<FloatType> search_direction;
      FloatType stp;

      more_thuente_1994()
      :
        xtol(1.e-16),
        ftol(1.e-4),
        gtol(0.9),
        stpmin(1.e-20),
        stpmax(1.e20),
        maxfev(20),
        info_code(0),
        info_meaning(0),
        nfev(0),
        stp(0)
      {
        raw.info_meaning = 0;
      }

      int
      start(
        af::ref<FloatType> const& x,
        FloatType const& functional,
        af::const_ref<FloatType> const& gradients,
        af::shared<FloatType> const& search_direction,
        FloatType const& initial_estimate_of_satisfactory_step_length)
      {
        SCITBX_ASSERT(gradients.size() == x.size());
        SCITBX_ASSERT(search_direction.size() == x.size());
        SCITBX_ASSERT(initial_estimate_of_satisfactory_step_length > 0);
        info_code = 0;
        info_meaning = 0;
        nfev = 0;
        this->search_direction = search_direction;
        stp = initial_estimate_of_satisfactory_step_length;
        raw.run(
          gtol,
          stpmin,
          stpmax,
          x.size(),
          x.begin(),
          functional,
          gradients.begin(),
          search_direction.begin(),
          stp,
          ftol,
          xtol,
          maxfev,
          info_code,
          nfev);
        info_meaning = raw.info_meaning;
        SCITBX_ASSERT(info_code == -1 || (info_code >= 1 && info_code <= 6));
        return info_code;
      }

      int
      next(
        af::ref<FloatType> const& x,
        FloatType const& functional,
        af::const_ref<FloatType> const& gradients)
      {
        SCITBX_ASSERT(info_code == -1);
        SCITBX_ASSERT(gradients.size() == x.size());
        SCITBX_ASSERT(search_direction.size() == x.size());
        raw.run(
          gtol,
          stpmin,
          stpmax,
          x.size(),
          x.begin(),
          functional,
          gradients.begin(),
          search_direction.begin(),
          stp,
          ftol,
          xtol,
          maxfev,
          info_code,
          nfev);
        info_meaning = raw.info_meaning;
        SCITBX_ASSERT(info_code == -1 || (info_code >= 1 && info_code <= 6));
        if (info_code != -1) raw.free_workspace();
        return info_code;
      }
  };

}} // namespace scitbx::line_search

#endif // SCITBX_LINE_SEARCH_MORE_THUENTE_1994_H
