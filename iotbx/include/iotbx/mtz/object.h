#ifndef IOTBX_MTZ_OBJECT_H
#define IOTBX_MTZ_OBJECT_H

#include <cmtzlib.h>
#include <ccp4_errno.h>
#include <cctbx/error.h>
#include <scitbx/array_family/shared.h>

namespace iotbx { namespace mtz {

  namespace af = scitbx::af;

  class crystal;

  class object
  {
    public:
      object()
      :
        ptr(CMtz::MtzMalloc(0, 0), ptr_deleter())
      {
        if (ptr.get() == 0) throw cctbx::error("MtzMalloc failed.");
        ptr->refs_in_memory = true;
      }

      object(af::const_ref<int> const& n_datasets_for_each_crystal)
      :
        ptr(CMtz::MtzMalloc(
            n_datasets_for_each_crystal.size(),
            const_cast<int*>(&*n_datasets_for_each_crystal.begin())),
          ptr_deleter())
      {
        if (ptr.get() == 0) throw cctbx::error("MtzMalloc failed.");
        ptr->refs_in_memory = true;
      }

      object(const char* file_name)
      :
        ptr(CMtz::MtzGet(file_name, true), ptr_deleter())
      {
        if (ptr.get() == 0) {
          throw cctbx::error(std::string("MTZ file read error: ") + file_name);
        }
      }

      int
      n_batches() const { return CMtz::MtzNbat(ptr.get()); }

      int
      n_reflections() const { return CMtz::MtzNref(ptr.get()); }

      int
      space_group_number() const
      {
        return CMtz::MtzSpacegroupNumber(ptr.get());
      }

      af::tiny<double, 2>
      max_min_resolution() const
      {
        if (ptr->nxtal == 0) return af::tiny<double, 2>(-1., -1.);
        float max_resolution;
        float min_resolution;
        CCTBX_ASSERT(
          CMtz::MtzResLimits(ptr.get(), &min_resolution, &max_resolution));
        return af::tiny<double, 2>(max_resolution, min_resolution);
      }

      int
      n_crystals() const { return CMtz::MtzNxtal(ptr.get()); }

      int
      n_active_crystals() const { return CMtz::MtzNumActiveXtal(ptr.get()); }

      af::shared<crystal>
      crystals() const;

    protected:
      boost::shared_ptr<CMtz::MTZ> ptr;

      struct ptr_deleter
      {
        void
        operator()(CMtz::MTZ* ptr) const
        {
          CCTBX_ASSERT(CMtz::MtzFree(ptr));
        }
      };
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_OBJECT_H
