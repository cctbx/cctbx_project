#include <iotbx/mtz/dataset.h>

namespace iotbx { namespace mtz {

  crystal::crystal(object const& mtz_object, int i_crystal)
  :
    mtz_object_(mtz_object),
    i_crystal_(i_crystal)
  {
    CCTBX_ASSERT(i_crystal >= 0);
    CCTBX_ASSERT(i_crystal < mtz_object.n_crystals());
  }

  CMtz::MTZXTAL*
  crystal::ptr() const
  {
    CCTBX_ASSERT(mtz_object_.n_crystals() > i_crystal_);
    return CMtz::MtzIxtal(mtz_object_.ptr(), i_crystal_);
  }

  crystal&
  crystal::set_id(int id)
  {
    if (ptr()->xtalid != id) {
      CMtz::MTZ* p = mtz_object().ptr();
      CCTBX_ASSERT(p->refs_in_memory);
      for(int i=0;i<p->nxtal;i++) {
        CCTBX_ASSERT(p->xtal[i]->xtalid != id);
      }
      ptr()->xtalid = id;
    }
    return *this;
  }

  af::small<double, 6>
  crystal::unit_cell_parameters() const
  {
    af::small<double, 6> result;
    float* cell = ptr()->cell;
    for(std::size_t i=0;i<6;i++) result.push_back(cell[i]);
    return result;
  }

  af::shared<dataset>
  crystal::datasets() const
  {
    af::shared<dataset> result((af::reserve(n_datasets())));
    for(int i_dataset=0;i_dataset<n_datasets();i_dataset++) {
      result.push_back(dataset(*this, i_dataset));
    }
    return result;
  }

  crystal&
  crystal::set_unit_cell_parameters(af::small<double, 6> const& parameters)
  {
    float* cell = ptr()->cell;
    for(std::size_t i=0;i<6;i++) cell[i] = static_cast<float>(parameters[i]);
    return *this;
  }

  dataset
  crystal::add_dataset(
    const char *name,
    double wavelength)
  {
    int i_dataset = n_datasets();
    CMtz::MTZSET* dataset_ptr = CMtz::MtzAddDataset(
      mtz_object().ptr(), ptr(), name, static_cast<float>(wavelength));
    CCTBX_ASSERT(dataset_ptr != 0);
    CCTBX_ASSERT(n_datasets() == i_dataset+1);
    dataset result(*this, i_dataset);
    CCTBX_ASSERT(result.ptr() == dataset_ptr);
    return result;
  }

}} // namespace iotbx::mtz
