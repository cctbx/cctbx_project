#ifndef IOTBX_MTZ_BATCH_H
#define IOTBX_MTZ_BATCH_H

#include <iotbx/mtz/object.h>

namespace iotbx { namespace mtz {

  inline
  bool
  string_is_null_terminated(const char* str, std::size_t str_size)
  {
    for(std::size_t i=0;i<str_size;i++) {
      if (str[i] == '\0') return true;
    }
    return false;
  }

#define IOTBX_MTZ_BATCH_GET_SET(type, name) \
      type \
      name() const { return ptr()->name; } \
       \
      batch& \
      set_##name(type const& value) \
      { \
        ptr()->name = value; \
        return *this; \
      }

#define IOTBX_MTZ_BATCH_GET_SET_ARRAY(type, name, sz) \
      af::shared<type> \
      name() const \
      { \
        CMtz::MTZBAT* batch_ptr = ptr(); \
        return af::shared<type>(batch_ptr->name, batch_ptr->name+sz); \
      } \
       \
      batch& \
      set_##name(af::const_ref<type> const& values) \
      { \
        if (values.size() != sz) { \
          throw cctbx::error("Wrong number of values."); \
        } \
        std::copy(values.begin(), values.end(), ptr()->name); \
        return *this; \
      }

  class batch
  {
    public:
      batch() {}

      batch(object const& mtz_object, int i_batch)
      :
        mtz_object_(mtz_object),
        i_batch_(i_batch)
      {
        CCTBX_ASSERT(i_batch >= 0);
        CCTBX_ASSERT(i_batch < mtz_object.n_batches());
      }

      object
      mtz_object() const { return mtz_object_; }

      int
      i_batch() const { return i_batch_; }

      CMtz::MTZBAT*
      ptr() const
      {
        CCTBX_ASSERT(mtz_object_.n_batches() > i_batch_);
        CMtz::MTZBAT* p = mtz_object_.ptr()->batch;
        for(int i=0;i<i_batch_;i++) {
          if (p == 0) break;
          p = p->next;
        }
        CCTBX_ASSERT(p != 0);
        return p;
      }

      IOTBX_MTZ_BATCH_GET_SET(int, num)

      const char*
      title() const
      {
        CMtz::MTZBAT* batch_ptr = ptr();
        CCTBX_ASSERT(string_is_null_terminated(
          batch_ptr->title, sizeof(batch_ptr->title)));
        return batch_ptr->title;
      }

      batch&
      set_title(const char* value)
      {
        CCTBX_ASSERT(value != 0);
        CMtz::MTZBAT* batch_ptr = ptr();
        strncpy(batch_ptr->title, value, sizeof(batch_ptr->title));
        batch_ptr->title[sizeof(batch_ptr->title)-1] = '\0';
        return *this;
      }

      af::shared<std::string>
      gonlab() const
      {
        CMtz::MTZBAT* batch_ptr = ptr();
        af::shared<std::string> result((af::reserve(3)));
        for(int i=0;i<3;i++) {
          CCTBX_ASSERT(string_is_null_terminated(
            batch_ptr->gonlab[i], sizeof(batch_ptr->gonlab)/3));
          result.push_back(batch_ptr->gonlab[i]);
        }
        return result;
      }

      batch&
      set_gonlab(af::const_ref<std::string> const& values)
      {
        CCTBX_ASSERT(values.size() == 3);
        CMtz::MTZBAT* batch_ptr = ptr();
        const std::size_t n = sizeof(batch_ptr->gonlab)/3;
        for(int i=0;i<3;i++) {
          strncpy(batch_ptr->gonlab[i], values[i].c_str(), n-1);
          batch_ptr->gonlab[i][n-2] = '\0';
          if (strchr(batch_ptr->gonlab[i], ' ') != 0) {
            throw cctbx::error(
              "MTZ batch \"gonlab\" values must not contain spaces.");
          }
        }
        return *this;
      }

      IOTBX_MTZ_BATCH_GET_SET(int, iortyp)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(int, lbcell, 6)
      IOTBX_MTZ_BATCH_GET_SET(int, misflg)
      IOTBX_MTZ_BATCH_GET_SET(int, jumpax)
      IOTBX_MTZ_BATCH_GET_SET(int, ncryst)
      IOTBX_MTZ_BATCH_GET_SET(int, lcrflg)
      IOTBX_MTZ_BATCH_GET_SET(int, ldtype)
      IOTBX_MTZ_BATCH_GET_SET(int, jsaxs)
      IOTBX_MTZ_BATCH_GET_SET(int, nbscal)
      IOTBX_MTZ_BATCH_GET_SET(int, ngonax)
      IOTBX_MTZ_BATCH_GET_SET(int, lbmflg)

      int
      ndet() const { return ptr()->ndet; }

      batch&
      set_ndet(int const& value)
      {
        CCTBX_ASSERT(value >= 0 && value <= 2);
        ptr()->ndet = value;
        return *this;
      }

      IOTBX_MTZ_BATCH_GET_SET(int, nbsetid)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, cell, 6)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, umat, 9)

      af::shared<float>
      phixyz() const
      {
        af::shared<float> result((af::reserve(6)));
        CMtz::MTZBAT* batch_ptr = ptr();
        for(int i=0;i<2;i++)
        for(int j=0;j<3;j++) {
          result.push_back(batch_ptr->phixyz[i][j]);
        }
        return result;
      }

      batch&
      set_phixyz(af::const_ref<float> const& values)
      {
        CCTBX_ASSERT(values.size() == 6);
        CMtz::MTZBAT* batch_ptr = ptr();
        int iv = 0;
        for(int i=0;i<2;i++)
        for(int j=0;j<3;j++,iv++) {
          batch_ptr->phixyz[i][j] = values[iv];
        }
        return *this;
      }

      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, crydat, 12)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, datum, 3)
      IOTBX_MTZ_BATCH_GET_SET(float, phistt)
      IOTBX_MTZ_BATCH_GET_SET(float, phiend)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, scanax, 3)
      IOTBX_MTZ_BATCH_GET_SET(float, time1)
      IOTBX_MTZ_BATCH_GET_SET(float, time2)
      IOTBX_MTZ_BATCH_GET_SET(float, bscale)
      IOTBX_MTZ_BATCH_GET_SET(float, bbfac)
      IOTBX_MTZ_BATCH_GET_SET(float, sdbscale)
      IOTBX_MTZ_BATCH_GET_SET(float, sdbfac)
      IOTBX_MTZ_BATCH_GET_SET(float, phirange)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, e1, 3)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, e2, 3)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, e3, 3)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, source, 3)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, so, 3)
      IOTBX_MTZ_BATCH_GET_SET(float, alambd)
      IOTBX_MTZ_BATCH_GET_SET(float, delamb)
      IOTBX_MTZ_BATCH_GET_SET(float, delcor)
      IOTBX_MTZ_BATCH_GET_SET(float, divhd)
      IOTBX_MTZ_BATCH_GET_SET(float, divvd)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, dx, 2)
      IOTBX_MTZ_BATCH_GET_SET_ARRAY(float, theta, 2)

      af::shared<float>
      detlm() const
      {
        af::shared<float> result((af::reserve(8)));
        CMtz::MTZBAT* batch_ptr = ptr();
        for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
        for(int k=0;k<2;k++) {
          result.push_back(batch_ptr->detlm[i][j][k]);
        }
        return result;
      }

      batch&
      set_detlm(af::const_ref<float> const& values)
      {
        CCTBX_ASSERT(values.size() == 8);
        CMtz::MTZBAT* batch_ptr = ptr();
        int iv = 0;
        for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
        for(int k=0;k<2;k++, iv++) {
          batch_ptr->detlm[i][j][k] = values[iv];
        }
        return *this;
      }

    protected:
      object mtz_object_;
      int i_batch_;
  };

  inline
  af::shared<batch>
  object::batches() const
  {
    af::shared<batch> result((af::reserve(n_batches())));
    for(int i_batch=0;i_batch<n_batches();i_batch++) {
      result.push_back(batch(*this, i_batch));
    }
    return result;
  }

  inline
  batch
  object::add_batch()
  {
    CMtz::MTZBAT* p = ptr()->batch;
    CMtz::MTZBAT* p_tail = p;
    int max_batch_number = 0;
    int i_batch = 0;
    for(;;i_batch++) {
      if (p == 0) break;
      max_batch_number = std::max(max_batch_number, p->num);
      p_tail = p;
      p = p->next;
    }
    std::vector<float> buf(
      static_cast<std::size_t>(NBATCHINTEGERS + NBATCHREALS),
      static_cast<float>(0));
    CCTBX_ASSERT(sizeof(float) == sizeof(int));
    CCTBX_ASSERT(CMtz::ccp4_lwbat(
      ptr(), 0, max_batch_number+1, &*buf.begin(), "") == 1);
    p = (p_tail == 0 ? ptr()->batch : p_tail->next);
    CCTBX_ASSERT(p != 0);
    CCTBX_ASSERT(p->next == 0);
    CCTBX_ASSERT(p->num == max_batch_number+1);
    return batch(*this, i_batch);
  }

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_BATCH_H
