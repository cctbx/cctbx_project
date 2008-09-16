#ifndef IOTBX_MTZ_BATCH_H
#define IOTBX_MTZ_BATCH_H

#include <iotbx/mtz/dataset.h>
#include <iotbx/error.h>

namespace iotbx { namespace mtz {

  //! Helper function to facilitate friendly assertions.
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

  //! Safe access to CMtz::MTZBAT* owned by an iotbx::mtz::object .
  /*! A batch object contains (owns) an iotbx::mtz::object and
      an integer index into the list of batches owned by the
      object.

      All members of the CMtz::MTZBAT struct are accessible through
      both a read-only interface, and a write interface, e.g.
      title() and set_title(). The names of the access methods
      are determined strictly by the names in the CMtz::MTZBAT struct.
      The access methods are not documented individually.
      Please refer to the injected show() method in
          $IOTBX_DIST/iotbx/mtz/__init__.py
      to see all names.
   */
  class batch
  {
    public:
      //! Not available from Python.
      batch() {}

      /*! \brief Initialization given a mtz_object and an integer
          index into the list of batches owned by the mtz_object.
       */
      /*! An exception is thrown if i_batch is out of range.
       */
      batch(object const& mtz_object, int i_batch)
      :
        mtz_object_(mtz_object),
        i_batch_(i_batch)
      {
        IOTBX_ASSERT(i_batch >= 0);
        IOTBX_ASSERT(i_batch < mtz_object.n_batches());
      }

      //! The contained iotbx::mtz::object instance.
      object
      mtz_object() const { return mtz_object_; }

      //! Integer index into the list of batches owned by mtz_object() .
      int
      i_batch() const { return i_batch_; }

      //! Raw C pointer. Not available from Python.
      /*! The pointer is obtained by traversing the linked list
          of batches.
          An exception is thrown if i_batch() is out of range.
       */
      CMtz::MTZBAT*
      ptr() const
      {
        IOTBX_ASSERT(mtz_object_.n_batches() > i_batch_);
        CMtz::MTZBAT* p = mtz_object_.ptr()->batch;
        for(int i=0;i<i_batch_;i++) {
          if (p == 0) break;
          p = p->next;
        }
        IOTBX_ASSERT(p != 0);
        return p;
      }

      IOTBX_MTZ_BATCH_GET_SET(int, num)

      const char*
      title() const
      {
        CMtz::MTZBAT* batch_ptr = ptr();
        IOTBX_ASSERT(string_is_null_terminated(
          batch_ptr->title, sizeof(batch_ptr->title)));
        return batch_ptr->title;
      }

      batch&
      set_title(const char* value)
      {
        IOTBX_ASSERT(value != 0);
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
          IOTBX_ASSERT(string_is_null_terminated(
            batch_ptr->gonlab[i], sizeof(batch_ptr->gonlab)/3));
          result.push_back(batch_ptr->gonlab[i]);
        }
        return result;
      }

      batch&
      set_gonlab(af::const_ref<std::string> const& values)
      {
        IOTBX_ASSERT(values.size() == 3);
        CMtz::MTZBAT* batch_ptr = ptr();
        const std::size_t n = sizeof(batch_ptr->gonlab)/3;
        for(int i=0;i<3;i++) {
          strncpy(batch_ptr->gonlab[i], values[i].c_str(), n-1);
          batch_ptr->gonlab[i][n-2] = '\0';
          if (std::strchr(batch_ptr->gonlab[i], ' ') != 0) {
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
        IOTBX_ASSERT(value >= 0 && value <= 2);
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
        IOTBX_ASSERT(values.size() == 6);
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
        IOTBX_ASSERT(values.size() == 8);
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

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_BATCH_H
