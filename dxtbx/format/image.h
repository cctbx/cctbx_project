/*
 * image.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_IMAGE_H
#define DXTBX_FORMAT_IMAGE_H

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace dxtbx { namespace format {

  /**
   * Base class for reading images
   */
  class ImageReader {
  public:

    typedef scitbx::af::tiny<std::size_t,2> size_type;

    ImageReader(const char *filename)
      : filename_(filename) {}

    /**
     * Get the data as integer
     */
    virtual
    scitbx::af::versa< int, scitbx::af::c_grid<2> > as_int() const {
      throw std::runtime_error("Override");
    }

    /**
     * Get the data as double
     */
    virtual
    scitbx::af::versa< double, scitbx::af::c_grid<2> > as_double() const {
      throw std::runtime_error("Override");
    }

    /**
     * Return the image filename
     */
    std::string filename() const {
      return filename_;
    }

    /**
     * Return the size of the image
     */
    size_type size() const {
      return size_;
    }

    /**
     * Return is an int
     */
    bool is_int() const {
      return
        type_ == "int16" ||
        type_ == "int32" ||
        type_ == "uint16" ||
        type_ == "uint32" ;
    }

    /**
     * Return if is a float
     */
    bool is_float() const {
      return type_ == "float32" || type_ == "float64";
    }

    /**
     * Return the data type
     */
    std::string type() const {
      return type_;
    }

  protected:

    std::string filename_;
    size_type size_;
    std::string type_;
  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_IMAGE_H
