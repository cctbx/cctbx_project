/*
 * image_reader.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_IMAGE_READER_H
#define DXTBX_FORMAT_IMAGE_READER_H

#include <vector>
#include <boost/variant.hpp>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace dxtbx { namespace format {

  /**
   * Base class for reading images
   */
  class ImageReader {
  public:

    typedef scitbx::af::c_grid<2> accessor_type;

    typedef scitbx::af::versa< short,          accessor_type > int16_type;
    typedef scitbx::af::versa< int,            accessor_type > int32_type;
    typedef scitbx::af::versa< unsigned short, accessor_type > uint16_type;
    typedef scitbx::af::versa< unsigned int,   accessor_type > uint32_type;
    typedef scitbx::af::versa< float,          accessor_type > float32_type;
    typedef scitbx::af::versa< double,         accessor_type > float64_type;

    typedef boost::variant<
      int16_type,
      int32_type,
      uint16_type,
      uint32_type,
      float32_type,
      float64_type
    > variant_type;

    /**
     * Initialise with the filename
     */
    ImageReader(const char *filename)
      : filename_(filename) {}

    /**
     * Return the filename
     */
    std::string filename() const {
      return filename_;
    }

    /**
     * Return the image
     */
    ImageBuffer image() const {
      return buffer_;
    }

  protected:

    /* /** */
    /*  * Helper function to convert a list of tile variants to an image buffer */
    /*  *1/ */
    /* ImageBuffer tile_list_to_image_buffer( */
    /*     scitbx::af::const_ref<variant_type> &tiles, */
    /*     scitbx::af::const_ref<std::string> &names) { */
    /*   if (boost::apply_visitor(IsDoubleVisitor(), tiles)) { */
    /*     return ImageBuffer(boost::apply_visitor(ImageConstructor<double>(names), tiles); */
    /*   } */
    /*   return ImageBuffer(boost::apply_visitor(ImageConstructor<int>(names), tiles); */
    /* } */

    std::string filename_;
    ImageBuffer buffer_;
  };


  class MultiImageReader {
  public:

    virtual ImageBuffer image(std::size_t index) const = 0;
    virtual std::size_t size() const = 0;

  };

  /**
   * Class to read from a list of images
   */
  template <typename ImageReaderType>
  class ImageListReader : public MultiImageReader {
  public:

    typedef ImageReaderType image_reader_type;
    typedef typename ImageReaderType::int16_type int16_type;
    typedef typename ImageReaderType::int32_type int32_type;
    typedef typename ImageReaderType::uint16_type uint16_type;
    typedef typename ImageReaderType::uint32_type uint32_type;
    typedef typename ImageReaderType::float64_type float32_type;
    typedef typename ImageReaderType::float64_type float64_type;
    // FIXME typedef typename ImageReaderType::variant_type variant_type;

    /**
     * Initialise with the filename
     */
    ImageListReader(const scitbx::af::const_ref<std::string> &filenames)
      : filenames_(filenames.begin(), filenames.end()) {}

    /**
     * Return the filename
     */
    scitbx::af::shared<std::string> filenames() const {
      return filenames_;
    }

    /**
     * Return the number of images
     */
    std::size_t size() const {
      return filenames_.size();
    }

    /**
     * Return the image
     */
    ImageBuffer image(std::size_t index) const {
      DXTBX_ASSERT(index < filenames_.size());
      ImageReaderType reader(filenames_[index].c_str());
      return reader.image();
    }

  protected:

    scitbx::af::shared<std::string> filenames_;

  };

}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_IMAGE_READER_H
