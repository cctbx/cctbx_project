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

    typedef Image::int16_type int16_type;
    typedef Image::int32_type int32_type;
    typedef Image::uint16_type uint16_type;
    typedef Image::uint32_type uint32_type;
    typedef Image::float64_type float32_type;
    typedef Image::float64_type float64_type;
    typedef Image::variant_type variant_type;

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
    Image image() const {
      DXTBX_ASSERT(tiles_.size() == names_.size());
      Image result;
      for (std::size_t i = 0; i < tiles_.size(); ++i) {
        result.push_back(ImageTile(tiles_[i], names_[i].c_str()));
      }
      return result;
    }

  protected:

    std::string filename_;
    std::vector<variant_type> tiles_;
    std::vector<std::string> names_;
  };


  class MultiImageReader {
  public:

    virtual Image image(std::size_t index) const = 0;
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
    typedef typename ImageReaderType::variant_type variant_type;

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
    Image image(std::size_t index) const {
      DXTBX_ASSERT(index < filenames_.size());
      ImageReaderType reader(filenames_[index].c_str());
      return reader.image();
    }

  protected:

    scitbx::af::shared<std::string> filenames_;

  };

}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_IMAGE_READER_H
