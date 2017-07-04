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

#include <vector>
#include <boost/variant.hpp>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace dxtbx { namespace format {


  /**
   * A class to represent an image tile
   */
  class ImageTile {
  public:

    // Variant data typedefs
    typedef scitbx::af::versa< short,          scitbx::af::c_grid<2> > int16_type;
    typedef scitbx::af::versa< int,            scitbx::af::c_grid<2> > int32_type;
    typedef scitbx::af::versa< unsigned short, scitbx::af::c_grid<2> > uint16_type;
    typedef scitbx::af::versa< unsigned int,   scitbx::af::c_grid<2> > uint32_type;
    typedef scitbx::af::versa< float,          scitbx::af::c_grid<2> > float32_type;
    typedef scitbx::af::versa< double,         scitbx::af::c_grid<2> > float64_type;

    typedef int32_type int_array_type;
    typedef float64_type double_array_type;

    typedef boost::variant <
      int16_type,
      int32_type,
      uint16_type,
      uint32_type,
      float32_type,
      float64_type
    > variant_type;


    /**
     * A visitor class to convert from/to different types.
     * Data is copied except when the to/from types are the same
     */
    template <typename ArrayType>
    class ConverterVisitor : public boost::static_visitor<ArrayType> {
    public:

      ArrayType operator()(ArrayType &v) const {
        return v;
      }

      template <typename OtherArrayType>
      ArrayType operator()(OtherArrayType &v) const {
        ArrayType result(v.accessor());
        std::copy(v.begin(), v.end(), result.begin());
        return result;
      }

    };

    /**
     * Is the data a double type
     */
    class IsDoubleVisitor : public boost::static_visitor<bool> {
    public:

      bool operator()(float32_type &v) const {
        return true;
      }

      bool operator()(float64_type &v) const {
        return true;
      }

      template <typename OtherArrayType>
      bool operator()(OtherArrayType &v) const {
        return false;
      }

    };


    /**
     * Initialize the class
     */
    ImageTile(variant_type data, const char *name)
      : data_(data),
        name_(name) {}

    /**
     * Is the data integer
     */
    bool is_int() const {
      return !boost::apply_visitor(IsDoubleVisitor(), data_);
    }

    /**
     * Is the data double
     */
    bool is_double() const {
      return boost::apply_visitor(IsDoubleVisitor(), data_);
    }

    /**
     * Get the data as integer
     */
    int_array_type as_int() const {
      return boost::apply_visitor(ConverterVisitor<int_array_type>(), data_);
    }

    /**
     * Get the data as double
     */
    double_array_type as_double() const {
      return boost::apply_visitor(ConverterVisitor<double_array_type>(), data_);
    }

    /**
     * Get the tile name
     */
    std::string name() {
      return name_;
    }

  protected:

    variant_type data_;
    std::string name_;
  };


  /**
   * A class to represent a multi-tile image
   */
  class Image {
  public:
    typedef ImageTile::int16_type int16_type;
    typedef ImageTile::int32_type int32_type;
    typedef ImageTile::uint16_type uint16_type;
    typedef ImageTile::uint32_type uint32_type;
    typedef ImageTile::float64_type float32_type;
    typedef ImageTile::float64_type float64_type;
    typedef ImageTile::variant_type variant_type;

    /**
     * Add a tile
     */
    void push_back(const variant_type &tile, const char *name) {
      tiles_.push_back(tile);
      names_.push_back(name);
    }

    /**
     * Get the tile names
     */
    scitbx::af::shared< std::string > tile_names() const {
      return scitbx::af::shared< std::string >(
          &names_[0],
          &names_[0] + names_.size());
    }

    /**
     * Get an image tile
     */
    ImageTile tile(std::size_t index) {
      DXTBX_ASSERT(index < n_tiles());
      return ImageTile(tiles_[index], names_[index].c_str());
    }

    /**
     * Get the number of tiles
     */
    std::size_t n_tiles() const {
      return tiles_.size();
    }

  protected:

    std::vector< variant_type > tiles_;
    std::vector< std::string > names_;
  };


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
        result.push_back(tiles_[i], names_[i].c_str());
      }
      return result;
    }

  protected:

    std::string filename_;
    std::vector<variant_type> tiles_;
    std::vector<std::string> names_;
  };

  
  /**
   * Class to read from a list of images
   */
  template <typename ImageReaderType>
  class ImageListReader {
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

#endif // DXTBX_FORMAT_IMAGE_H
