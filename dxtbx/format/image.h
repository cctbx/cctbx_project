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

#include <dxtbx/error.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace dxtbx { namespace format {


  /**
   * An image tile containing data from a single detector panel
   */
  template <typename T>
  class ImageTile {
  public:

    typedef scitbx::af::c_grid<2> accessor_type;
    typedef scitbx::af::versa< T, accessor_type > array_type;

    /**
     * Initialize the class
     */
    ImageTile(array_type data)
      : data_(data),
        name_("") {}

    /**
     * Initialize the class
     */
    ImageTile(array_type data, const char *name)
      : data_(data),
        name_(name) {}

    /**
     * Get the image data
     */
    array_type data() const {
      return data_;
    }

    /**
     * Get the tile name
     */
    std::string name() const {
      return name_;
    }

    /**
     * Is the array empty
     */
    bool empty() const {
      return data_.empty();
    }

    /**
     * Get the accessor
     */
    accessor_type accessor() const {
      return data_.accessor();
    }

  protected:

    array_type data_;
    std::string name_;

  };


  /**
   * An image containing data from all detector panels
   */
  template <typename T>
  class Image {
  public:

    typedef T data_type;
    typedef ImageTile<T> tile_type;
    typedef typename tile_type::array_type array_type;
    typedef typename scitbx::af::shared<tile_type> tile_array_type;
    typedef typename tile_array_type::iterator iterator;

    /**
     * Construct empty
     */
    Image() {

    }

    /**
     * Construct with a single tile
     */
    Image(const ImageTile<T> &tile) {
      tiles_.push_back(tile);
    }

    /**
     * Add a tile
     */
    void push_back(const ImageTile<T> &tile) {
      tiles_.push_back(tile);
    }

    /**
     * Get the tile names
     */
    scitbx::af::shared< std::string > tile_names() const {
      scitbx::af::shared<std::string> names;
      for (std::size_t i = 0; i < tiles_.size(); ++i) {
        names.push_back(tiles_[i].name());
      }
      return names;
    }

    /**
     * Get an image tile
     */
    ImageTile<T> tile(std::size_t index) const {
      DXTBX_ASSERT(index < n_tiles());
      return tiles_[index];
    }

    /**
     * Get the number of tiles
     */
    std::size_t n_tiles() const {
      return tiles_.size();
    }

    /**
     * Is the image empty
     */
    bool empty() const {
      return tiles_.empty();
    }

    /**
     * Get the begin iterator
     */
    iterator begin() {
      return tiles_.begin();
    }

    /**
     * Get the end iterator
     */
    iterator end() {
      return tiles_.end();
    }

  protected:

    scitbx::af::shared< ImageTile<T> > tiles_;

  };


  /**
   * A class to hold image data which can be either int or double
   */
  class ImageBuffer {
  public:

    typedef int empty_type;
    typedef Image<int> int_image_type;
    typedef Image<double> double_image_type;

    // The variant type
    typedef boost::variant <
      empty_type,
      int_image_type,
      double_image_type
    > variant_type;

    /**
     * A visitor class to convert from/to different types.
     * Data is copied except when the to/from types are the same
     */
    template <typename ImageType>
    class ConverterVisitor : public boost::static_visitor<ImageType> {
    public:

      ImageType operator()(const empty_type &) const {
        throw DXTBX_ERROR("ImageBuffer is empty");
        return ImageType();
      }

      ImageType operator()(const ImageType &v) const {
        return v;
      }

      template <typename OtherImageType>
      ImageType operator()(const OtherImageType &v) const {
        ImageType result;
        for (std::size_t i = 0; i < v.n_tiles(); ++i) {
          typedef typename ImageType::tile_type ImageTileType;
          typedef typename ImageType::array_type ArrayType;
          ArrayType data(
              v.tile(i).accessor(),
              scitbx::af::init_functor_null<typename ArrayType::value_type>());
          std::uninitialized_copy(v.tile(i).data().begin(),
                                  v.tile(i).data().end(), data.begin());
          result.push_back(ImageTileType(data));
        }
        return result;
      }

    };

    /**
     * Is the buffer empty
     */
    class IsEmptyVisitor : public boost::static_visitor<bool> {
    public:

      bool operator()(const empty_type &v) const {
        return true;
      }

      template <typename OtherImageType>
      bool operator()(const OtherImageType &v) const {
        return false;
      }

    };

    /**
     * Is the data an int type
     */
    class IsIntVisitor : public boost::static_visitor<bool> {
    public:

      bool operator()(const int_image_type &v) const {
        return true;
      }

      template <typename OtherImageType>
      bool operator()(const OtherImageType &v) const {
        return false;
      }

    };

    /**
     * Is the data a double type
     */
    class IsDoubleVisitor : public boost::static_visitor<bool> {
    public:

      bool operator()(const double_image_type &v) const {
        return true;
      }

      template <typename OtherImageType>
      bool operator()(const OtherImageType &v) const {
        return false;
      }

    };

    /**
     * Construct an empty buffer
     */
    ImageBuffer() {}

    /**
     * Construct a buffer with data
     */
    ImageBuffer(variant_type data)
      : data_(data) {}

    /**
     * @returns Is the buffer empty
     */
    bool is_empty() const {
      return boost::apply_visitor(IsEmptyVisitor(), data_);
    }

    /**
     * @returns Is the buffer an int
     */
    bool is_int() const {
      return boost::apply_visitor(IsIntVisitor(), data_);
    }

    /**
     * @returns Is the buffer a double
     */
    bool is_double() const {
      return boost::apply_visitor(IsDoubleVisitor(), data_);
    }

    /**
     * @returns The buffer as an int image
     */
    Image<int> as_int() const {
      return boost::apply_visitor(ConverterVisitor< Image<int> >(), data_);
    }

    /**
     * @returns The buffer as a double image
     */
    Image<double> as_double() const {
      return boost::apply_visitor(ConverterVisitor< Image<double> >(), data_);
    }

  protected:

    variant_type data_;

  };

}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_IMAGE_H
