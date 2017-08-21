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

    typedef scitbx::af::c_grid<2> accessor_type;

    // Variant data typedefs
    typedef scitbx::af::versa< bool,           accessor_type > bool_type;
    typedef scitbx::af::versa< short,          accessor_type > int16_type;
    typedef scitbx::af::versa< int,            accessor_type > int32_type;
    typedef scitbx::af::versa< unsigned short, accessor_type > uint16_type;
    typedef scitbx::af::versa< unsigned int,   accessor_type > uint32_type;
    typedef scitbx::af::versa< float,          accessor_type > float32_type;
    typedef scitbx::af::versa< double,         accessor_type > float64_type;

    typedef bool_type bool_array_type;
    typedef int32_type int_array_type;
    typedef float64_type double_array_type;

    typedef boost::variant <
      bool_type,
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

      ArrayType operator()(const ArrayType &v) const {
        return v;
      }

      template <typename OtherArrayType>
      ArrayType operator()(const OtherArrayType &v) const {
        ArrayType result(v.accessor());
        std::copy(v.begin(), v.end(), result.begin());
        return result;
      }

    };
    
    /**
     * Is the data a bool type
     */
    class IsBoolVisitor : public boost::static_visitor<bool> {
    public:

      bool operator()(const bool_type &v) const {
        return true;
      }

      template <typename OtherArrayType>
      bool operator()(const OtherArrayType &v) const {
        return false;
      }

    };

    /**
     * Is the data an int type
     */
    class IsIntVisitor : public boost::static_visitor<bool> {
    public:
      
      bool operator()(const int16_type &v) const {
        return true;
      }

      bool operator()(const int32_type &v) const {
        return true;
      }

      bool operator()(const uint16_type &v) const {
        return true;
      }

      bool operator()(const uint32_type &v) const {
        return true;
      }

      template <typename OtherArrayType>
      bool operator()(const OtherArrayType &v) const {
        return false;
      }

    };

    /**
     * Is the data a double type
     */
    class IsDoubleVisitor : public boost::static_visitor<bool> {
    public:

      bool operator()(const float32_type &v) const {
        return true;
      }

      bool operator()(const float64_type &v) const {
        return true;
      }

      template <typename OtherArrayType>
      bool operator()(const OtherArrayType &v) const {
        return false;
      }

    };
    
    /**
     * Is the array empty
     */
    class IsEmptyVisitor : public boost::static_visitor<bool> {
    public:

      template <typename OtherArrayType>
      bool operator()(const OtherArrayType &v) const {
        return v.empty();
      }

    };
    
    /**
     * Get the accessor
     */
    class AccessorVisitor : public boost::static_visitor<accessor_type> {
    public:

      template <typename OtherArrayType>
      accessor_type operator()(const OtherArrayType &v) const {
        return v.accessor();
      }

    };

    /**
     * Initialize the class
     */
    ImageTile(variant_type data)
      : data_(data),
        name_("") {}

    /**
     * Initialize the class
     */
    ImageTile(variant_type data, const char *name)
      : data_(data),
        name_(name) {}
    
    /**
     * Is the data boolean
     */
    bool is_bool() const {
      return boost::apply_visitor(IsBoolVisitor(), data_);
    }

    /**
     * Is the data integer
     */
    bool is_int() const {
      return boost::apply_visitor(IsIntVisitor(), data_);
    }

    /**
     * Is the data double
     */
    bool is_double() const {
      return boost::apply_visitor(IsDoubleVisitor(), data_);
    }

    /**
     * Get the data as bool
     */
    bool_array_type as_bool() const {
      DXTBX_ASSERT(is_bool());
      return boost::apply_visitor(ConverterVisitor<bool_array_type>(), data_);
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
    std::string name() const {
      return name_;
    }

    /**
     * Is the array empty
     */
    bool empty() const {
      return boost::apply_visitor(IsEmptyVisitor(), data_);
    }

    /**
     * Get the accessor
     */
    accessor_type accessor() const {
      return boost::apply_visitor(AccessorVisitor(), data_);
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
    typedef ImageTile::bool_type bool_type;
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
    void push_back(const ImageTile &tile) {
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
    ImageTile tile(std::size_t index) {
      DXTBX_ASSERT(index < n_tiles());
      return tiles_[index];
    }

    /**
     * Get the number of tiles
     */
    std::size_t n_tiles() const {
      return tiles_.size();
    }

  protected:

    scitbx::af::shared<ImageTile> tiles_;
  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_IMAGE_H
