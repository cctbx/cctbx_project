/*
 * smv.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_SMV_READER_H
#define DXTBX_FORMAT_SMV_READER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <dxtbx/format/image.h>
#include <dxtbx/error.h>

#ifdef _WIN32  // use the unistd.h file from cbflib_adaptbx since VC++9.0 hasn't got one
#if _MSC_VER < 1600
#include <unistd.h>
#endif
#endif

namespace dxtbx { namespace format {

  namespace detail {

    /**
     * Check is a string starts with a substring
     */
    inline
    bool startswith(const std::string s, const std::string prefix) {
      return s.substr(0, prefix.size()) == prefix;
    }

    /**
     * Get a header item
     */
    inline
    bool get_header_item(std::string item, std::string &name, std::string &value) {
      std::size_t pos1 = item.find("=");
      if (pos1 != std::string::npos) {
        name = item.substr(0, pos1);
        item = item.substr(pos1+1,item.size());
        std::size_t pos2 = item.find(";");
        if (pos2 != std::string::npos) {
          value = item.substr(0, pos2);
          return true;
        }
      }
      return false;
    }

    /**
     * Get a value from a string
     */
    template <typename T>
    T get_value(std::string item) {
      T result;
      std::istringstream stream(item);
      stream >> result;
      return result;
    }

    /**
     * Check if the system is big endian
     */
    inline
    bool is_big_endian()
    {
      uint32_t val = 1;
      char * buff = (char *) & val;
      return (buff[0] == 0);
    }

    /**
     * Check if the byte order string is big endian
     */
    inline
    bool is_big_endian(std::string byte_order) {
      if (byte_order == "big_endian") {
        return true;
      } else if (byte_order != "little_endian") {
        DXTBX_ERROR("Unknown byte order");
      }
      return false;
    }

    /**
     * Swap bytes for an array of shorts
     */
    inline
    void swap_bytes(short *first, short *last) {
      for (short *it = first; it != last; ++it) {
        short x = *it;
        *it = x >> 8 | x << 8;
      }
    }

    /**
     * Swap bytes for an array of unsigned shorts
     */
    inline
    void swap_bytes(unsigned short *first, unsigned short *last) {
      for (unsigned short *it = first; it != last; ++it) {
        unsigned short x = *it;
        *it = x >> 8 | x << 8;
      }
    }

    /**
     * Swap bytes for an array of longs
     */
    inline
    void swap_bytes(int *first, int *last) {
      for (int *it = first; it != last; ++it) {
        int x = *it;
        *it = (x << 24) |
              (x << 8 & 0xff0000) |
              (x >> 8 & 0xff00) |
              (x >> 24);
      }
    }

    /**
     * Swap bytes for an array of unsigned longs
     */
    inline
    void swap_bytes(unsigned int *first, unsigned int *last) {
      for (unsigned int *it = first; it != last; ++it) {
        unsigned int x = *it;
        *it = (x << 24) |
              (x << 8 & 0xff0000) |
              (x >> 8 & 0xff00) |
              (x >> 24);
      }
    }

  }


  /**
   * A class to read an SMV Image
   */
  class SMVReader : public ImageReader {
  public:

    /**
     * Construct the class with the filename
     */
    SMVReader(const char *filename)
      : ImageReader(filename),
        slow_size_(0),
        fast_size_(0) {
      std::ifstream handle(filename, std::ifstream::binary);
      read_header(handle);
      read_data(handle);
    }

    /**
     * Return the byte order
     */
    std::string byte_order() const {
      return byte_order_;
    }


  protected:

    // Helper structs to get types
    template <typename T>
    struct array_type {
      typedef scitbx::af::versa< T, scitbx::af::c_grid<2> > type;
    };

    /**
     * Read the number of bytes in the header
     */
    std::size_t read_header_nbytes(std::ifstream &handle) {

      // Read first couple of lines and check
      std::string line;
      std::getline(handle, line);
      DXTBX_ASSERT(detail::startswith(line, "{"));
      std::getline(handle, line);
      DXTBX_ASSERT(detail::startswith(line, "HEADER_BYTES="));

      // Get the number of header bytes
      std::string name;
      std::string value;
      detail::get_header_item(line, name, value);
      std::size_t nbytes = detail::get_value<std::size_t>(value);

      // Go to the beginning
      handle.seekg(0, std::ios_base::beg);
      DXTBX_ASSERT(handle.good());

      // Return the number of bytes
      return nbytes;
    }

    /**
     * Read the image header
     */
    void read_header(std::ifstream &handle) {

      // Open the file and check is open
      DXTBX_ASSERT(handle.is_open());

      // Get the number of bytes in the header
      std::size_t nbytes = read_header_nbytes(handle);

      // Read the bytes into a string stream
      std::vector<char> buffer(nbytes);
      handle.read(&buffer[0], nbytes);
      handle.seekg(nbytes, std::ios_base::beg);
      DXTBX_ASSERT(handle.good());

      // Create a string stream
      std::stringstream header(&buffer[0]);
      std::string line;
      while (std::getline(header, line)) {

        // End of the header
        if (detail::startswith(line, "}")) {
          break;
        }

        // Get the header item
        std::string name;
        std::string value;
        if (detail::get_header_item(line, name, value)) {
          if (name == "DIM") {
            DXTBX_ASSERT(detail::get_value<std::size_t>(value) == 2);
          } else if (name == "TYPE") {
            if (value == "signed_short") {
              type_ = "int16";
            } else if (value == "unsigned_short") {
              type_ = "uint16";
            } else if (value == "signed short int") {
              type_ = "int16";
            } else if (value == "unsigned short int") {
              type_ = "uint16";
            }
          } else if (name == "Data_type") {
            if (value == "signed_short") {
              type_ = "int16";
            } else if (value == "unsigned_short") {
              type_ = "uint16";
            } else if (value == "signed short int") {
              type_ = "int16";
            } else if (value == "unsigned short int") {
              type_ = "uint16";
            }
          } else if (name == "BYTE_ORDER") {
            DXTBX_ASSERT(
                value == "big_endian" ||
                value == "little_endian");
            byte_order_ = value;
                  } else if (name == "SIZE1") {
            slow_size_ = detail::get_value<std::size_t>(value);
          } else if (name == "SIZE2") {
            fast_size_ = detail::get_value<std::size_t>(value);
          } else if (name == "COMPRESSION") {
            DXTBX_ASSERT(value == "None");
          }
        }
      }

      // Check we have everything
      DXTBX_ASSERT(type_ != "");
      DXTBX_ASSERT(byte_order_ != "");
      DXTBX_ASSERT(slow_size_ > 0);
      DXTBX_ASSERT(fast_size_ > 0);
    }

    /**
     * Read the image data
     */
    void read_data(std::ifstream &handle) {

      // Get the element size
      if (type_ == "int16") {
        read_data_detail<int, short>(handle);
      } else if (type_ == "uint16") {
        read_data_detail<int, unsigned short>(handle);
      } else if (type_ == "int32") {
        read_data_detail<int, int>(handle);
      } else if (type_ == "uint32") {
        read_data_detail<int, unsigned int>(handle);
      } else {
        DXTBX_ASSERT("Unsupported type");
      }
    }

    template <typename OutputType, typename InputType>
    void read_data_detail(std::ifstream &handle) {

      typedef typename array_type<InputType>::type input_array_data_type;
      typedef typename array_type<OutputType>::type output_array_data_type;

      // The image grid
      DXTBX_ASSERT(slow_size_ > 0);
      DXTBX_ASSERT(fast_size_ > 0);
      scitbx::af::c_grid<2> grid(slow_size_, fast_size_);

      // Allocate the array
      std::size_t element_size = sizeof(InputType);
      std::size_t nbytes = element_size * slow_size_ * fast_size_;
      input_array_data_type data(grid);

      // Read the data
      handle.read(reinterpret_cast<char*>(&data[0]), nbytes);

      // Swap bytes if necessary
      if (detail::is_big_endian(byte_order_) != detail::is_big_endian()) {
        detail::swap_bytes(data.begin(), data.end());
      }

      // Check handle is still good
      DXTBX_ASSERT(handle.good());

      // Add to the tiles list
      output_array_data_type output(data.accessor());
      std::copy(data.begin(), data.end(), output.begin());
      buffer_ = ImageBuffer(Image<OutputType>(ImageTile<OutputType>(output, "")));
    }

    std::size_t slow_size_;
    std::size_t fast_size_;
    std::string type_;
    std::string byte_order_;
  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_SMV_READER_H
