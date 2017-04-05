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
#ifndef DXTBX_FORMAT_SMV_H
#define DXTBX_FORMAT_SMV_H

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
    void swap_bytes(long *first, long *last) {
      for (long *it = first; it != last; ++it) {
        long x = *it;
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
    void swap_bytes(unsigned long *first, unsigned long *last) {
      for (unsigned long *it = first; it != last; ++it) {
        unsigned long x = *it;
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

    typedef ImageReader::size_type size_type;

    /**
     * Construct the class with the filename
     */
    SMVReader(const char *filename)
      : ImageReader(filename) {
      std::ifstream handle(filename, std::ifstream::binary);
      read_header(handle);
      read_data(handle);
    }

    /**
     * Get the image data as the desired type
     */
    template<typename T>
    scitbx::af::versa< T, scitbx::af::c_grid<2> > as_type() const {
      scitbx::af::c_grid<2> grid(size_[0], size_[1]);
      scitbx::af::versa< T, scitbx::af::c_grid<2> > result(grid);

      // Get the number of elements
      std::size_t n_elements = size_[0] * size_[1];

      // Copy the data into the output array
      if (type_ == "int16") {
        const short *buffer = reinterpret_cast<const short*>(&data_[0]);
        DXTBX_ASSERT(data_.size() == n_elements * sizeof(short));
        std::copy(buffer, buffer + n_elements, result.begin());
      } else if (type_ == "uint16") {
        const unsigned short *buffer = reinterpret_cast<const unsigned short*>(&data_[0]);
        DXTBX_ASSERT(data_.size() == n_elements * sizeof(unsigned short));
        std::copy(buffer, buffer + n_elements, result.begin());
      } else if (type_ == "int32") {
        const long *buffer = reinterpret_cast<const long*>(&data_[0]);
        DXTBX_ASSERT(data_.size() == n_elements * sizeof(long));
        std::copy(buffer, buffer + n_elements, result.begin());
      } else if (type_ == "uint32") {
        const unsigned long *buffer = reinterpret_cast<const unsigned long*>(&data_[0]);
        DXTBX_ASSERT(data_.size() == n_elements * sizeof(unsigned long));
        std::copy(buffer, buffer + n_elements, result.begin());
      }

      // Return the result
      return result;
    }

    /**
     * Get the image data as an integer array
     */
    scitbx::af::versa< int, scitbx::af::c_grid<2> > as_int() const {
      return as_type<int>();
    }

    /**
     * Get the image data as a double array
     */
    scitbx::af::versa< double, scitbx::af::c_grid<2> > as_double() const {
      return as_type<double>();
    }

  protected:

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
            size_[0] = detail::get_value<std::size_t>(value);
          } else if (name == "SIZE2") {
            size_[1] = detail::get_value<std::size_t>(value);
          } else if (name == "COMPRESSION") {
            DXTBX_ASSERT(value == "None");
          }
        }
      }

      // Check we have everything
      DXTBX_ASSERT(type_ != "");
      DXTBX_ASSERT(byte_order_ != "");
      DXTBX_ASSERT(size_.all_gt(0));
    }

    /**
     * Read the image data
     */
    void read_data(std::ifstream &handle) {

      // Get the element size
      std::size_t element_size = 0;
      if (type_ == "int16") {
        element_size = sizeof(short);
      } else if (type_ == "uint16") {
        element_size = sizeof(unsigned short);
      } else if (type_ == "int32") {
        element_size = sizeof(long);
      } else if (type_ == "uint32") {
        element_size = sizeof(unsigned long);
      } else {
        DXTBX_ASSERT("Unsupported type");
      }

      // The number of bytes
      std::size_t n_elements = size_[0] * size_[1];
      std::size_t nbytes = n_elements * element_size;
      DXTBX_ASSERT(nbytes > 0);

      // Allocate the data buffer
      data_.resize(nbytes);

      // Read the data
      handle.read(&data_[0], nbytes);
      DXTBX_ASSERT(handle.good());

      // If we need to swap bytes
      if (detail::is_big_endian(byte_order_) != detail::is_big_endian()) {
        if (type_ == "int16") {
          short *buffer = reinterpret_cast<short*>(&data_[0]);
          detail::swap_bytes(buffer, buffer + n_elements);
        } else if (type_ == "uint16") {
          unsigned short *buffer = reinterpret_cast<unsigned short*>(&data_[0]);
          detail::swap_bytes(buffer, buffer + n_elements);
        } else if (type_ == "int32") {
          long *buffer = reinterpret_cast<long*>(&data_[0]);
          detail::swap_bytes(buffer, buffer + n_elements);
        } else if (type_ == "uint32") {
          unsigned long *buffer = reinterpret_cast<unsigned long*>(&data_[0]);
          detail::swap_bytes(buffer, buffer + n_elements);
        }
      }
    }

    std::vector<char> data_;
  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_SMV_H
