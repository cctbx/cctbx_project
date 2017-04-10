/*
 * cbf.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_CBF_H
#define DXTBX_FORMAT_CBF_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <include/cbf.h>
#include <dxtbx/format/image.h>
#include <dxtbx/error.h>
#include <cbflib_adaptbx/detectors/buffer_based_service.h>

#define cbf_failnez(x) { int err; err = (x); if (err) { \
  std::cout<<"error code "<<err<<std::endl; DXTBX_ERROR("CBFlib error in " #x " "); }}

namespace dxtbx { namespace format {

  namespace detail {

    /**
     * trim from end of string (right)
     */
    inline
    std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v") {
      s.erase(s.find_last_not_of(t) + 1);
      return s;
    }

    /**
     * trim from beginning of string (left)
     */
    inline
    std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v") {
      s.erase(0, s.find_first_not_of(t));
      return s;
    }

    /**
     * trim from both ends of string (left & right)
     */
    inline
    std::string& trim(std::string& s, const char* t = " \t\n\r\f\v") {
      return ltrim(rtrim(s, t), t);
    }

    /**
     * Get a header item
     */
    inline
    bool get_cbf_header_item(std::string item, std::string &name, std::string &value) {
      std::size_t pos1 = item.find(":");
      if (pos1 != std::string::npos) {
        name = item.substr(0, pos1);
        value = item.substr(pos1+1, item.size());
        trim(name);
        trim(value);
        return true;
      }
      return false;
    }

    /**
     * Get a value from a string
     */
    template <typename T>
    T get_cbf_header_value(std::string item) {
      T result;
      std::istringstream stream(item);
      stream >> result;
      return result;
    }
  }


  /**
   * A class to read a CBF Image (quickly)
   */
  class CBFFastReader : public ImageReader {
  public:

    typedef ImageReader::size_type size_type;

    /**
     * Construct the class with the filename
     */
    CBFFastReader(const char *filename)
      : ImageReader(filename) {
      read_data();
    }

    /**
     * Get the image data as the desired type
     */
    template<typename T>
    scitbx::af::versa< T, scitbx::af::c_grid<2> > as_type() const {
      scitbx::af::c_grid<2> grid(size_[0], size_[1]);
      scitbx::af::versa< T, scitbx::af::c_grid<2> > result(grid);
      DXTBX_ASSERT(data_.accessor().all_eq(result.accessor()));
      std::copy(data_.begin(), data_.end(), result.begin());
      return result;
    }

    /**
     * Get the image data as an integer array
     */
    scitbx::af::versa< int, scitbx::af::c_grid<2> > as_int() const {
      return data_;
    }

    /**
     * Get the image data as a double array
     */
    scitbx::af::versa< double, scitbx::af::c_grid<2> > as_double() const {
      return as_type<double>();
    }

  protected:

    /**
     * Read the image data
     */
    void read_data() {

      // Open the file handle
      std::ifstream handle(filename_, std::ifstream::binary);

      // Read all the data into a buffer
      std::string buffer;
      handle.seekg(0, std::ios::end);
      std::ifstream::pos_type file_size = handle.tellg();
      buffer.resize(file_size);
      handle.seekg(0, std::ios::beg);
      handle.read(&buffer[0], file_size);

      // Find the data offset
      std::string start_tag = "\x0c\x1a\x04\xd5";
      std::size_t data_offset = buffer.find(start_tag) + 4;
      DXTBX_ASSERT(data_offset != std::string::npos);

      // Get a stream to the header information
      std::stringstream header(buffer.substr(0, data_offset));

      // Initialise some info
      type_ = "int32";
      std::size_t length = 0;
      std::size_t data_size = 0;
      bool byte_offset = false;

      // Parse the header
      std::string line;
      while (std::getline(header, line)) {
        std::string name;
        std::string value;
        if (detail::get_cbf_header_item(line, name, value)) {
          if (name == "X-Binary-Size-Fastest-Dimension") {
            size_[1] = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Size-Second-Dimension") {
            size_[0] = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Number-of-Elements") {
            length = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Size") {
            data_size = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Element-Type") {
            if (value == "signed 16-bit integer") {
              type_ = "int16";
            } else if (value == "signed 32-bit integer") {
              type_ = "int32";
            } else {
              DXTBX_ERROR("Can only handle signed 16/32-bit integer data");
            }
          } else if (name == "X-Binary-Element-Byte-Order") {
            DXTBX_ASSERT(value == "LITTLE_ENDIAN");
          }
        } else if (
          line.find("conversions") != std::string::npos &&
          line.find("x-CBF_BYTE_OFFSET")) {
            byte_offset=true;
        }
      }

      // Check the input
      DXTBX_ASSERT(size_[0] > 0);
      DXTBX_ASSERT(size_[1] > 0);
      DXTBX_ASSERT(size_[1] * size_[0] == length);
      DXTBX_ASSERT(byte_offset == true);
      DXTBX_ASSERT(data_size >= length);

      // Resize the array
      scitbx::af::c_grid<2> grid(size_);
      data_.resize(grid);

      // Uncompress the data
      iotbx::detectors::buffer_uncompress(
          buffer.c_str() + data_offset,
          data_size,
          &data_[0]);
    }

    scitbx::af::versa< int, scitbx::af::c_grid<2> > data_;
  };


  /**
   * A class to read a CBF Image
   */
  class CBFReader : public ImageReader {
  public:

    typedef ImageReader::size_type size_type;

    /**
     * Construct the class with the filename
     */
    CBFReader(const char *filename)
      : ImageReader(filename) {
      read_data();
    }

  protected:

    void read_data() {

      // Open the cbf handle
      cbf_handle cbf_h;
      cbf_failnez(cbf_make_handle(&cbf_h));
      DXTBX_ASSERT(cbf_h);

      FILE *handle = std::fopen(filename_.c_str(), "rb");
      DXTBX_ASSERT(handle);

      // Read the file
      cbf_failnez(cbf_read_widefile(cbf_h, handle, MSG_DIGEST));

      cbf_find_category(cbf_h, "array_structure");
      cbf_find_column(cbf_h, "encoding_type");
      cbf_select_row(cbf_h, 0);

      std::string array_data_type;
      unsigned int nrows = 0;
      cbf_count_rows(cbf_h, &nrows);
      for (std::size_t i = 0; i < nrows; ++i) {
        const char *type = 0;
        cbf_get_value(cbf_h, &type);
        if (i == 0) {
          array_data_type = type;
          DXTBX_ASSERT(
              array_data_type == "signed 32-bit integer" ||
              array_data_type == "signed 64-bit real IEEE");
        } else {
          DXTBX_ASSERT(array_data_type == type);
        }
        cbf_next_row(cbf_h);
      }

      cbf_find_category(cbf_h, "array_data");
      cbf_count_rows(cbf_h, &nrows);
      for (std::size_t i = 0; i < nrows; ++i) {
        cbf_find_column(cbf_h, "array_id");
        const char *name = 0;
        cbf_get_value(cbf_h, &name);

        cbf_find_column(cbf_h, "data");

        const char *data_type = 0;
        cbf_get_typeofvalue(cbf_h, &data_type);
        DXTBX_ASSERT(std::string(data_type) == "bnry");

        if (array_data_type == "signed 32-bit integer") {
          unsigned int compression = 0;
          int binary_id = 0;
          size_t elsize = 0;
          int elsigned = 0;
          int elunsigned = 0;
          size_t elements = 0;
          int minelement = 0;
          int maxelement = 0;
          const char *byteorder = 0;
          size_t dimfast = 0;
          size_t dimmid = 0;
          size_t dimslow = 0;
          size_t padding = 0;
          cbf_get_integerarrayparameters_wdims_fs(
              cbf_h,
              &compression,
              &binary_id,
              &elsize,
              &elsigned,
              &elunsigned,
              &elements,
              &minelement,
              &maxelement,
              &byteorder,
              &dimfast,
              &dimmid,
              &dimslow,
              &padding);
          DXTBX_ASSERT(std::string(byteorder) == "little_endian");
          DXTBX_ASSERT(dimfast * dimmid * dimslow == elements);
          DXTBX_ASSERT(elements > 0);

          std::vector<int> buffer(elements);

          size_t elements_read;
          cbf_get_integerarray(
              cbf_h,
              &binary_id,
              &buffer[0],
              elsize,
              elsigned,
              elements,
              &elements_read);
          DXTBX_ASSERT(elements_read == elements);
        } else if (array_data_type == "signed 64-bit real IEEE") {
            unsigned int compression = 0;
            int binary_id = 0;
            size_t elsize = 0;
            size_t elements = 0;
            const char *byteorder = 0;
            size_t dimfast = 0;
            size_t dimmid = 0;
            size_t dimslow = 0;
            size_t padding = 0;
          cbf_get_realarrayparameters_wdims_fs(
            cbf_h,
            &compression,
            &binary_id,
            &elsize,
            &elements,
            &byteorder,
            &dimfast,
            &dimmid,
            &dimslow,
            &padding);
          DXTBX_ASSERT(std::string(byteorder) == "little_endian");
          DXTBX_ASSERT(dimfast * dimmid * dimslow == elements);
          DXTBX_ASSERT(elements > 0);

          std::vector<double> buffer(elements);

          size_t elements_read;
          cbf_get_realarray(
              cbf_h,
              &binary_id,
              &buffer[0],
              elsize,
              elements,
              &elements_read);

          DXTBX_ASSERT(elements_read == elements);
        }

        cbf_next_row(cbf_h);
      }

      cbf_find_category(cbf_h, "array_structure_list_section");
      if (true) {

        cbf_count_rows(cbf_h, &nrows);
        for (std::size_t i = 0; i < nrows; ++i) {
          cbf_find_column(cbf_h, "id");
          const char *section_name = 0;
          cbf_get_value(cbf_h, &section_name);

          cbf_find_column(cbf_h, "array_id");

          int axis_index;
          int axis_start;
          int axis_end;

          cbf_find_column(cbf_h, "index");
          cbf_get_integervalue(cbf_h, &axis_index);

          cbf_find_column(cbf_h, "start");
          cbf_get_integervalue(cbf_h, &axis_start);

          cbf_find_column(cbf_h, "end");
          cbf_get_integervalue(cbf_h, &axis_end);

          cbf_next_row(cbf_h);
        }

      }

      // Free the CBF handle
      cbf_failnez(cbf_free_handle(cbf_h));
    }

  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_CBF_H
