/*
 * cbf_reader.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_CBF_READER_H
#define DXTBX_FORMAT_CBF_READER_H

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <include/cbf.h>
#include <dxtbx/format/image.h>
#include <dxtbx/error.h>
#include <cbflib_adaptbx/detectors/buffer_based_service.h>

#define cbf_check(x) DXTBX_ASSERT((x) == 0)

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


    template <typename T>
    struct cbf_array_buffer {};

    /**
     * Helper struct to read an integer cbf buffer
     */
    template<>
    struct cbf_array_buffer<int> {

      std::vector<int> data;
      std::size_t dimfast;
      std::size_t dimmid;
      std::size_t dimslow;

      cbf_array_buffer(cbf_handle cbf_h)
        : dimfast(0),
          dimmid(0),
          dimslow(0) {

        // Get the parameters
        unsigned int compression = 0;
        int binary_id = 0;
        size_t elsize = 0;
        int elsigned = 0;
        int elunsigned = 0;
        size_t elements = 0;
        int minelement = 0;
        int maxelement = 0;
        const char *byteorder = 0;
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
        DXTBX_ASSERT(elsigned == 1);
        if (dimslow == 0) dimslow = 1;
        DXTBX_ASSERT(std::string(byteorder) == "little_endian");
        DXTBX_ASSERT(dimfast * dimmid * dimslow == elements);
        DXTBX_ASSERT(elements > 0);

        // Resize
        data.resize(elements);

        // Read the data
        size_t elements_read;
        cbf_get_integerarray(
          cbf_h,
          &binary_id,
          &data[0],
          sizeof(int),
          elsigned,
          elements,
          &elements_read);
        DXTBX_ASSERT(elements_read == elements);
      }

    };

    /**
     * Helper struct to read a double cbf buffer
     */
    template<>
    struct cbf_array_buffer<double> {

      std::vector<double> data;
      std::size_t dimfast;
      std::size_t dimmid;
      std::size_t dimslow;

      cbf_array_buffer(cbf_handle cbf_h)
        : dimfast(0),
          dimmid(0),
          dimslow(0) {

        // Get parameters
        unsigned int compression = 0;
        int binary_id = 0;
        size_t elsize = 0;
        size_t elements = 0;
        const char *byteorder = 0;
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
        if (dimslow == 0) dimslow = 1;
        DXTBX_ASSERT(elsize = sizeof(double));
        DXTBX_ASSERT(std::string(byteorder) == "little_endian");
        DXTBX_ASSERT(dimfast * dimmid * dimslow == elements);
        DXTBX_ASSERT(elements > 0);

        // Resize the data
        data.resize(elements);

        // Read the data
        size_t elements_read;
        cbf_get_realarray(
          cbf_h,
          &binary_id,
          &data[0],
          sizeof(double),
          elements,
          &elements_read);
        DXTBX_ASSERT(elements_read == elements);
      }

    };

  }


  /**
   * A class to read a CBF Image (quickly)
   * Only works for signed 16/32 bit int data
   */
  class CBFFastReader : public ImageReader {
  public:

    /**
     * Construct the class with the filename
     */
    CBFFastReader(const char *filename)
      : ImageReader(filename) {
      read_data();
    }

  protected:

    /**
     * Read the image data
     */
    void read_data() {

      // Open the file handle
      std::ifstream handle(filename_.c_str(), std::ifstream::binary);

      // Read all the data into a buffer
      std::string buffer;
      std::ifstream::pos_type file_start = handle.tellg();
      handle.seekg(0, std::ios::end);
      std::ifstream::pos_type file_end = handle.tellg();
      DXTBX_ASSERT(file_end > file_start);
      std::size_t file_size = file_end - file_start;
      DXTBX_ASSERT(file_size > 0);
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
      std::string type = "int32";
      std::size_t length = 0;
      std::size_t data_size = 0;
      std::size_t fast_size = 0;
      std::size_t slow_size = 0;
      bool byte_offset = false;

      // Parse the header
      std::string line;
      while (std::getline(header, line)) {
        std::string name;
        std::string value;
        if (detail::get_cbf_header_item(line, name, value)) {
          if (name == "X-Binary-Size-Fastest-Dimension") {
            fast_size = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Size-Second-Dimension") {
            slow_size = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Number-of-Elements") {
            length = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Size") {
            data_size = detail::get_cbf_header_value<std::size_t>(value);
          } else if (name == "X-Binary-Element-Type") {
            if (value == "signed 16-bit integer") {
              type = "int16";
            } else if (value == "signed 32-bit integer") {
              type = "int32";
            } else {
              DXTBX_ERROR("Can only handle signed 16/32-bit integer data");
            }
          } else if (name == "X-Binary-Element-Byte-Order") {
            DXTBX_ASSERT(value == "LITTLE_ENDIAN");
          }
        } else if (
          line.find("conversions") != std::string::npos &&
          line.find("x-CBF_BYTE_OFFSET")) {
            byte_offset = true;
        }
      }

      // Check the input
      DXTBX_ASSERT(fast_size > 0);
      DXTBX_ASSERT(slow_size > 0);
      DXTBX_ASSERT(slow_size * fast_size == length);
      DXTBX_ASSERT(byte_offset == true);
      DXTBX_ASSERT(data_size >= length);

      // Resize the array
      scitbx::af::c_grid<2> grid(slow_size, fast_size);
      scitbx::af::versa< int, scitbx::af::c_grid<2> > data(grid);

      // Uncompress the data
      iotbx::detectors::buffer_uncompress(
          buffer.c_str() + data_offset,
          data_size,
          &data[0]);

      // Add to the tiles
      buffer_ = ImageBuffer(Image<int>(ImageTile<int>(data)));
    }

  };


  /**
   * A class to read a CBF Image
   */
  class CBFReader : public ImageReader {
  public:

    /**
     * Construct the class with the filename
     */
    CBFReader(const char *filename)
      : ImageReader(filename) {
      read_data();
    }

  protected:

    /**
     * Helper struct for multi-tile CBF files
     */
    struct Section {
      std::string name;
      int slow_start;
      int slow_end;
      int mid_start;
      int mid_end;
      int fast_start;
      int fast_end;

      Section():
        slow_start(-1),
        slow_end(-1),
        mid_start(-1),
        mid_end(-1),
        fast_start(-1),
        fast_end(-1) {}
    };

    /**
     * Read the data
     */
    void read_data() {

      // Open the cbf handle
      cbf_handle cbf_h;
      cbf_check(cbf_make_handle(&cbf_h));
      DXTBX_ASSERT(cbf_h);

      // Open the file
      FILE *handle = std::fopen(filename_.c_str(), "rb");
      DXTBX_ASSERT(handle);

      // Read the file
      cbf_check(cbf_read_widefile(cbf_h, handle, MSG_DIGEST));

      // Find the array structure
      if (cbf_find_category(cbf_h, "array_structure") == 0) {
        cbf_check(cbf_find_column(cbf_h, "encoding_type"));
        cbf_check(cbf_select_row(cbf_h, 0));

        // Get the array data type for each section
        std::string array_data_type;
        unsigned int nrows = 0;
        cbf_check(cbf_count_rows(cbf_h, &nrows));
        for (std::size_t i = 0; i < nrows; ++i) {
          const char *type = 0;
          cbf_check(cbf_get_value(cbf_h, &type));
          if (i == 0) {
            array_data_type = type;
            DXTBX_ASSERT(
                array_data_type == "signed 32-bit integer" ||
                array_data_type == "signed 64-bit real IEEE");
          } else {
            DXTBX_ASSERT(array_data_type == type);
          }
          cbf_check(cbf_next_row(cbf_h));
        }

        // Call appropriate function depending on the type
        if (array_data_type == "signed 32-bit integer") {
          read_multi_tile_data_detail<int>(cbf_h);
        } else if (array_data_type == "signed 64-bit real IEEE") {
          read_multi_tile_data_detail<double>(cbf_h);
        } else {
          DXTBX_ERROR("Unsupported data type");
        }
      } else {
        read_single_tile_data_detail(cbf_h);
      }

      // Free the CBF handle
      cbf_check(cbf_free_handle(cbf_h));
    }

    /**
     * Read single tile data
     */
    void read_single_tile_data_detail(cbf_handle cbf_h) {

      // Get the data
      cbf_check(cbf_find_category(cbf_h, "array_data"));
      cbf_check(cbf_find_column(cbf_h, "data"));

      const char *data_type = 0;
      cbf_check(cbf_get_typeofvalue(cbf_h, &data_type));
      DXTBX_ASSERT(data_type != 0);
      DXTBX_ASSERT(std::string(data_type) == "bnry");

      // Read the data
      detail::cbf_array_buffer<int> buffer(cbf_h);
      DXTBX_ASSERT(buffer.dimslow == 1);

      // Allocate and copy the data
      scitbx::af::c_grid<2> grid(buffer.dimmid, buffer.dimfast);
      scitbx::af::versa< int, scitbx::af::c_grid<2> > data(grid);
      DXTBX_ASSERT(buffer.data.size() == data.size());
      std::copy(buffer.data.begin(), buffer.data.end(), data.begin());

      // Add to the tiles
      buffer_ = ImageBuffer(Image<int>(ImageTile<int>(data)));
    }

    /**
     * Read multi-tile data
     */
    template <typename T>
    void read_multi_tile_data_detail(cbf_handle &cbf_h) {

      std::map < std::string, std::vector<Section> > section_lookup;
      bool has_sections = false;

      // Check if the data has sections
      if (cbf_find_category(cbf_h, "array_structure_list_section") == 0) {

        has_sections = true;

        // Count the number of sections
        unsigned int nsections = 0;
        cbf_check(cbf_count_rows(cbf_h, &nsections));
        for (std::size_t i = 0; i < nsections; ++i) {

          // Find the section ID and get the name
          cbf_check(cbf_find_column(cbf_h, "id"));
          const char *section_name = 0;
          cbf_check(cbf_get_value(cbf_h, &section_name));
          DXTBX_ASSERT(section_name != 0);

          // Find the array id
          cbf_check(cbf_find_column(cbf_h, "array_id"));
          const char *array_id = 0;
          cbf_check(cbf_get_value(cbf_h, &array_id));
          DXTBX_ASSERT(array_id != 0);

          // Get the section and append to the list
          int axis_index;
          int axis_start;
          int axis_end;

          cbf_check(cbf_find_column(cbf_h, "index"));
          cbf_check(cbf_get_integervalue(cbf_h, &axis_index));

          cbf_check(cbf_find_column(cbf_h, "start"));
          cbf_check(cbf_get_integervalue(cbf_h, &axis_start));

          cbf_check(cbf_find_column(cbf_h, "end"));
          cbf_check(cbf_get_integervalue(cbf_h, &axis_end));

          // Add the section
          std::vector<Section> &sections = section_lookup[array_id];
          int index = -1;
          for (std::size_t j = 0; j < sections.size(); ++j) {
            if (sections[j].name == section_name) {
              index = j;
              break;
            }
          }
          if (index >= 0) {
            Section &s = sections[index];
            if (axis_index == 3) {
              s.slow_start = axis_start - 1;
              s.slow_end = axis_end;
            } else if (axis_index == 2) {
              s.mid_start = axis_start - 1;
              s.mid_end = axis_end;
            } else if (axis_index == 1) {
              s.fast_start = axis_start - 1;
              s.fast_end = axis_end;
            } else {
              DXTBX_ASSERT("Bad axis index");
            }
          } else {
            Section s;
            s.name = section_name;
            if (axis_index == 3) {
              s.slow_start = axis_start - 1;
              s.slow_end = axis_end;
            } else if (axis_index == 2) {
              s.mid_start = axis_start - 1;
              s.mid_end = axis_end;
            } else if (axis_index == 1) {
              s.fast_start = axis_start - 1;
              s.fast_end = axis_end;
            } else {
              DXTBX_ASSERT("Bad axis index");
            }
            sections.push_back(s);
          }

          // Get the next section
          cbf_check(cbf_next_row(cbf_h));
        }
      }

      // Find the array data
      unsigned int nrows = 0;
      cbf_check(cbf_find_category(cbf_h, "array_data"));
      cbf_check(cbf_count_rows(cbf_h, &nrows));
      DXTBX_ASSERT(nrows > 0);
      Image<T> image;
      if (nrows == 1) {

        // Get the data
        cbf_check(cbf_find_column(cbf_h, "data"));

        const char *data_type = 0;
        cbf_check(cbf_get_typeofvalue(cbf_h, &data_type));
        DXTBX_ASSERT(data_type != 0);
        DXTBX_ASSERT(std::string(data_type) == "bnry");

        // Read the data
        detail::cbf_array_buffer<T> buffer(cbf_h);

        DXTBX_ASSERT(buffer.dimslow == 1);

        // Allocate and copy the data
        scitbx::af::c_grid<2> grid(buffer.dimmid, buffer.dimfast);
        scitbx::af::versa< T, scitbx::af::c_grid<2> > data(grid);
        DXTBX_ASSERT(buffer.data.size() == data.size());
        std::copy(buffer.data.begin(), buffer.data.end(), data.begin());

        // Add to the tiles
        image.push_back(ImageTile<T>(data));

      } else {
        for (std::size_t i = 0; i < nrows; ++i) {

          // Get the array name
          cbf_check(cbf_find_column(cbf_h, "array_id"));
          const char *name = 0;
          cbf_check(cbf_get_value(cbf_h, &name));
          DXTBX_ASSERT(name != 0);

          // Get the data
          cbf_check(cbf_find_column(cbf_h, "data"));

          const char *data_type = 0;
          cbf_check(cbf_get_typeofvalue(cbf_h, &data_type));
          DXTBX_ASSERT(data_type != 0);
          DXTBX_ASSERT(std::string(data_type) == "bnry");

          // Read the data
          detail::cbf_array_buffer<T> buffer(cbf_h);

          if (has_sections) {

            // Get the sections
            std::vector<Section> sections = section_lookup[name];

            // Loop through the sections
            for (std::size_t j = 0; j < sections.size(); ++j) {

              // Get the section info
              std::string section_name = sections[j].name;
              int slow_start = sections[j].slow_start;
              int slow_end = sections[j].slow_end;
              int mid_start = sections[j].mid_start;
              int mid_end = sections[j].mid_end;
              int fast_start = sections[j].fast_start;
              int fast_end = sections[j].fast_end;

              // Get the size of the sections
              int slow_size = slow_end - slow_start;
              int mid_size = mid_end - mid_start;
              int fast_size = fast_end - fast_start;

              // Check the size of the sections
              DXTBX_ASSERT(slow_size == 1);
              DXTBX_ASSERT(mid_size > 0);
              DXTBX_ASSERT(fast_size > 0);
              DXTBX_ASSERT(fast_start >= 0);
              DXTBX_ASSERT(slow_start >= 0);
              DXTBX_ASSERT(mid_start >= 0);
              DXTBX_ASSERT(fast_end <= buffer.dimfast);
              DXTBX_ASSERT(mid_end <= buffer.dimmid);
              DXTBX_ASSERT(slow_end <= buffer.dimslow);

              // Allocate
              scitbx::af::c_grid<2> grid(mid_size, fast_size);
              scitbx::af::versa< T, scitbx::af::c_grid<2> > data(grid);

              // Copy the data from the buffer
              std::size_t size1 = buffer.dimfast;
              std::size_t size2 = buffer.dimmid * buffer.dimfast;
              for (std::size_t y = 0; y < mid_size; ++y) {
                for (std::size_t x = 0; x < fast_size; ++x) {
                  std::size_t xx = x + fast_start;
                  std::size_t yy = y + mid_start;
                  std::size_t zz = slow_start;
                  std::size_t k = xx + yy*size1 + zz*size2;
                  DXTBX_ASSERT(k < buffer.data.size());
                  data(y,x) = buffer.data[k];
                }
              }

              // Add to the tiles
              image.push_back(ImageTile<T>(data, section_name.c_str()));
            }

          } else {
            DXTBX_ASSERT(buffer.dimslow == 1);

            // Allocate and copy the data
            scitbx::af::c_grid<2> grid(buffer.dimmid, buffer.dimfast);
            scitbx::af::versa< T, scitbx::af::c_grid<2> > data(grid);
            DXTBX_ASSERT(buffer.data.size() == data.size());
            std::copy(buffer.data.begin(), buffer.data.end(), data.begin());

            // Add to the tiles
            image.push_back(ImageTile<T>(data, name == 0 ? "" : name));
          }

          // Go to the next array
          cbf_check(cbf_next_row(cbf_h));
        }
      }

      buffer_ = ImageBuffer(image);
    }

  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_CBF_READER_H
