/*
 * tiff.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_TIFF_H
#define DXTBX_FORMAT_TIFF_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <tiffio.h>
#include <dxtbx/format/image.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace format {

  namespace detail {

    /**
     * Ignore TIFF warnings (will print stuff about unknown tags otherwise)
     */
    void tiff_warning_handler(const char* module, const char* fmt, va_list ap) {}
  }


  /**
   * A class to read an TIFF Image
   */
  class TIFFReader : public ImageReader {
  public:

    typedef ImageReader::size_type size_type;

    /**
     * Construct the class with the filename
     */
    TIFFReader(const char *filename)
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
      } else if (type_ == "float32") {
        const float *buffer = reinterpret_cast<const float*>(&data_[0]);
        DXTBX_ASSERT(data_.size() == n_elements * sizeof(float));
        std::copy(buffer, buffer + n_elements, result.begin());
      } else if (type_ == "float64") {
        const double *buffer = reinterpret_cast<const double*>(&data_[0]);
        DXTBX_ASSERT(data_.size() == n_elements * sizeof(double));
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

    void read_data() {

      // Set the warning handler
      TIFFSetWarningHandler(&detail::tiff_warning_handler);

      // Open the TIFF file
      TIFF *tiff = TIFFOpen(filename_.c_str(), "r");
      DXTBX_ASSERT(tiff);

      // Get the image length
      uint32 slow_size = 0;
      uint32 fast_size = 0;
      uint32 samples_per_pixel = 1;
      uint32 bits_per_sample = 1;
      uint32 sample_format = SAMPLEFORMAT_UINT;
      TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &slow_size);
      TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH,  &fast_size);
      TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL,  &samples_per_pixel);
      TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE,  &bits_per_sample);
      TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT,  &sample_format);
      DXTBX_ASSERT(samples_per_pixel == 1);
      size_[0] = slow_size;
      size_[1] = fast_size;

      // Check the sample format
      if (sample_format == SAMPLEFORMAT_INT) {
        if (bits_per_sample == 16) {
          type_ = "int16";
        } else if (bits_per_sample == 32) {
          type_ = "int32";
        } else {
          DXTBX_ERROR("Unsupported bits per sample");
        }
      } else if (sample_format == SAMPLEFORMAT_UINT) {
        if (bits_per_sample == 16) {
          type_ = "uint16";
        } else if (bits_per_sample == 32) {
          type_ = "uint32";
        } else {
          DXTBX_ERROR("Unsupported bits per sample");
        }
      } else if (sample_format == SAMPLEFORMAT_IEEEFP) {
        if (bits_per_sample == 32) {
          type_ = "float32";
        } else if (bits_per_sample == 64) {
          type_ = "float64";
        } else {
          DXTBX_ERROR("Unsupported bits per sample");
        }
      } else {
        DXTBX_ERROR("Unsupported format");
      }

      // Get the scan line length and allocate a buffer
      tsize_t scanline_size = TIFFScanlineSize(tiff);
      DXTBX_ASSERT(scanline_size == (bits_per_sample / 8) * fast_size);

      // Allocate the image array
      data_.resize(scanline_size * slow_size);

      // Loop through and read the data
      for (std::size_t j = 0; j < slow_size; ++j) {
        DXTBX_ASSERT(TIFFReadScanline(tiff, (void *)&data_[j * scanline_size], j) == 1);
      }

      // Close the tiff file
      TIFFClose(tiff);
    }

    std::vector<char> data_;
  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_TIFF_H
