/*
 * tiff_reader.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_TIFF_READER_H
#define DXTBX_FORMAT_TIFF_READER_H

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

    /**
     * Construct the class with the filename
     */
    TIFFReader(const char *filename)
      : ImageReader(filename),
        slow_size_(0),
        fast_size_(0) {
      read_data();
    }

  protected:

    // Helper structs to get types
    template <typename T>
    struct array_type {
      typedef scitbx::af::versa< T, scitbx::af::c_grid<2> > type;
    };

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
      slow_size_ = slow_size;
      fast_size_ = fast_size;

      // Check the sample format
      if (sample_format == SAMPLEFORMAT_INT) {
        if (bits_per_sample == 16) {
          DXTBX_ASSERT(8 * sizeof(short) == (bits_per_sample));
          read_data_detail<short>(tiff);
        } else if (bits_per_sample == 32) {
          DXTBX_ASSERT(8 * sizeof(int) == (bits_per_sample));
          read_data_detail<int>(tiff);
        } else {
          DXTBX_ERROR("Unsupported bits per sample");
        }
      } else if (sample_format == SAMPLEFORMAT_UINT) {
        if (bits_per_sample == 16) {
          DXTBX_ASSERT(8 * sizeof(unsigned short) == (bits_per_sample));
          read_data_detail<unsigned short>(tiff);
        } else if (bits_per_sample == 32) {
          DXTBX_ASSERT(8 * sizeof(unsigned int) == (bits_per_sample));
          read_data_detail<unsigned int>(tiff);
        } else {
          DXTBX_ERROR("Unsupported bits per sample");
        }
      } else if (sample_format == SAMPLEFORMAT_IEEEFP) {
        if (bits_per_sample == 32) {
          DXTBX_ASSERT(8 * sizeof(float) == (bits_per_sample));
          read_data_detail<float>(tiff);
        } else if (bits_per_sample == 64) {
          DXTBX_ASSERT(8 * sizeof(double) == (bits_per_sample));
          read_data_detail<double>(tiff);
        } else {
          DXTBX_ERROR("Unsupported bits per sample");
        }
      } else {
        DXTBX_ERROR("Unsupported format");
      }

      // Close the tiff file
      TIFFClose(tiff);
    }

    template <typename T>
    void read_data_detail(TIFF *tiff) {

      typedef typename array_type<T>::type array_data_type;

      // Get the scan line length and allocate a buffer
      tsize_t scanline_size = TIFFScanlineSize(tiff);
      std::size_t element_size = sizeof(T);
      DXTBX_ASSERT(scanline_size == element_size * fast_size_);

      // The image grid
      DXTBX_ASSERT(slow_size_ > 0);
      DXTBX_ASSERT(fast_size_ > 0);
      scitbx::af::c_grid<2> grid(slow_size_, fast_size_);

      // Allocate the array
      array_data_type data(grid);

      // Loop through and read the data
      for (std::size_t j = 0; j < slow_size_; ++j) {
        DXTBX_ASSERT(TIFFReadScanline(tiff, (void *)&data[j * fast_size_], j) == 1);
      }

      // Add to the tiles list
      tiles_.push_back(variant_type(data));
      names_.push_back("");
    }

    std::size_t slow_size_;
    std::size_t fast_size_;

  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_TIFF_READER_H
