/*
 * hdf5_reader.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_FORMAT_HDF5_READER_H
#define DXTBX_FORMAT_HDF5_READER_H

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <include/cbf.h>
#include <dxtbx/format/image.h>
#include <dxtbx/error.h>
#include <hdf5.h>

namespace dxtbx { namespace format {

  /**
   * A class to read a HDF5 Image
   */
  class HDF5Reader : public MultiImageReader {
  public:

    /**
     * Construct the class with the filename
     */
    HDF5Reader(
        hid_t handle,
        const scitbx::af::const_ref<std::string> &datasets)
      : handle_(handle),
        datasets_(datasets.begin(), datasets.end()) {
      generate_lookup();
      first_ = 0;
      last_ = lookup_.size();
    }

    /**
     * Construct the class with the filename
     */
    HDF5Reader(
        hid_t handle,
        const scitbx::af::const_ref<std::string> &datasets,
        std::size_t first,
        std::size_t last)
      : handle_(handle),
        datasets_(datasets.begin(), datasets.end()),
        first_(first),
        last_(last) {
      DXTBX_ASSERT(first < last);
      generate_lookup();
      DXTBX_ASSERT(last <= lookup_.size());
    }

    /**
     * Return the filename
     */
    std::string filename() const {
      char buffer[1024];
      int n = H5Fget_name(handle_, buffer, 1024);
      DXTBX_ASSERT(n > 0);
      return std::string(buffer, n);
    }

    /**
     * Get the file handle
     */
    hid_t handle() const {
      return handle_;
    }

    /**
     * Return the dataset list
     */
    scitbx::af::shared<std::string> datasets() const {
      return datasets_;
    }

    /**
     * Return the first index
     */
    std::size_t first() const {
      return first_;
    }

    /**
     * Return the last index
     */
    std::size_t last() const {
      return last_;
    }

    /**
     * Return the number of images
     */
    std::size_t size() const {
      return last_ - first_;
    }

    /**
     * Get the image at the given index
     */
    Image image(std::size_t index) const {
      DXTBX_ASSERT(index < size());
      Image result;
      result.push_back(read_data(index), "");
      return result;
    }

  protected:

    /**
     * An item in the lookup list
     */
    struct Item {
      std::string dataset;
      std::size_t index;
      Item(std::string d, std::size_t i)
        : dataset(d),
          index(i) {}
    };

    /**
     * Generate a list of dataset, index pairs that correspond to the location
     * of the image at each index. A HDF5 file can have data in multiple
     * datasets
     */
    void generate_lookup() {
      for (std::size_t i = 0; i < datasets_.size(); ++i) {
        std::size_t n = num_images(datasets_[i].c_str());
        for (std::size_t j = 0; j < n; ++j) {
          lookup_.push_back(Item(datasets_[i], j));
        }
      }
    }

    /**
     * Return the number of images in a dataset
     */
    std::size_t num_images(const char *dataset) const {
      hid_t dataset_id = H5Dopen(handle_, dataset, H5P_DEFAULT);
      hid_t file_space_id = H5Dget_space(dataset_id);
      std::size_t ndims = H5Sget_simple_extent_ndims(file_space_id);
      DXTBX_ASSERT(ndims == 3);
      std::vector<hsize_t> dataset_dims(ndims);
      H5Sget_simple_extent_dims(file_space_id, &dataset_dims[0], NULL);
      H5Dclose(dataset_id);
      H5Sclose(file_space_id);
      return dataset_dims[0];
    }

    /**
     * Read the data at the given index
     */
    scitbx::af::versa< int, scitbx::af::c_grid<2> > read_data(std::size_t index) const {
      Item item = lookup_[first_ + index];
      return read_data_detail(item.dataset.c_str(), item.index);
    }

    /**
     * Read the data in the dataset at the given index
     */
    scitbx::af::versa< int, scitbx::af::c_grid<2> > read_data_detail(
        const char *dataset, std::size_t index) const {

      // Get the file space
      hid_t dataset_id = H5Dopen(handle_, dataset, H5P_DEFAULT);
      hid_t file_space_id = H5Dget_space(dataset_id);
      std::size_t ndims = H5Sget_simple_extent_ndims(file_space_id);
      DXTBX_ASSERT(ndims == 3);
      std::vector<hsize_t> dataset_dims(ndims);
      H5Sget_simple_extent_dims(file_space_id, &dataset_dims[0], NULL);

      // Create the grid
      std::vector<hsize_t> start(ndims);
      std::vector<hsize_t> count(ndims);
      start[0] = index;
      start[1] = 0;
      start[2] = 0;
      count[0] = 1;
      count[1] = dataset_dims[1];
      count[2] = dataset_dims[2];
      DXTBX_ASSERT(index < dataset_dims[0]);

      // Create the data array
      scitbx::af::flex_grid<> grid(dataset_dims[1], dataset_dims[2]);
      scitbx::af::versa<int, scitbx::af::c_grid<2> > data(grid);

      // Create the dataspace id
      herr_t status1 = H5Sselect_hyperslab(
          file_space_id,
          H5S_SELECT_SET,
          &start[0],
          NULL,
          &count[0],
          NULL);
      DXTBX_ASSERT(status1 >= 0);

      // Create the memory space size
      hid_t mem_space_id = H5Screate_simple(ndims, &count[0], NULL);

      // Copy the data
      herr_t status2 = H5Dread(
          dataset_id,
          H5T_NATIVE_INT,
          mem_space_id,
          file_space_id,
          H5P_DEFAULT,
          &data[0]);
      DXTBX_ASSERT(status2 >= 0);

      // Close some stuff
      H5Sclose(mem_space_id);
      H5Sclose(file_space_id);
      H5Dclose(dataset_id);

      // Return the data
      return data;
    }

    hid_t handle_;
    scitbx::af::shared<std::string> datasets_;
    std::size_t first_;
    std::size_t last_;
    std::vector<Item> lookup_;

  };


}} // namespace dxtbx::format

#endif // DXTBX_FORMAT_HDF5_READER_H
