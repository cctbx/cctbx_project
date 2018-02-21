#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/shared_ptr.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <vector>
#include <dxtbx/imageset.h>
#include <dxtbx/model/pixel_to_millimeter.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace boost_python {

  using model::OffsetParallaxCorrectedPxMmStrategy;
  using model::OffsetPxMmStrategy;

  namespace detail {

    /**
     * Pickle a python object to a string
     */
    std::string pickle_dumps(boost::python::object x) {
      boost::python::object main = boost::python::import("__main__");
      boost::python::object global(main.attr("__dict__"));
      boost::python::object result = exec(
          "def dumps(x):import pickle; return pickle.dumps(x)",
          global, global);
      boost::python::object dumps = global["dumps"];
      return boost::python::extract<std::string>(dumps(x))();
    }

    /**
     * Unpickle a python object from a string
     */
    boost::python::object pickle_loads(std::string x) {
      if (x == "") {
        return boost::python::object();
      }
      boost::python::object main = boost::python::import("__main__");
      boost::python::object global(main.attr("__dict__"));
      boost::python::object result = exec(
          "def loads(x):import pickle; return pickle.loads(x)",
          global, global);
      boost::python::object loads = global["loads"];
      return loads(x);
    }

    /**
     * Tuple from list
     */
    boost::python::tuple list_to_tuple(boost::python::list x) {
      boost::python::object main = boost::python::import("__main__");
      boost::python::object global(main.attr("__dict__"));
      boost::python::object result = exec(
          "def func_list_to_tuple(x):return tuple(x)",
          global, global);
      boost::python::object func_list_to_tuple = global["func_list_to_tuple"];
      return boost::python::extract<boost::python::tuple>(func_list_to_tuple(x))();
    }
  }

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData> make_imageset_data(
      boost::python::object reader,
      boost::python::object masker,
      std::string filename_template,
      std::string vendor,
      boost::python::dict params,
      boost::python::object format) {

    // Create the pointer
    boost::shared_ptr<ImageSetData> self(new ImageSetData(reader, masker));

    // Set some stuff
    self->set_template(filename_template);
    self->set_vendor(vendor);
    self->set_params(detail::pickle_dumps(params));
    self->set_format(detail::pickle_dumps(format));

    // Return the imageset data
    return self;
  }

  /**
   * Set the parameters
   */
  boost::python::object ImageSetData_get_params(ImageSetData &self) {
    return detail::pickle_loads(self.get_params());
  }

  /**
   * Get the parameters
   */
  void ImageSetData_set_params(ImageSetData &self, boost::python::dict params) {
    self.set_params(detail::pickle_dumps(params));
  }

  /**
   * Set the format class
   */
  boost::python::object ImageSetData_get_format(ImageSetData &self) {
    return detail::pickle_loads(self.get_format());
  }

  /**
   * Get the format class
   */
  void ImageSetData_set_format(ImageSetData &self, boost::python::dict format) {
    self.set_format(detail::pickle_dumps(format));
  }


  /**
   * A constructor for the imageset class
   */
  boost::shared_ptr<ImageSet> make_imageset(
      const ImageSetData &data,
      boost::python::object indices) {

    if (indices == boost::python::object()) {
      return boost::shared_ptr<ImageSet>(new ImageSet(data));
    }

    return boost::shared_ptr<ImageSet>(
        new ImageSet(
          data,
          boost::python::extract<
            scitbx::af::const_ref<std::size_t> >(indices)()));
  }

  /**
   * Implement pickling for ImageSetData class
   */
  struct ImageSetDataPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getinitargs(ImageSetData obj) {
      return boost::python::make_tuple(
          obj.reader(),
          obj.masker());
    }

    static
    boost::shared_ptr<BeamBase> get_beam(const ImageSetData &self, std::size_t i) {
      return self.get_beam(i);
    }

    static
    boost::shared_ptr<Detector> get_detector(const ImageSetData &self, std::size_t i) {
      return self.get_detector(i);
    }

    static
    boost::shared_ptr<Goniometer> get_goniometer(const ImageSetData &self, std::size_t i) {
      return self.get_goniometer(i);
    }

    static
    boost::shared_ptr<Scan> get_scan(const ImageSetData &self, std::size_t i) {
      return self.get_scan(i);
    }

    template <typename Model, typename Func>
    static
    boost::python::tuple get_model_list(ImageSetData obj, Func get) {

      // Create a list of models and a list of indices
      std::vector< boost::shared_ptr<Model> > model_list;
      std::vector<std::size_t> index_list;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        boost::shared_ptr<Model> m = get(obj, i);
        std::size_t k = model_list.size();
        for (std::size_t j = 0; j < k; ++j) {
          if (m.get() == model_list[j].get()) {
            k = j;
            break;
          }
        }
        if (k == model_list.size()) {
          model_list.push_back(m);
        }
        index_list.push_back(k);
      }

      // Convert to python lists and return tuple
      boost::python::list models;
      boost::python::list indices;
      for (std::size_t i = 0; i < model_list.size(); ++i) {
        models.append(model_list[i]);
      }
      for (std::size_t i = 0; i < index_list.size(); ++i) {
        indices.append(index_list[i]);
      }
      return boost::python::make_tuple(models, indices);
    }

    static
    boost::python::tuple get_model_tuple(ImageSetData obj) {
      return boost::python::make_tuple(
          ImageSetDataPickleSuite::get_model_list<BeamBase>(
            obj, &ImageSetDataPickleSuite::get_beam),
          ImageSetDataPickleSuite::get_model_list<Detector>(
            obj, &ImageSetDataPickleSuite::get_detector),
          ImageSetDataPickleSuite::get_model_list<Goniometer>(
            obj, &ImageSetDataPickleSuite::get_goniometer),
          ImageSetDataPickleSuite::get_model_list<Scan>(
            obj, &ImageSetDataPickleSuite::get_scan));
    }

    static
    boost::python::tuple get_lookup_tuple(ImageSetData obj) {
      return boost::python::make_tuple(
          boost::python::make_tuple(
            obj.external_lookup().mask().get_filename(),
            obj.external_lookup().mask().get_data()),
          boost::python::make_tuple(
            obj.external_lookup().gain().get_filename(),
            obj.external_lookup().gain().get_data()),
          boost::python::make_tuple(
            obj.external_lookup().pedestal().get_filename(),
            obj.external_lookup().pedestal().get_data()),
          boost::python::make_tuple(
            obj.external_lookup().dx().get_filename(),
            obj.external_lookup().dx().get_data()),
          boost::python::make_tuple(
            obj.external_lookup().dy().get_filename(),
            obj.external_lookup().dy().get_data()));
    }

    static
    boost::python::tuple getstate(ImageSetData obj) {
      return boost::python::make_tuple(
          ImageSetDataPickleSuite::get_model_tuple(obj),
          ImageSetDataPickleSuite::get_lookup_tuple(obj),
          obj.get_template(),
          obj.get_vendor(),
          obj.get_params(),
          obj.get_format());
    }

    template <typename Model, typename Func>
    static
    void set_model_list(ImageSetData &obj, boost::python::tuple data, Func set) {

      // Extract to python lists
      boost::python::list models = boost::python::extract<
        boost::python::list>(data[0])();
      boost::python::list indices = boost::python::extract<
        boost::python::list>(data[1])();

      // Convert to c++ vectors
      std::vector< boost::shared_ptr<Model> > model_list;
      std::vector< std::size_t > index_list;
      for (std::size_t i = 0; i < boost::python::len(models); ++i) {
        model_list.push_back(
            boost::python::extract<
              boost::shared_ptr<Model> >(models[i])());
      }
      for (std::size_t i = 0; i < boost::python::len(indices); ++i) {
        index_list.push_back(boost::python::extract<std::size_t>(indices[i])());
      }

      // Set the models
      DXTBX_ASSERT(index_list.size() == obj.size());
      for (std::size_t i = 0; i < index_list.size(); ++i) {
        DXTBX_ASSERT(index_list[i] < model_list.size());
        ((&obj)->*set)(model_list[index_list[i]], i);
      }
    }

    static
    void set_model_tuple(ImageSetData &obj, boost::python::tuple models) {
      DXTBX_ASSERT(boost::python::len(models) == 4);
      ImageSetDataPickleSuite::set_model_list<BeamBase>(
          obj,
          boost::python::extract<boost::python::tuple>(models[0])(),
          &ImageSetData::set_beam);
      ImageSetDataPickleSuite::set_model_list<Detector>(
          obj,
          boost::python::extract<boost::python::tuple>(models[1]),
          &ImageSetData::set_detector);
      ImageSetDataPickleSuite::set_model_list<Goniometer>(
          obj,
          boost::python::extract<boost::python::tuple>(models[2]),
          &ImageSetData::set_goniometer);
      ImageSetDataPickleSuite::set_model_list<Scan>(
          obj,
          boost::python::extract<boost::python::tuple>(models[3]),
          &ImageSetData::set_scan);
    }

    template <typename Data, typename Func>
    static
    void set_lookup_item(ImageSetData &obj, boost::python::tuple lookup, Func item) {
      DXTBX_ASSERT(boost::python::len(lookup) == 2);

      // Extract
      std::string filename = boost::python::extract<std::string>(lookup[0])();
      Data data = boost::python::extract<Data>(lookup[1])();

      // Set the filename and data
      ((&obj.external_lookup())->*item)().set_filename(filename);
      ((&obj.external_lookup())->*item)().set_data(data);
    }

    static
    void set_lookup_tuple(ImageSetData &obj, boost::python::tuple lookup) {
      DXTBX_ASSERT(boost::python::len(lookup) == 5);
      ImageSetDataPickleSuite::set_lookup_item< Image<bool> >(
          obj,
          boost::python::extract<boost::python::tuple>(lookup[0])(),
          &ExternalLookup::mask);
      ImageSetDataPickleSuite::set_lookup_item< Image<double> >(
          obj,
          boost::python::extract<boost::python::tuple>(lookup[1])(),
          &ExternalLookup::gain);
      ImageSetDataPickleSuite::set_lookup_item< Image<double> >(
          obj,
          boost::python::extract<boost::python::tuple>(lookup[2])(),
          &ExternalLookup::pedestal);
      ImageSetDataPickleSuite::set_lookup_item< Image<double> >(
          obj,
          boost::python::extract<boost::python::tuple>(lookup[3])(),
          &ExternalLookup::dx);
      ImageSetDataPickleSuite::set_lookup_item< Image<double> >(
          obj,
          boost::python::extract<boost::python::tuple>(lookup[4])(),
          &ExternalLookup::dy);
    }

    static
    void setstate(ImageSetData &obj, boost::python::tuple state) {

      DXTBX_ASSERT(boost::python::len(state) == 6);

      // Set the models
      ImageSetDataPickleSuite::set_model_tuple(
          obj,
          boost::python::extract<boost::python::tuple>(state[0])());

      // Set the lookup
      ImageSetDataPickleSuite::set_lookup_tuple(
          obj,
          boost::python::extract<boost::python::tuple>(state[1])());

      // Set the properties
      obj.set_template(boost::python::extract<std::string>(state[2])());
      obj.set_vendor(boost::python::extract<std::string>(state[3])());
      obj.set_params(boost::python::extract<std::string>(state[4])());
      obj.set_format(boost::python::extract<std::string>(state[5])());
    }
  };

  /**
   * Implement pickling for ImageSet class
   */
  struct ImageSetPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getinitargs(ImageSet obj) {
      return boost::python::make_tuple(
          obj.data(),
          obj.indices());
    }

  };

  /**
   * Implement pickling for ImageGrid class
   */
  struct ImageGridPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getinitargs(ImageGrid obj) {
      return boost::python::make_tuple(
          obj.data(),
          obj.indices(),
          obj.get_grid_size());
    }

  };

  /**
   * Implement pickling for ImageSweep class
   */
  struct ImageSweepPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getinitargs(ImageSweep obj) {
      return boost::python::make_tuple(
          obj.data(),
          obj.indices(),
          obj.get_beam(),
          obj.get_detector(),
          obj.get_goniometer(),
          obj.get_scan());
    }

  };

  /**
   * Get the external lookup data as a tuple of images
   */
  template<typename T>
  boost::python::object ExternalLookupItem_get_data(
      const ExternalLookupItem<T> &obj) {

    // Get the image data
    Image<T> data = obj.get_data();

    // If empty then return None
    if (data.empty()) {
      return boost::python::object();
    }

    // Otherwise, put image data in a list
    boost::python::list result;
    for (std::size_t i = 0; i < data.n_tiles(); ++i) {
      result.append(data.tile(i).data());
    }

    // Return the image
    return result;
  }

  /**
   * Set the external lookup data
   */
  template <typename T>
  void ExternalLookupItem_set_data(
      ExternalLookupItem<T> &obj,
      boost::python::object item) {

    typedef typename scitbx::af::flex<T>::type flex_type;

    Image<T> data;

    // If item is not None
    if (item != boost::python::object()) {

      std::string name = boost::python::extract<std::string>(
          item.attr("__class__").attr("__name__"))();

      if (name == "tuple") {

        // If we have a tuple then add items of tuple to image data
        for (std::size_t i = 0; i < boost::python::len(item); ++i) {
          flex_type a = boost::python::extract<flex_type>(item)();
          data.push_back(
            ImageTile<T>(
              scitbx::af::versa<T, scitbx::af::c_grid<2> >(
                a.handle(),
                scitbx::af::c_grid<2>(a.accessor()))));
        }
      } else {

        try {

          // If we have a single array then add
          flex_type a = boost::python::extract<flex_type>(item)();
          data.push_back(
            ImageTile<T>(
              scitbx::af::versa<T, scitbx::af::c_grid<2> >(
                a.handle(),
                scitbx::af::c_grid<2>(a.accessor()))));

        } catch (boost::python::error_already_set) {

          data = boost::python::extract< Image<T> >(item)();
          boost::python::handle_exception();
        }
      }
    }

    // Set the image data
    obj.set_data(data);
  }


  template <typename T>
  boost::python::tuple image_as_tuple(const Image<T> &image) {
    boost::python::list result;
    for (std::size_t i = 0; i < image.n_tiles(); ++i) {
      result.append(image.tile(i).data());
    }
    return detail::list_to_tuple(result);
  }

  boost::python::tuple ImageSet_get_raw_data(ImageSet &self, std::size_t index) {
    boost::python::tuple result;
    ImageBuffer buffer = self.get_raw_data(index);
    if (buffer.is_int()) {
      result = image_as_tuple<int>(buffer.as_int());
    } else if (buffer.is_double()) {
      result = image_as_tuple<double>(buffer.as_double());
    } else {
      throw DXTBX_ERROR("Problem reading raw data");
    }
    return result;
  }

  boost::python::tuple ImageSet_get_corrected_data(ImageSet &self, std::size_t index) {
    return image_as_tuple<double>(self.get_corrected_data(index));
  }

  boost::python::tuple ImageSet_get_gain(ImageSet &self, std::size_t index) {
    return image_as_tuple<double>(self.get_gain(index));
  }

  boost::python::tuple ImageSet_get_pedestal(ImageSet &self, std::size_t index) {
    return image_as_tuple<double>(self.get_pedestal(index));
  }

  boost::python::tuple ImageSet_get_mask(ImageSet &self, std::size_t index) {
    return image_as_tuple<bool>(self.get_mask(index));
  }

  /**
   * Wrapper for the external lookup items
   */
  template<typename T>
  void external_lookup_item_wrapper(const char *name) {
    using namespace boost::python;

    class_< ExternalLookupItem<T> >(name)
      .add_property("filename",
          &ExternalLookupItem<T>::get_filename,
          &ExternalLookupItem<T>::set_filename)
      .add_property("data",
          &ExternalLookupItem<T>::get_data,
          &ExternalLookupItem<T>::set_data)
      ;

  }


  /**
   * If we have offset arrays set in the imageset then update the pixel to
   * millimeter strategy to use them
   */
  void ImageSet_update_detector_px_mm_data(ImageSet &self) {
    Image<double> dx = self.external_lookup().dx().get_data();
    Image<double> dy = self.external_lookup().dy().get_data();
    DXTBX_ASSERT(dx.empty() == dy.empty());
    if (dx.empty() && dy.empty()) {
      return;
    }
    for (std::size_t i = 0; i < self.size(); ++i) {
      ImageSet::detector_ptr detector = self.get_detector_for_image(i);
      DXTBX_ASSERT(dx.n_tiles() == detector->size());
      DXTBX_ASSERT(dy.n_tiles() == detector->size());
      for (std::size_t i = 0; i < detector->size(); ++i) {
        Panel& panel = detector->operator[](i);
        if (panel.get_px_mm_strategy()->name() == "ParallaxCorrectedPxMmStrategy" ||
            panel.get_px_mm_strategy()->name() == "OffsetParallaxCorrectedPxMmStrategy") {
          boost::shared_ptr<OffsetParallaxCorrectedPxMmStrategy> strategy = boost::make_shared<
            OffsetParallaxCorrectedPxMmStrategy>(
              panel.get_mu(),
              panel.get_thickness(),
              dx.tile(i).data(),
              dy.tile(i).data());
          panel.set_px_mm_strategy(strategy);
        } else if (
            panel.get_px_mm_strategy()->name() == "SimplePxMmStrategy" ||
            panel.get_px_mm_strategy()->name() == "OffsetPxMmStrategy") {
          boost::shared_ptr<OffsetPxMmStrategy> strategy = boost::make_shared<OffsetPxMmStrategy>(
            dx.tile(i).data(),
            dy.tile(i).data());
          panel.set_px_mm_strategy(strategy);
        }
      }

    }
  }

  /**
   * If we have offset arrays set in the imageset then update the pixel to
   * millimeter strategy to use them
   */
  void ImageSweep_update_detector_px_mm_data(ImageSweep &self) {
    ImageSweep::detector_ptr detector = self.get_detector();
    Image<double> dx = self.external_lookup().dx().get_data();
    Image<double> dy = self.external_lookup().dy().get_data();
    DXTBX_ASSERT(dx.empty() == dy.empty());
    if (dx.empty() && dy.empty()) {
      return;
    }
    DXTBX_ASSERT(dx.n_tiles() == detector->size());
    DXTBX_ASSERT(dy.n_tiles() == detector->size());
    for (std::size_t i = 0; i < detector->size(); ++i) {
      Panel& panel = detector->operator[](i);
      if (panel.get_px_mm_strategy()->name() == "ParallaxCorrectedPxMmStrategy" ||
          panel.get_px_mm_strategy()->name() == "OffsetParallaxCorrectedPxMmStrategy") {
        boost::shared_ptr<OffsetParallaxCorrectedPxMmStrategy> strategy = boost::make_shared<
          OffsetParallaxCorrectedPxMmStrategy>(
          panel.get_mu(),
          panel.get_thickness(),
          dx.tile(i).data(),
          dy.tile(i).data());
        panel.set_px_mm_strategy(strategy);
      } else if (
          panel.get_px_mm_strategy()->name() == "SimplePxMmStrategy" ||
          panel.get_px_mm_strategy()->name() == "OffsetPxMmStrategy") {
        boost::shared_ptr<OffsetPxMmStrategy> strategy = boost::make_shared<OffsetPxMmStrategy>(
          dx.tile(i).data(),
          dy.tile(i).data());
        panel.set_px_mm_strategy(strategy);
      }
    }
  }

  /**
   * Export the imageset classes
   */
  void export_imageset() {
    using namespace boost::python;

    external_lookup_item_wrapper<double>("ExternalLookupItemDouble");
    external_lookup_item_wrapper<bool>("ExternalLookupItemBool");

     class_<ExternalLookup>("ExternalLookup")
      .add_property("mask",
          make_function(
            &ExternalLookup::mask,
            return_internal_reference<>()))
      .add_property("gain",
          make_function(
            &ExternalLookup::gain,
            return_internal_reference<>()))
      .add_property("pedestal",
          make_function(
            &ExternalLookup::pedestal,
            return_internal_reference<>()))
      .add_property("dx",
          make_function(
            &ExternalLookup::dx,
            return_internal_reference<>()))
      .add_property("dy",
          make_function(
            &ExternalLookup::dy,
            return_internal_reference<>()))
      ;

    class_<ImageSetData, boost::shared_ptr<ImageSetData> >("ImageSetData", no_init)
      .def(init<
          boost::python::object,
          boost::python::object>())
      .def("__init__",
          make_constructor(
            &make_imageset_data,
            default_call_policies(), (
              arg("reader"),
              arg("masker"),
              arg("template") = "",
              arg("vendor") = "",
              arg("params") = boost::python::object(),
              arg("format") = boost::python::object())))
      .def("reader", &ImageSetData::reader)
      .def("masker", &ImageSetData::masker)
      .def("get_data", &ImageSetData::get_data)
      .def("get_mask", &ImageSetData::get_mask)
      .def("has_single_file_reader", &ImageSetData::has_single_file_reader)
      .def("get_path", &ImageSetData::get_path)
      .def("get_master_path", &ImageSetData::get_master_path)
      .def("get_image_identifier", &ImageSetData::get_image_identifier)
      .def("mark_for_rejection", &ImageSetData::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetData::is_marked_for_rejection)
      .def("get_beam", &ImageSetData::get_beam)
      .def("get_detector", &ImageSetData::get_detector)
      .def("get_goniometer", &ImageSetData::get_goniometer)
      .def("get_scan", &ImageSetData::get_scan)
      .def("set_beam", &ImageSetData::set_beam)
      .def("set_detector", &ImageSetData::set_detector)
      .def("set_goniometer", &ImageSetData::set_goniometer)
      .def("set_scan", &ImageSetData::set_scan)
      .def("get_template", &ImageSetData::get_template)
      .def("set_template", &ImageSetData::set_template)
      .def("get_vendor", &ImageSetData::get_vendor)
      .def("set_vendor", &ImageSetData::set_vendor)
      .def("get_params", &ImageSetData_get_params)
      .def("set_params", &ImageSetData_set_params)
      .def("get_format_class", &ImageSetData_get_format)
      .def("set_format_class", &ImageSetData_set_format)
      .add_property("external_lookup",
          make_function(
            &ImageSetData::external_lookup,
            return_internal_reference<>()))
      .def_pickle(ImageSetDataPickleSuite())
      ;

    class_<ImageSet>("ImageSet", no_init)
      .def("__init__",
          make_constructor(
            &make_imageset,
            default_call_policies(), (
              arg("data"),
              arg("indices") = boost::python::object())))
      .def("data", &ImageSet::data)
      .def("indices", &ImageSet::indices)
      .def("size", &ImageSet::size)
      .def("__len__", &ImageSet::size)
      .def("has_dynamic_mask", &ImageSet::has_dynamic_mask)
      .def("get_raw_data", &ImageSet_get_raw_data)
      .def("get_corrected_data", &ImageSet_get_corrected_data)
      .def("get_gain", &ImageSet_get_gain)
      .def("get_pedestal", &ImageSet_get_pedestal)
      .def("get_mask", &ImageSet_get_mask)
      .def("get_beam",
          &ImageSet::get_beam_for_image, (
            arg("index") = 0))
      .def("get_detector",
          &ImageSet::get_detector_for_image, (
            arg("index") = 0))
      .def("get_goniometer",
          &ImageSet::get_goniometer_for_image, (
            arg("index") = 0))
      .def("get_scan",
          &ImageSet::get_scan_for_image, (
            arg("index") = 0))
      .def("set_beam",
          &ImageSet::set_beam_for_image, (
            arg("index") = 0))
      .def("set_detector",
          &ImageSet::set_detector_for_image, (
            arg("index") = 0))
      .def("set_goniometer",
          &ImageSet::set_goniometer_for_image, (
            arg("index") = 0))
      .def("set_scan",
          &ImageSet::set_scan_for_image, (
            arg("index") = 0))
      .def("get_path", &ImageSet::get_path)
      .def("get_image_identifier", &ImageSet::get_image_identifier)
      .def("mark_for_rejection", &ImageSet::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSet::is_marked_for_rejection)
      .def("as_imageset", &ImageSet::as_imageset)
      .def("complete_set", &ImageSet::complete_set)
      .def("partial_set", &ImageSet::partial_set)
      .def("__eq__", &ImageSet::operator==)
      .def("__ne__", &ImageSet::operator!=)
      .def("update_detector_px_mm_data",
          &ImageSet_update_detector_px_mm_data)
      .add_property("external_lookup",
          make_function(
            &ImageSet::external_lookup,
            return_internal_reference<>()))
      .def_pickle(ImageSetPickleSuite())
      ;

    class_<ImageGrid, bases<ImageSet> >("ImageGrid", no_init)
      .def(init<
          const ImageSetData &,
          int2
          >((arg("data"),
             arg("grid_size"))))
      .def(init<
          const ImageSetData &,
          const scitbx::af::const_ref<std::size_t> &,
          int2
          >((arg("data"),
             arg("indices"),
             arg("grid_size"))))
      .def("get_grid_size", &ImageGrid::get_grid_size)
      .def("from_imageset", &ImageGrid::from_imageset)
      .def_pickle(ImageGridPickleSuite())
      ;

    class_<ImageSweep, bases<ImageSet> >("ImageSweep", no_init)
      .def(init<
          const ImageSetData &,
          const ImageSweep::beam_ptr &,
          const ImageSweep::detector_ptr &,
          const ImageSweep::goniometer_ptr &,
          const ImageSweep::scan_ptr &
          >((
              arg("data"),
              arg("beam"),
              arg("detector"),
              arg("goniometer"),
              arg("scan"))))
      .def(init<
          const ImageSetData &,
          const scitbx::af::const_ref<std::size_t> &,
          const ImageSweep::beam_ptr &,
          const ImageSweep::detector_ptr &,
          const ImageSweep::goniometer_ptr &,
          const ImageSweep::scan_ptr &
          >((
              arg("data"),
              arg("indices"),
              arg("beam"),
              arg("detector"),
              arg("goniometer"),
              arg("scan"))))
      .def("get_beam", &ImageSweep::get_beam_for_image)
      .def("get_detector", &ImageSweep::get_detector_for_image)
      .def("get_goniometer", &ImageSweep::get_goniometer_for_image)
      .def("get_scan", &ImageSweep::get_scan_for_image)
      .def("set_beam", &ImageSweep::set_beam_for_image)
      .def("set_detector", &ImageSweep::set_detector_for_image)
      .def("set_goniometer", &ImageSweep::set_goniometer_for_image)
      .def("set_scan", &ImageSweep::set_scan_for_image)
      .def("get_beam", &ImageSweep::get_beam)
      .def("get_detector", &ImageSweep::get_detector)
      .def("get_goniometer", &ImageSweep::get_goniometer)
      .def("get_scan", &ImageSweep::get_scan)
      .def("set_beam", &ImageSweep::set_beam)
      .def("set_detector", &ImageSweep::set_detector)
      .def("set_goniometer", &ImageSweep::set_goniometer)
      .def("set_scan", &ImageSweep::set_scan)
      .def("get_array_range", &ImageSweep::get_array_range)
      .def("complete_set", &ImageSweep::complete_sweep)
      .def("partial_set", &ImageSweep::partial_sweep)
      .def("update_detector_px_mm_data",
          &ImageSweep_update_detector_px_mm_data)
      .def_pickle(ImageSweepPickleSuite())
      ;
  }


  BOOST_PYTHON_MODULE(dxtbx_imageset_ext)
  {
    export_imageset();
  }

}} //namespace dxtbx::boost_python
