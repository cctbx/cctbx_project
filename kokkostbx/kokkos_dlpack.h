#ifndef KOKKOS_DLPACK_H
#define KOKKOS_DLPACK_H
#include <Kokkos_Core.hpp>
#include <dlpack/dlpack.h>

namespace kokkostbx {

template<typename DataType, typename SpaceType>
DLDataTypeCode getDLPackTypeCode() {
  using ValueType = typename Kokkos::View<DataType, SpaceType>::value_type;
  if (std::is_same<ValueType, float>::value) {
    return kDLFloat;
  } else if (std::is_same<ValueType, double>::value) {
    return kDLFloat;
  } else if (std::is_same<ValueType, int>::value) {
    return kDLInt;
  } else if (std::is_same<ValueType, unsigned int>::value) {
    return kDLUInt;
  // } else if (std::is_same<ValueType, bool>::value) {
    // return kDLBool;
  } else {
    // Unsupported data type
    throw std::runtime_error("Unsupported data type for DLPack conversion");
  }
}

template<typename DataType, typename SpaceType>
DLDataType getDLPackDataType() {
  DLDataType dtype;
  dtype.code = getDLPackTypeCode<DataType, SpaceType>();
  dtype.bits = sizeof(typename Kokkos::View<DataType, SpaceType>::value_type) * 8;
  dtype.lanes = 1;
  return dtype;
}

template<typename SpaceType>
DLDevice getDLPackDevice() {
  DLDevice dl_device;
  if (std::is_same<SpaceType, Kokkos::HostSpace>::value) {
      dl_device = {kDLCPU, 0};
  }
#ifdef KOKKOS_ENABLE_CUDA
  else if (std::is_same<SpaceType, Kokkos::CudaSpace>::value) {
      dl_device = {kDLCUDA, 0};
  } else if (std::is_same<SpaceType, Kokkos::CudaUVMSpace>::value) {
      dl_device = {kDLCUDAManaged, 0};
  } else if (std::is_same<SpaceType, Kokkos::CudaHostPinnedSpace>::value) {
      dl_device = {kDLCUDAHost, 0};
  }
#endif
#ifdef KOKKOS_ENABLE_HIP
  else if (std::is_same<SpaceType, Kokkos::HIPSpace>::value) {
      dl_device = {kDLROCM, 0};
  } else if (std::is_same<SpaceType, Kokkos::HIPHostPinnedSpace>::value) {
      dl_device = {kDLROCMHost, 0};
  } 
#endif
  else {
      // Extend to other device types as needed
      throw std::runtime_error("Unsupported Kokkos device type for DLPack conversion.");
  }
  return dl_device;
}

template<typename DataType, typename SpaceType>
DLManagedTensor* view_to_dlpack(Kokkos::View<DataType, SpaceType>& view) {
  // Get the Kokkos view size and dimensions
  constexpr size_t numDims = Kokkos::View<DataType, SpaceType>::rank;
  int64_t* shape = new int64_t[numDims];
  for (size_t i = 0; i < numDims; i++) {
    shape[i] = view.extent(i);
  }

  // Create a DLPack tensor
  DLManagedTensor* dlpackTensor = new DLManagedTensor;
  dlpackTensor->dl_tensor.data = view.data();
  // dlpackTensor->dl_tensor.device.device_type = DLDeviceType::kDLCPU;
  // dlpackTensor->dl_tensor.device.device_id = 0;
  dlpackTensor->dl_tensor.device = getDLPackDevice<SpaceType>();
  dlpackTensor->dl_tensor.ndim = numDims;    
  dlpackTensor->dl_tensor.dtype = getDLPackDataType<DataType, SpaceType>();
  dlpackTensor->dl_tensor.shape = shape;
  dlpackTensor->dl_tensor.strides = nullptr;
  dlpackTensor->dl_tensor.byte_offset = 0;
  dlpackTensor->manager_ctx = nullptr;
  dlpackTensor->deleter = [](DLManagedTensor* tensor) {
      std::cout << "Blob" << std::endl;
      delete[] tensor->dl_tensor.shape;
  };
  return dlpackTensor;
}

// template <typename ViewType>
// class KokkosViewToDLPack {
// public:
//   KokkosViewToDLPack(ViewType view) : view_(view) {}

//   torch::Tensor convertToDLPack() {
//     // Convert the Kokkos view to DLPack
//     DLManagedTensor* dlpackTensor = convertToDLPack();

//     // Convert the DLPack tensor to PyTorch
//     torch::Tensor tensor = torch::from_dlpack(dlpackTensor);

//     // Free the DLPack tensor memory
//     delete[] dlpackTensor->dl_tensor.shape;
//     delete dlpackTensor;

//     return tensor;
//   }
  
// private:
//   ViewType view_;

//   DLManagedTensor* convertToDLPack() {
//     // Get the Kokkos view size and dimensions
//     size_t numDims = ViewType::rank;
//     size_t* shape = new size_t[numDims];
//     for (size_t i = 0; i < numDims; i++) {
//       shape[i] = view_.extent(i);
//     }

//     // Create a DLPack tensor
//     DLManagedTensor* dlpackTensor = new DLManagedTensor;
//     dlpackTensor->dl_tensor.data = view_.data();
//     dlpackTensor->dl_tensor.ctx = const_cast<void*>(view_.impl_map().template device_data<void>());
//     dlpackTensor->dl_tensor.ndim = numDims;
//     dlpackTensor->dl_tensor.dtype = getDLPackDataType();
//     dlpackTensor->dl_tensor.shape = shape;
//     dlpackTensor->dl_tensor.strides = nullptr;
//     dlpackTensor->dl_tensor.byte_offset = 0;
//     dlpackTensor->manager_ctx = nullptr;
//     dlpackTensor->deleter = [](DLManagedTensor* tensor) { delete[] tensor->dl_tensor.shape; };

//     return dlpackTensor;
//   }

//   DLDataType getDLPackDataType() {
//     DLDataType dtype;
//     dtype.code = getDLPackTypeCode();
//     dtype.bits = sizeof(typename ViewType::value_type) * 8;
//     dtype.lanes = 1;
//     return dtype;
//   }

//   DLDataTypeCode getDLPackTypeCode() {
//     using ValueType = typename ViewType::value_type;
//     if (std::is_same<ValueType, float>::value) {
//       return kDLFloat;
//     } else if (std::is_same<ValueType, double>::value) {
//       return kDLFloat;
//     } else if (std::is_same<ValueType, int>::value) {
//       return kDLInt;
//     } else if (std::is_same<ValueType, unsigned int>::value) {
//       return kDLUInt;
//     } else if (std::is_same<ValueType, bool>::value) {
//       return kDLBool;
//     } else {
//       // Unsupported data type
//       throw std::runtime_error("Unsupported data type for DLPack conversion");
//     }
//   }
// };

}

#endif  // KOKKOS_DLPACK_H
