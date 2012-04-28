#include <boost/python.hpp>

#include <scitbx/array_family/shared.h>

#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>

namespace ipctbx {

namespace af = scitbx::af;

template <typename ElementType>
int
copy_to_new_segment(
  af::const_ref<ElementType> const& data)
{
  int shmid = shmget(
    IPC_PRIVATE,
    data.size() * sizeof(ElementType),
    IPC_CREAT|IPC_EXCL|S_IRUSR|S_IWUSR);
  if (shmid == -1){
    std::string msg = strerror(errno);
    throw std::runtime_error("shmget() failure: " + msg);
  }
  void* void_ptr = shmat(shmid, 0, 0);
  if (void_ptr == reinterpret_cast<void *>(-1)) {
    std::string msg = strerror(errno);
    throw std::runtime_error("shmat() failure: " + msg);
  }
  ElementType* data_ptr = static_cast<ElementType*>(void_ptr);
  std::copy(data.begin(), data.end(), data_ptr);
  return shmid;
}

template <typename ElementType>
af::shared<ElementType>
copy_from_segment(
  int shmid)
{
  void* void_ptr = shmat(shmid, 0, SHM_RDONLY);
  if (void_ptr == reinterpret_cast<void *>(-1)) {
    std::string msg = strerror(errno);
    throw std::runtime_error("shmat() failure: " + msg);
  }
  struct shmid_ds shm_stat;
  if (shmctl(shmid, IPC_STAT, &shm_stat) != 0) {
    std::string msg = strerror(errno);
    throw std::runtime_error("shmctl(IPC_STAT) failure: " + msg);
  }
  std::size_t data_size = shm_stat.shm_segsz / sizeof(ElementType);
  if (data_size * sizeof(ElementType) != shm_stat.shm_segsz) {
    throw std::runtime_error(
      "copy_from_segment() failure: unexpected segment size");
  }
  ElementType const* data_ptr = static_cast<ElementType const*>(void_ptr);
  return af::shared<ElementType>(data_ptr, data_ptr+data_size);
}

boost::python::list
list_segment_ids(
  int attached_status=0,
  bool ignore_errors=false)
{
  boost::python::list result;
  struct shmid_ds shm_info;
  int max_id = shmctl(0, SHM_INFO, &shm_info);
  if (max_id < 0) {
    if (ignore_errors) return result;
    std::string msg = strerror(errno);
    throw std::runtime_error("shmctl(SHM_INFO) failure: " + msg);
  }
  for (int i=0;i<=max_id;i++) {
    struct shmid_ds shm_stat;
    int shmid = shmctl(i, SHM_STAT, &shm_stat);
    if (shmid < 0) {
      if (i == 0 && max_id == 0) break;
      if (ignore_errors) continue;
      std::string msg = strerror(errno);
      throw std::runtime_error("shmctl(SHM_STAT) failure: " + msg);
    }
    if (attached_status == 0
        || ((attached_status < 0) == (shm_stat.shm_nattch == 0))) {
      result.append(shmid);
    }
  }
  return result;
}

bool
shmctl_rmid(
  int shmid,
  bool ignore_errors=false)
{
  if (shmctl(shmid, IPC_RMID, 0) == -1) {
    if (ignore_errors) return false;
    std::string msg = strerror(errno);
    throw std::runtime_error("shmctl(IPC_RMID) failure: " + msg);
  }
  return true;
}

void
wrap_all()
{
  using namespace boost::python;
  def("copy_to_new_segment", copy_to_new_segment<int>, (arg("data")));
  def("copy_to_new_segment", copy_to_new_segment<double>, (arg("data")));
  def("copy_from_segment_int", copy_from_segment<int>, (arg("shmid")));
  def("copy_from_segment_double", copy_from_segment<double>, (arg("shmid")));
  def("list_segment_ids", list_segment_ids, (
    arg("attached_status")=0,
    arg("ignore_errors")=false));
  def("shmctl_rmid", shmctl_rmid, (
    arg("shmid"),
    arg("ignore_errors")=false));
}

} // namespace ipctbx

BOOST_PYTHON_MODULE(ipctbx_ext)
{
  ipctbx::wrap_all();
}
