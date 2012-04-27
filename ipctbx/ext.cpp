#include <boost/python.hpp>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>

int
exercise_create_segment(
  std::size_t number_of_integers)
{
  int shmid = shmget(
    IPC_PRIVATE,
    number_of_integers * sizeof(int),
    IPC_CREAT|IPC_EXCL|S_IRUSR|S_IWUSR);
  if (shmid == -1){
    std::string msg = strerror(errno);
    throw std::runtime_error("shmget() failure: " + msg);
  }
  void* ptr = shmat(shmid, 0, 0);
  if (ptr == reinterpret_cast<void *>(-1)) {
    std::string msg = strerror(errno);
    throw std::runtime_error("shmat() failure: " + msg);
  }
  int* data = static_cast<int*>(ptr);
  for (std::size_t i=0;i<number_of_integers;i++) {
    data[i] = static_cast<int>(i % 13);
  }
  return shmid;
}

std::size_t
exercise_read_segment(
  int shmid,
  std::size_t number_of_integers)
{
  std::size_t result = 0;
  void* ptr = shmat(shmid, 0, SHM_RDONLY);
  if (ptr == reinterpret_cast<void *>(-1)) {
    std::string msg = strerror(errno);
    throw std::runtime_error("shmat() failure: " + msg);
  }
  int const* data = static_cast<int const*>(ptr);
  for (std::size_t i=0;i<number_of_integers;i++){
    if (data[i] == static_cast<int>(i % 13)) {
      result++;
    }
  }
  return result;
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
    struct shmid_ds shm_segment;
    int shmid = shmctl(i, SHM_STAT, &shm_segment);
    if (shmid < 0) {
      if (i == 0 && max_id == 0) break;
      if (ignore_errors) continue;
      std::string msg = strerror(errno);
      throw std::runtime_error("shmctl(SHM_STAT) failure: " + msg);
    }
    if (attached_status == 0
        || ((attached_status < 0) == (shm_segment.shm_nattch == 0))) {
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

BOOST_PYTHON_MODULE(ipctbx_ext)
{
  using namespace boost::python;
  def("exercise_create_segment", exercise_create_segment);
  def("exercise_read_segment", exercise_read_segment);
  def("list_segment_ids", list_segment_ids, (
    arg("attached_status")=0,
    arg("ignore_errors")=false));
  def("shmctl_rmid", shmctl_rmid, (
    arg("shmid"),
    arg("ignore_errors")=false));
}
