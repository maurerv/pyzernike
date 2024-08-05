#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ZernikeDescriptor.h"

namespace py = pybind11;

template <typename T> class PyZernikeDescriptor {
private:
  std::unique_ptr<ZernikeDescriptor<T, T>> descriptor;

public:
  PyZernikeDescriptor(
      py::array_t<T, py::array::c_style | py::array::forcecast> array,
      int order) {
    py::buffer_info buf = array.request();
    if (buf.ndim != 3)
      throw std::runtime_error("Number of dimensions must be 3");
    int dim = buf.shape[0];
    if (buf.shape[1] != dim || buf.shape[2] != dim)
      throw std::runtime_error("Input must be a cube");

    T *ptr = static_cast<T *>(buf.ptr);
    descriptor = std::make_unique<ZernikeDescriptor<T, T>>(ptr, dim, order);
  }

  py::array_t<T> GetDescriptors() const {
    const std::vector<T> &invariants = descriptor->GetInvariants();
    auto result = py::array_t<T>(invariants.size());
    py::buffer_info buf = result.request();
    std::copy(invariants.begin(), invariants.end(), static_cast<T *>(buf.ptr));
    return result;
  }

  py::array_t<std::complex<T>> Reconstruct(int dim, int minN = 0,
                                           int maxN = 100, int minL = 0,
                                           int maxL = 100) {
    using ComplexT = std::complex<T>;
    std::vector<std::vector<std::vector<ComplexT>>> grid(
        dim,
        std::vector<std::vector<ComplexT>>(dim, std::vector<ComplexT>(dim)));
    descriptor->Reconstruct(grid, minN, maxN, minL, maxL);

    auto result = py::array_t<ComplexT>({dim, dim, dim});
    py::buffer_info buf = result.request();
    ComplexT *ptr = static_cast<ComplexT *>(buf.ptr);

    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        for (int k = 0; k < dim; ++k) {
          ptr[(i * dim + j) * dim + k] = grid[i][j][k];
        }
      }
    }

    return result;
  }

  void SaveInvariants(const char *fname) { descriptor->SaveInvariants(fname); }
};

template <typename T>
void declare_zernike_descriptor(py::module &m, const std::string &type_suffix) {
  using Class = PyZernikeDescriptor<T>;
  std::string pyclass_name = std::string("ZernikeDescriptor") + type_suffix;

  py::class_<Class>(m, pyclass_name.c_str())
      .def(py::init<py::array_t<T, py::array::c_style | py::array::forcecast>,
                    int>(),
           py::arg("data"), py::arg("order"))
      .def("get_descriptors", &Class::GetDescriptors)
      .def("save_invariants", &Class::SaveInvariants)
      .def("reconstruct", &Class::Reconstruct, py::arg("dim"),
           py::arg("min_n") = 0, py::arg("max_n") = 100, py::arg("min_l") = 0,
           py::arg("max_l") = 100);
}

PYBIND11_MODULE(_core, m) {
  declare_zernike_descriptor<float>(m, "Float");
  declare_zernike_descriptor<double>(m, "Double");
};
