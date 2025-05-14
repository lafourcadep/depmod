#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "config.h"

namespace py = pybind11;

void init_config(py::module_ &m) {
  
  auto subm = m.def_submodule("config");

  subm.doc() = "depmod config module";
  
  py::class_<Configuration>(subm, "_lib_Configuration")
    .def(py::init<double, double, double, size_t, size_t>())
    .def_readonly("N", &Configuration::N);
}
