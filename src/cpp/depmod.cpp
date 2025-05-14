#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>

namespace py = pybind11;

// Forward declaration of the submodules initialization functions.
void init_config(py::module_&);
void init_deformation(py::module_&);
void init_core(py::module_&);

// actual definition of the python module.
PYBIND11_MODULE(_lib, m) {
  m.doc() = "depmod c++ module";
  
  // init the different submodules.
  init_config(m);
  init_deformation(m);
  init_core(m);
};
