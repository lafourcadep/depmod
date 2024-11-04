#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>

namespace py = pybind11;

void init_config(py::module_&);
void init_lammps(py::module_&);
void init_deformation(py::module_&);

PYBIND11_MODULE(_lib, m) {
  m.doc() = "DEPMOD C++ lib";
  
  init_config(m);
  init_deformation(m);
  init_lammps(m);

};
