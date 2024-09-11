#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

#include "lammps.h"

void init_lammps(py::module_ &m) {
  
  auto subm = m.def_submodule("lammps");

  subm.doc() = "Lammps sub-module";

  subm.def("_lib_lammps_generate_box_evolution_data_brute", &lammps_generate_box_evolution_data_brute);
}
