#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "core.h"

namespace py = pybind11;

void init_core(py::module_ &m) {
  
  auto subm = m.def_submodule("core");

  subm.doc() = "depmod core sub-module";

  py::enum_<INTEGRATION_METHOD>(subm, "_lib_IntegrationMethod")
    .value("BRUTE", INTEGRATION_METHOD::BRUTE)
    .value("CENTER", INTEGRATION_METHOD::CENTER)
    .value("BRUZY", INTEGRATION_METHOD::BRUZY)
    .value("KPTS", INTEGRATION_METHOD::KPTS);

  // Main functions to be called from python.
  subm.def("_lib_lammps_generate_box_evolution_data", &lammps_generate_box_evolution_data);
  subm.def("_lib_exastamp_generate_box_evolution_data", &exastamp_generate_box_evolution_data);
}
