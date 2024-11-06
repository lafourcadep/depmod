#include "core.h"
#include "lammps.h"
#include "exastamp.h"


void lammps_generate_box_evolution_data(INTEGRATION_METHOD method, Configuration* config, Deformation* deform, py::array_t<double> t, py::array_t<double> H, py::array_t<double> G, py::array_t<double> F) {
  if (config->K > 1)
    method = (INTEGRATION_METHOD) (method | INTEGRATION_METHOD::KPTS);

  if (! lammps_methods().count(method))
    throw std::runtime_error("Invalid integration method bitflag");
  lammps_methods().at(method)(config, deform, t, H, G, F);
}

void exastamp_generate_box_evolution_data(INTEGRATION_METHOD method, Configuration* config, Deformation* deform, py::array_t<double> t, py::array_t<double> F) {
  if (config->K > 1)
    method = (INTEGRATION_METHOD) (method | INTEGRATION_METHOD::KPTS);
  
  if (! exastamp_methods().count(method))
    throw std::runtime_error("Invalid integration method bitflag");
  exastamp_methods().at(method)(config, deform, t, F);
}
