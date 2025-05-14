#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <Eigen/Dense>

#include "deformation.h"
#include "eigen.h"
#include "linalg.h"
#include "parsers.h"

namespace py = pybind11;

template <bool LAMMPS>
void generate_box_evolution_data(DeformationPath* defpath, py::array_t<double> buffer) {

  constexpr size_t index_H = 2;
  constexpr size_t index_G = 11;
  constexpr size_t index_F = 20;

  size_t imax = defpath->imax();
  size_t kmax = defpath->kpts();
  double dt = defpath->dt();

  Matrix3d L = defpath->L(); // velocity_gradient

  auto buf = buffer.mutable_unchecked<2>();

  Map<const Matrix3d> H0(buf.data(0, index_H), 3, 3);
  Map<const Matrix3d> F0(buf.data(0, index_F), 3, 3);

  Map<Matrix3d> Hi(NULL, 3, 3);
  Map<Matrix3d> Gi(NULL, 3, 3);
  Map<Matrix3d> Fi(NULL, 3, 3); // F at time t

  Matrix3d Fj = F0;

  for (size_t i = 1; i < imax; ++i) {

    new (&Hi) Map<Matrix3d>(buf.mutable_data(i, index_H), 3, 3);
    new (&Gi) Map<Matrix3d>(buf.mutable_data(i, index_G), 3, 3);
    new (&Fi) Map<Matrix3d>(buf.mutable_data(i, index_F), 3, 3);

    for (size_t k = 0; k < kmax; ++k) {
      // Runge-Kutta 4 integration scheme to solve ODE
      Matrix3d k1 = L * Fj;
      Matrix3d k2 = L * (Fj.array() + 0.5 * dt * k1.array()).matrix();
      Matrix3d k3 = L * (Fj.array() + 0.5 * dt * k2.array()).matrix();
      Matrix3d k4 = L * (Fj.array() + 1.0 * dt * k3.array()).matrix();
      Fi =
          (Fj.array() + (dt / 6.0) * (k1.array() + 2.0 * k2.array() + 2. * k3.array() + k4.array()))
              .matrix();
      Fj = Fi;
    }

    Hi = Fi * H0;

    if constexpr (LAMMPS) {
      align_to_lammps_convention(Gi, Hi);
    } else {
      Gi = Hi;
    }

    buf(i, 0) = i;
    buf(i, 1) = i * kmax * dt;
  }
};

void read_lattice_from_file(py::array_t<double>& x, const std::string& file,
                            const std::string& format, const std::string& compression) {
  auto cell = x.mutable_unchecked<2>();
  io::Context ctx = io::read_atom_file(file, format, compression);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      cell(i, j) = ctx.cell[i][j];
    }
  }
}

PYBIND11_MODULE(_lib, m) {
  m.doc() = "depmod c++ module";

  py::class_<Deformation>(m, "_Deformation")
      .def(py::init<const Matrix3d&>())
      .def(py::init<const Matrix3d&, double>())
      .def("S", &Deformation::S);

  py::class_<DeformationPath>(m, "_DeformationPath")
      .def(py::init<Deformation, double, double, double, size_t, size_t>())
      .def("strain_rate", &DeformationPath::strain_rate)
      .def("tmin", &DeformationPath::tmin)
      .def("tmax", &DeformationPath::tmax)
      .def("npts", &DeformationPath::npts)
      .def("kpts", &DeformationPath::kpts)
      .def("velocity_gradient", &DeformationPath::L);

  m.def("read_lattice_from_file", &read_lattice_from_file);
  m.def("lammps_generate_box_data", &generate_box_evolution_data<true>);
  m.def("exastamp_generate_box_data", &generate_box_evolution_data<false>);
}
