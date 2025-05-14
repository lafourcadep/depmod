#include "lammps.h"
#include "core.h"
#include "linalg.h"
#include <iostream>


static const std::map<uint32_t, lmp_function_t> LAMMPS_METHODS = {
  { INTEGRATION_METHOD::BRUTE,                            &_lammps_box_gen_brute_n },
  { INTEGRATION_METHOD::BRUTE | INTEGRATION_METHOD::KPTS, &_lammps_box_gen_brute_k },
};

const std::map<uint32_t, lmp_function_t>& lammps_methods() {
  return LAMMPS_METHODS;
};

void _lammps_box_gen_brute_n(Configuration* config, Deformation* deform, py::array_t<double> t, py::array_t<double> H, py::array_t<double> G, py::array_t<double> F) {

  size_t i_max = config->N + 1;

  // Mutable memory of py::array_t with no bounds check
  auto T = t.mutable_unchecked<1>();
  auto Ht = H.mutable_unchecked<2>(); 
  auto Gt = G.mutable_unchecked<2>();
  auto Ft = F.mutable_unchecked<2>();

  // Eigen::Map pointing to py::array_t data
  Map<const Matrix3d> H0(Ht.data(0, 0), 3, 3);
  Map<const Matrix3d> F0(Ft.data(0, 0), 3, 3);

  Map<Matrix3d> Hi(NULL, 3, 3);
  Map<Matrix3d> Gi(NULL, 3, 3);
  Map<Matrix3d> Fi(NULL, 3, 3); // F at time t
  // Map<Matrix3d> Fj(NULL, 3, 3); // F at time t - dt
  
  Matrix3d Fj = F0;

  // Deformation factor
  Matrix3d K = deform->compute_K(config->dt, config->gammadot);

  size_t i, j;

  for (i = 1, j = 0; i < i_max; i++, j++) {
    T(i) = i * config->dt;

    // Memory map
    // In theory memory is swapped without call to allocator
    new (&Hi) Map<Matrix3d>(Ht.mutable_data(i, 0), 3, 3);
    new (&Gi) Map<Matrix3d>(Gt.mutable_data(i, 0), 3, 3);
    new (&Fi) Map<Matrix3d>(Ft.mutable_data(i, 0), 3, 3);
    // new (&Fj) Map<Matrix3d>(Ft.mutable_data(j, 0), 3, 3);

    // F(i) = F(i - 1) + dF; dF = (K * F(i - 1))
    Fi = (Fj.array() + (K * Fj).array()).matrix();
    Fj = Fi;

    // H(i) = F(i) * H(0)
    Hi = Fi * H0;

    // Align to lammps convention (i.e. upper triangular matrix)
    align_to_lammps_convention(Gi, Hi);
  }
};

void _lammps_box_gen_brute_k(Configuration* config, Deformation* deform, py::array_t<double> t, py::array_t<double> H, py::array_t<double> G, py::array_t<double> F) {

  size_t i_max = config->N + 1;
  size_t k_max = config->K;

  // Mutable memory of py::array_t with no bounds check
  auto T = t.mutable_unchecked<1>();
  auto Ht = H.mutable_unchecked<2>(); 
  auto Gt = G.mutable_unchecked<2>();
  auto Ft = F.mutable_unchecked<2>();

  // Eigen::Map pointing to py::array_t data
  Map<const Matrix3d> H0(Ht.data(0, 0), 3, 3);
  Map<const Matrix3d> F0(Ft.data(0, 0), 3, 3);

  Map<Matrix3d> Hi(NULL, 3, 3);
  Map<Matrix3d> Gi(NULL, 3, 3);
  Map<Matrix3d> Fi(NULL, 3, 3); // F at time t
  
  // Temp matrix
  Matrix3d tempFi = Matrix3d::Zero();
  Matrix3d tempFj = F0;

  // Deformation factor
  Matrix3d K = deform->compute_K(config->dt, config->gammadot);

  size_t i, k;

  for (i = 1; i < i_max; i++) {

    // Memory map
    new (&Hi) Map<Matrix3d>(Ht.mutable_data(i, 0), 3, 3);
    new (&Gi) Map<Matrix3d>(Gt.mutable_data(i, 0), 3, 3);
    new (&Fi) Map<Matrix3d>(Ft.mutable_data(i, 0), 3, 3);

    for (k = 0; k < k_max; k++) {
      tempFi = (tempFj.array() + (K * tempFj).array()).matrix();
      tempFj = tempFi;
    }
    
    T(i) = k_max * i * config->dt;

    Fi = tempFi;
    Hi = Fi * H0;

    align_to_lammps_convention(Gi, Hi);
  }
};
