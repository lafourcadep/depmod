#include "lammps.h"
#include "linalg.h"

#include <iostream>

void lammps_generate_box_evolution_data_brute(
  Configuration* config,
  Deformation* deform,
  py::array_t<double> t,
  py::array_t<double> H,
  py::array_t<double> G,
  py::array_t<double> F
) {

  if (config->K == 1) {
    lmp_box_gen_brute_npts(
      config, 
      deform,
      t,
      H,
      G,
      F
    );
  } else {
    lmp_box_gen_brute_kpts(
      config, 
      deform,
      t,
      H,
      G,
      F
    );
  }
};

void lmp_box_gen_brute_npts(
  Configuration* config,
  Deformation* deform,
  py::array_t<double> t,
  py::array_t<double> H,
  py::array_t<double> G,
  py::array_t<double> F
) {

  size_t imax = config->N + 1;

  auto T = t.mutable_unchecked<1>();
  auto Ht = H.mutable_unchecked<2>();
  auto Gt = G.mutable_unchecked<2>();
  auto Ft = F.mutable_unchecked<2>();

  Map<const Matrix3d> H0(Ht.data(0, 0), 3, 3);

  Map<Matrix3d> Hi(NULL, 3, 3);
  Map<Matrix3d> Gi(NULL, 3, 3);
  Map<Matrix3d> Fi(NULL, 3, 3);
  Map<Matrix3d> Fj(NULL, 3, 3);

  Matrix3d K = deform->compute_K(config->dt, config->gammadot);
  
  size_t i, j; 

  for (i = 1, j = 0; i < imax; i++, j++) {

    T(i) = i * config->dt;

    // In principle memory is swapped with no call to allocator
    // TODO: Might be memory bound. Find a better solution?
    new (&Hi) Map<Matrix3d>(Ht.mutable_data(i, 0), 3, 3);
    new (&Gi) Map<Matrix3d>(Gt.mutable_data(i, 0), 3, 3);
    new (&Fi) Map<Matrix3d>(Ft.mutable_data(i, 0), 3, 3);
    new (&Fj) Map<Matrix3d>(Ft.mutable_data(j, 0), 3, 3);

    // F(i) = F(i-1) + dF
    Fi = Fj.array() + (K * Fj).array();

    // H(i) = F(i)*H(0)
    Hi = Fi * H0; 

    // Rotate back to lammps covention (i.e. triangular matrix)
    align_to_lammps_convention(Gi, Hi);
  };

};

void lmp_box_gen_brute_kpts(
  Configuration* config,
  Deformation* deform,
  py::array_t<double> t,
  py::array_t<double> H,
  py::array_t<double> G,
  py::array_t<double> F
) {

  size_t imax = config->N + 1;
  size_t kmax = config->K;

  auto T = t.mutable_unchecked<1>();
  auto Ht = H.mutable_unchecked<2>();
  auto Gt = G.mutable_unchecked<2>();
  auto Ft = F.mutable_unchecked<2>();

  Map<const Matrix3d> H0(Ht.data(0, 0), 3, 3);
  Map<const Matrix3d> F0(Ft.data(0, 0), 3, 3);

  Map<Matrix3d> Hi(NULL, 3, 3);
  Map<Matrix3d> Gi(NULL, 3, 3);
  Map<Matrix3d> Fi(NULL, 3, 3);

  Matrix3d Fi_tmp = Matrix3d::Zero();
  Matrix3d Fj_tmp = F0;

  Matrix3d K = deform->compute_K(config->dt, config->gammadot);

  size_t i, k;

  for (i = 1; i < imax; i++) {
    
    new (&Fi) Map<Matrix3d>(Ft.mutable_data(i, 0));
    new (&Hi) Map<Matrix3d>(Ht.mutable_data(i, 0));
    new (&Gi) Map<Matrix3d>(Gt.mutable_data(i, 0));
  
    for (k = 0; k < kmax; k++) {
      Fi_tmp = Fj_tmp.array() + (K * Fj_tmp).array();
      Fj_tmp = Fi_tmp;
    }

    T(i) = i * kmax * config->dt;

    Fi = Fi_tmp;
    Hi = Fi * H0;
    align_to_lammps_convention(Gi, Hi);
  };
};
