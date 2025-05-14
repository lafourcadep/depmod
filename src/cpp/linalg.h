#pragma once

#include "dtype.h"

static constexpr const double RAD_TO_DEG = 180. / M_PI;
static constexpr const double DEG_TO_RAD = M_PI / 180.;

inline void normalize(Ref<Vector3d> v) {
  v.normalize();
};

inline Vector3d unitvec(Ref<Vector3d> v) {
  return v.normalized();
};

Matrix3d rotation_matrix_around_axis(Ref<Vector3d> axis, double angle);

Matrix3d rotation_matrix_from_vectors(Ref<Vector3d> u, Ref<Vector3d> v);

void align_to_lammps_convention(Ref<Matrix3d> G, Ref<Matrix3d> H);

void align_to_lammps_convention(Ref<Matrix3d> G, Ref<Matrix3d> H, Ref<Vector3d> zaxis);
