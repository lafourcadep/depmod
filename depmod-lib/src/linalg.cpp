#include "linalg.h"


Matrix3d rotation_matrix_from_vectors(Ref<Vector3d> u, Ref<Vector3d> v) {

  Vector3d a = unitvec(u);
  Vector3d b = unitvec(v);

  Vector3d x = a.cross(b);
  
  double c = a.dot(b);
  double s = x.norm(); 

  Matrix3d K {{0., -x(2), x(1)}, {x(2), 0., -x(0)}, {-x(1), x(0), 0.}};
  Matrix3d Q = Matrix3d::Identity().array() + K.array() + ((K * K).array() * ((1 - c) / (s * s)));

  return Q;
};


Matrix3d rotation_matrix_around_axis(Ref<Vector3d> axis, double angle) {

  double angle_rad = DEG_TO_RAD * angle;

  Vector3d a = unitvec(axis);

  Matrix3d K {{0., a(2), -a(1)}, {-a(2), 0., a(0)}, {a(1), -a(0), 0.}};
  Matrix3d Q = Matrix3d::Identity().array() - K.array() * sin(angle_rad) + ((K * K).array() * (1 - cos(angle_rad)));

  return Q;
};

void align_to_lammps_convention(Ref<Matrix3d> G, Ref<Matrix3d> H) {
  Vector3d zaxis {0., 0., 1.};

  align_to_lammps_convention(G, H, zaxis);
}

void align_to_lammps_convention(Ref<Matrix3d> G, Ref<Matrix3d> H, Ref<Vector3d> zaxis) {

  Vector3d a = H.col(0);
  Vector3d b = H.col(1);

  Vector3d axb = a.cross(b).normalized();
  
  Matrix3d temp = Matrix3d::Zero();
  Matrix3d Q = Matrix3d::Identity();

  double res = (axb.array() - zaxis.array()).matrix().norm();

  // 1: Bring back a and b in (x, y) plane by aligning a^b with z-axis
  if (res < 1.0e-5) {
    temp = H;
  } else {
    Q = rotation_matrix_from_vectors(axb, zaxis);
    temp = Q * H;
  };

  // 2: bring back a along x-axis by rotation around z-axis
  double xdir =  RAD_TO_DEG * atan(H(1, 0) / H(0, 0));

  Q = rotation_matrix_around_axis(zaxis, -1. * xdir);
  G = Q * temp;
};
