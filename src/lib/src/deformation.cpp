#include "deformation.h"

/* Deformation */
Deformation::Deformation(const Matrix3d& S): m_S(S) {};


/* Traction - Compression */

TractionCompression::TractionCompression(const Matrix3d& S, bool cflag): Deformation(S), m_cflag(cflag) {};

TractionCompression::TractionCompression(const Matrix3d& S): TractionCompression(S, false) {};

Matrix3d TractionCompression::compute_K(double dt, double gammadot) {
  double c = (m_cflag) ? -1. : 1.; // swith for compression / traction
  return (c * dt * gammadot * m_S.array()).matrix();
};


/* Traction - Compression Isochore */

TractionCompressionIsoV::TractionCompressionIsoV(const Matrix3d& S, bool cflag): Deformation(S), m_cflag(cflag) {};

TractionCompressionIsoV::TractionCompressionIsoV(const Matrix3d& S): TractionCompressionIsoV(S, false) {};

Matrix3d TractionCompressionIsoV::compute_K(double dt, double gammadot) {
  
  double c = (m_cflag) ? -1. : 1.; // swith for compression / traction
  
  double gamma = dt * c * gammadot;
  double alpha = sqrt(1 / (1 + gamma)) - 1;
  
  Matrix3d F_gamma = gamma * m_S.array();
  Matrix3d F_alpha = alpha * (Matrix3d::Identity().array() - m_S.array());

  return (F_gamma + F_alpha).matrix();
};


/* Pure Shear */

PureShear::PureShear(const Matrix3d& S): Deformation(S) {};

Matrix3d PureShear::compute_K(double dt, double gammadot) {
  return (dt * gammadot * m_S.array()).matrix();
};
