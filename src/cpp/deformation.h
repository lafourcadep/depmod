#pragma once

#include "eigen.h"

// Class that represent a deformation
class Deformation {
public:
  Deformation(const Matrix3d& S) : m_S(S), m_cfac(1.0) {}
  Deformation(const Matrix3d& S, double cfac) : m_S(S), m_cfac(cfac) {}
  inline const Matrix3d S(void) { return m_S; }
  inline double cfac(void) { return m_cfac; }

protected:
  Matrix3d m_S;
  double m_cfac;
};

// Class that represent a deformation path (i.e. a deformation + strain_rate)
class DeformationPath {
public:
  DeformationPath(Deformation deformation, double strain_rate, double tmin, double tmax,
                  size_t npts, size_t kpts)
      : m_deformation(deformation), m_strain_rate(strain_rate), m_tmin(tmin), m_tmax(tmax),
        m_npts(npts), m_kpts(kpts) {

    m_dt = (m_tmax - m_tmin) / static_cast<double>(kpts * npts);
  }

  inline double strain_rate() const { return m_strain_rate; }
  inline double tmin() const { return m_tmin; }
  inline double tmax() const { return m_tmax; }
  inline double imax() const { return npts() + 1; }
  inline double dt() const { return m_dt; }

  inline size_t npts() const { return m_npts; }
  inline size_t kpts() const { return m_kpts; }

  inline Matrix3d L() {
    return (m_deformation.cfac() * m_strain_rate * m_deformation.S().array()).matrix();
  }

private:
  Deformation m_deformation;
  double m_strain_rate, m_tmin, m_tmax, m_dt;
  size_t m_npts, m_kpts;
};
