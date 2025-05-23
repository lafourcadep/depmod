#pragma once
#include <memory>

#include "eigen.h"
#include "exprtk.hpp"

struct StrainRate {
public:
  StrainRate() = default;
  StrainRate(double strain_rate) : m_strain_rate(strain_rate) {};
  virtual ~StrainRate() = default;
  virtual inline double operator()(double) const { return m_strain_rate; };

protected:
  double m_strain_rate;
};

class CustomExpressionStrainRate : public StrainRate {
public:
  CustomExpressionStrainRate(const std::string& expr_str) {
    m_symbol_table.add_variable("t", m_t);
    m_symbol_table.add_variable("T", m_t);
    m_symbol_table.add_constants();
    m_expression.register_symbol_table(m_symbol_table);

    exprtk::parser<double> parser;
    if (!parser.compile(expr_str, m_expression)) {
      throw std::runtime_error("Invalid expression: " + expr_str);
    }
  }

  inline double operator()(double t) const override {
    const_cast<double&>(m_t) = t;
    return m_expression.value();
  }

private:
  exprtk::expression<double> m_expression;
  exprtk::symbol_table<double> m_symbol_table;
  double m_t;
};

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
  DeformationPath(Deformation deformation, std::shared_ptr<StrainRate>& strain_rate, double tmin,
                  double tmax, size_t npts, size_t kpts)
      : m_deformation(deformation), m_tmin(tmin), m_tmax(tmax), m_npts(npts), m_kpts(kpts),
        m_strain_rate(strain_rate) {

    m_dt = (m_tmax - m_tmin) / static_cast<double>(kpts * npts);
  }

  inline const std::shared_ptr<StrainRate>& strain_rate() const { return m_strain_rate; }
  inline double tmin() const { return m_tmin; }
  inline double tmax() const { return m_tmax; }
  inline double imax() const { return npts() + 1; }
  inline double dt() const { return m_dt; }

  inline size_t npts() const { return m_npts; }
  inline size_t kpts() const { return m_kpts; }

  inline Matrix3d Lt(double t) {
    return (m_deformation.cfac() * (*m_strain_rate)(t)*m_deformation.S().array()).matrix();
  }

private:
  Deformation m_deformation;
  double m_tmin, m_tmax, m_dt;
  size_t m_npts, m_kpts;
  std::shared_ptr<StrainRate> m_strain_rate = nullptr;
};
