#include "config.h"

Configuration::Configuration(
  double gammadot,
  double t_min,
  double t_max,
  size_t n,
  size_t k
):
  gammadot(gammadot),
  t_min(t_min),
  t_max(t_max),
  N(n),
  K(k)
{
  M = N * K;
  dt = (t_max - t_min) / (double) M;
};
