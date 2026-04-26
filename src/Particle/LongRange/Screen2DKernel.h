#ifndef QMCPLUSPLUS_SCREEN2D_KERNEL_H
#define QMCPLUSPLUS_SCREEN2D_KERNEL_H

#include <cmath>

namespace qmcplusplus
{
/** Screened defect potential between two grounded plates separated by 2*dgate.
 *
 *  V(r) = sum_{m=-inf}^{inf} (-1)^m / sqrt(r^2 + (2*m*dgate)^2)
 *       = (2/dgate) * sum_{k=0}^{inf} K_0((2k+1)*pi*r/(2*dgate))
 *
 *  References:
 *    Throckmorton & Vafek, PRB 86, 115447 (2012)
 *    Valenti et al., arXiv:2307.15119, eq. (B4)
 *
 *  Implementation:
 *    r < dgate/4 : real-space image sum truncated at 10*mimg
 *    otherwise   : Bessel K_0 series truncated at mimg (rapidly convergent)
 */
template<typename T>
inline T screen2DKernel(T r, T dgate, int mimg)
{
  T v = 0;
  if (r < dgate / 4)
  {
    const int msr = 10 * mimg;
    for (int m = -msr; m <= msr; ++m)
    {
      const T d = T(2) * dgate * m;
      v += std::pow(T(-1), m) / std::sqrt(r * r + d * d);
    }
    return v;
  }
  const T arg = r / (T(2) * dgate);
  for (int k = 0; k <= mimg; ++k)
    v += std::cyl_bessel_k(0, (T(2) * k + T(1)) * M_PI * arg);
  return T(2) * v / dgate;
}
} // namespace qmcplusplus
#endif
