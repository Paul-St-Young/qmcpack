#include "ScreenedDefect.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"
#include "LongRange/Screen2DKernel.h"

#include <cmath>
#include <vector>

namespace qmcplusplus
{

ScreenedDefect::Return_t ScreenedDefect::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  const auto& d_ab(P.getDistTableAB(itab));
  for (int i=0; i<nelec; i++)
  {
    const auto& dists = d_ab.getDistRow(i);
    for (int j=0; j<ndefect; j++)
    {
      const RealType r = dists[j];
      if (r >= rmax_spline) continue; // tail of screened potential is zero
      value_ += rVspline->splint(r) / r;
    }
  }
  value_ = vconst + target_charge*charge*value_;
  return value_; // !!!! this does NOT go to scalar.dat
}

bool ScreenedDefect::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(dgate, "dgate");
  attrib.add(mimg, "max_image");
  attrib.put(cur);
  vconst = 0.5*nelec * (-1.0*2*M_PI*dgate/area);

  // Tabulate r * V_screened(r) on a linear grid in [0, rmax_spline].
  // Spline to accelerate evaluate.
  rmax_spline      = 15.0 * dgate;
  const int ngrid  = 1024;
  myGrid           = std::make_shared<LinearGrid<RealType>>();
  myGrid->set(0.0, rmax_spline, ngrid);

  std::vector<RealType> rv(ngrid);
  for (int ig = 1; ig < ngrid - 1; ++ig)
  {
    const RealType r = (*myGrid)[ig];
    rv[ig]           = r * screen2DKernel<RealType>(r, dgate, mimg);
  }
  // enforce limits
  rv[0]         = 1.0;
  rv[ngrid - 1] = 0.0;
  const RealType deriv0  = (rv[1] - rv[0]) / ((*myGrid)[1] - (*myGrid)[0]);
  rVspline               = std::make_shared<OneDimCubicSpline<RealType>>(myGrid->makeClone(), rv);
  rVspline->spline(0, deriv0, ngrid - 1, 0.0);
  return true;
}

bool ScreenedDefect::get(std::ostream& os) const
{
  os << "Screened defect potential" << std::endl;
  os << "  ndefect = " <<  ndefect << std::endl;
  os << "  charge  = " <<  charge << std::endl;
  os << "  nelec   = " <<  nelec << std::endl;
  os << "  dgate   = " <<  dgate << " bohr*" << std::endl;
  os << "  mimg    = " <<  mimg << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> ScreenedDefect::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<ScreenedDefect>(*this);
}

} // qmcplusplus
