#include "ScreenedDefect.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"

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
      auto r = dists[j];
      auto r2 = r*r;
      for (int m=-mimg;m<=mimg;m++)
      {
        value_ += std::pow(-1, m)/std::sqrt(
          r2 + (2*dgate*m)*(2*dgate*m)
        );
      }
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
