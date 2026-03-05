#include "ScreenedDefect.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

ScreenedDefect::Return_t ScreenedDefect::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  return value_;
}

bool ScreenedDefect::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(dgate, "dgate");
  attrib.put(cur);
  return true;
}

bool ScreenedDefect::get(std::ostream& os) const
{
  os << "Screened defect potential" << std::endl;
  os << "  dgate = " <<  dgate << " bohr*"<< std::endl;
  return true;
}

std::unique_ptr<OperatorBase> ScreenedDefect::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<ScreenedDefect>(*this);
}

} // qmcplusplus
