//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, Hofstra University
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, Hofstra University
//////////////////////////////////////////////////////////////////////////////////////

#include "CosinePotential.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

CosinePotential::Return_t CosinePotential::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  const size_t Nelec = P.getTotalNum();
  for (size_t iel = 0; iel < Nelec; iel++)
  {
    const auto& r = (P.getActivePtcl() == iel) ? P.activeR(iel) : P.R[iel];
    value_ += 2*vq*std::cos(dot(qvec, r));
  }
  return value_;
}

bool CosinePotential::put(xmlNodePtr cur)
{
  // read inputs
  OhmmsAttributeSet attrib;
  attrib.add(vq, "vq");
  attrib.add(qvec, "qvec");
  attrib.put(cur);
  auto sc = myP.getSimulationCell();
  auto cell = sc.getLattice();
  auto gidx = cell.k_unit(qvec);
  for (size_t l=0;l<gidx.size();l++)
  {
    if (std::abs(std::round(gidx[l])-gidx[l]) > 1e-6)
    {
      std::ostringstream msg;
      msg << " invalid qvec input " << gidx << " - not integer" << std::endl;
      throw std::runtime_error(msg.str());
    }
  }
  return true;
}

bool CosinePotential::get(std::ostream& os) const
{
  os << "External cosine potential" << std::endl;
  os << "  vq = " <<  vq << " Ha" << std::endl;
  os << "  qvec = " <<  qvec << " 1/Bohr" << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> CosinePotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<CosinePotential>(*this);
}

} // qmcplusplus
