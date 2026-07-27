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

#ifndef QMCPLUSPLUS_COSINE_POTENTIAL_H
#define QMCPLUSPLUS_COSINE_POTENTIAL_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class CosinePotential : public OperatorBase
{
public:
  CosinePotential(ParticleSet& P) : myP(P)
  {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
  };
  ~CosinePotential(){};
  Return_t evaluate(ParticleSet& P) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) override;
  // ---- begin required overrides
  void resetTargetParticleSet(ParticleSet& P) override {APP_ABORT("not implemented");};
  std::string getClassName() const override {return "cosine";};
  // required overrides end ----
private:
  RealType vq;
  TinyVector<RealType, 3> qvec;
  bool spinsus = false;
  ParticleSet& myP;
};
} // qmcplusplus
#endif
