#ifndef QMCPLUSPLUS_SCREENED_DEFECT_H
#define QMCPLUSPLUS_SCREENED_DEFECT_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class ScreenedDefect : public OperatorBase
{
public:
  ScreenedDefect(ParticleSet& P, ParticleSet& sP)
  : source(sP) {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
  };
  ~ScreenedDefect(){};
  Return_t evaluate(ParticleSet& P) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) override;
  // ---- begin required overrides
  void resetTargetParticleSet(ParticleSet& P) override {APP_ABORT("not implemented");};
  std::string getClassName() const override {return "moire";};
  // required overrides end ----
private:
  RealType dgate;
  ParticleSet& source;
};
} // qmcplusplus
#endif
