#ifndef QMCPLUSPLUS_SCREENED_DEFECT_H
#define QMCPLUSPLUS_SCREENED_DEFECT_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class ScreenedDefect : public OperatorBase
{
public:
  ScreenedDefect(ParticleSet& P, ParticleSet& sP)
  : itab(P.addTable(sP)), ndefect(sP.getTotalNum()), nelec(P.getTotalNum()), mimg(100), dgate(1.0), charge(1.0), vconst(0.0) {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
    area = P.getLattice().Volume / P.getLattice().R(2, 2);
    const SpeciesSet& target_species = P.getSpeciesSet();
    const SpeciesSet& source_species = sP.getSpeciesSet();
    const int source_charge_idx      = source_species.findAttribute("charge");
    const int target_charge_idx      = target_species.findAttribute("charge");
    target_charge = target_species(target_charge_idx, 0);
    charge = source_species(source_charge_idx, 0);
  };
  ~ScreenedDefect(){};
  Return_t evaluate(ParticleSet& P) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) override;
  // ---- begin required overrides
  void resetTargetParticleSet(ParticleSet& P) override {APP_ABORT("not implemented");};
  std::string getClassName() const override {return "defect";};
  // required overrides end ----
private:
  const int itab;
  int ndefect, nelec, mimg;
  RealType dgate, charge, vconst, area;
  RealType target_charge;
};
} // qmcplusplus
#endif
