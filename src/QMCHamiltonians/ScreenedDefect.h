#ifndef QMCPLUSPLUS_SCREENED_DEFECT_H
#define QMCPLUSPLUS_SCREENED_DEFECT_H

#include <memory>
#include <vector>

#include "QMCHamiltonians/OperatorBase.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimGridBase.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{
class ScreenedDefect : public OperatorBase
{
public:
  ScreenedDefect(ParticleSet& P, ParticleSet& sP)
  : itab(P.addTable(sP)), ndefect(sP.getTotalNum()), nelec(P.getTotalNum()), mimg(100), mlat(0), dgate(1.0), charge(1.0), vconst(0.0), rmax_spline(0.0) {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
    const auto& cell = P.getLattice();
    area   = P.getLattice().Volume / cell.R(2, 2);
    lattR_ = cell.R;
    for (size_t l = 0; l < OHMMS_DIM; ++l)
      bcs_[l] = cell.BoxBConds[l];
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
  int ndefect, nelec, mimg, mlat;
  RealType dgate, charge, vconst, area;
  RealType target_charge;
  RealType rmax_spline;
  std::shared_ptr<LinearGrid<RealType>> myGrid;
  std::shared_ptr<OneDimCubicSpline<RealType>> rVspline;
  Tensor<RealType, OHMMS_DIM> lattR_;
  TinyVector<int, OHMMS_DIM> bcs_;
  std::vector<TinyVector<RealType, OHMMS_DIM>> lats_;
};
} // qmcplusplus
#endif
