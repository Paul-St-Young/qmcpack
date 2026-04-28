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
    const auto& disps = d_ab.getDisplRow(i);
    for (int j=0; j<ndefect; j++)
    {
      for (const auto& dlat : lats_)
      {
        const auto dr  = disps[j] + dlat;
        const RealType r = std::sqrt(dot(dr, dr));
        if (r >= rmax_spline) continue; // screened potential tail -> 0
        value_ += rVspline->splint(r) / r;
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
  const RealType vsr_k0 = 2 * M_PI * dgate;
  vconst = -(nelec*target_charge) * (ndefect*charge) * vsr_k0 / area;

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

  // build lattice translation list: include images that can come within rmax_spline
  // of any in-cell pair displacement (similar pattern to CoulombPBCAA::evalSR).
  // Only sum images along periodic directions (BoxBConds[i] == 1).
  lats_.clear();
  RealType max_a = 0;
  RealType min_a = std::numeric_limits<RealType>::max();
  bool any_periodic = false;
  for (size_t i = 0; i < OHMMS_DIM; ++i)
  {
    if (!bcs_[i]) continue;
    any_periodic = true;
    RealType ai_sq = 0;
    for (size_t l = 0; l < OHMMS_DIM; ++l)
      ai_sq += lattR_(i, l) * lattR_(i, l);
    const RealType ai = std::sqrt(ai_sq);
    max_a             = std::max(max_a, ai);
    min_a             = std::min(min_a, ai);
  }
  if (!any_periodic)
  {
    mlat = 0;
    lats_.push_back(TinyVector<RealType, OHMMS_DIM>(0));
    return true;
  }
  // safe upper bound on MIC pair displacement; mlat captures any image whose
  // translation can bring an out-of-rmax_spline pair within range.
  const RealType disp_bound = max_a;
  mlat                      = static_cast<int>(std::ceil((rmax_spline + disp_bound) / min_a));
  int nlat[OHMMS_DIM];
  for (size_t i = 0; i < OHMMS_DIM; ++i)
    nlat[i] = bcs_[i] ? mlat : 0;
  for (int ix = -nlat[0]; ix <= nlat[0]; ++ix)
    for (int iy = -nlat[1]; iy <= nlat[1]; ++iy)
      for (int iz = -nlat[2]; iz <= nlat[2]; ++iz)
      {
        TinyVector<RealType, OHMMS_DIM> R = 0;
        for (size_t l = 0; l < OHMMS_DIM; ++l)
          R[l] = ix * lattR_(0, l) + iy * lattR_(1, l) + iz * lattR_(2, l);
        // prune images that are too far for any in-cell pair to reach
        if (std::sqrt(dot(R, R)) - disp_bound > rmax_spline)
          continue;
        lats_.push_back(R);
      }
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
  os << "  mlat    = " <<  mlat << " (" << lats_.size() << " lattice images)" << std::endl;
  os << "  rmax    = " <<  rmax_spline << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> ScreenedDefect::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<ScreenedDefect>(*this);
}

} // qmcplusplus
