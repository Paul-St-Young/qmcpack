//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/LongRange/Screen2D.h"
#include "Particle/LongRange/Screen2DKernel.h"

namespace qmcplusplus
{

Screen2D::Screen2D(ParticleSet& ref, mRealType kc_in)
  : LRHandlerBase(kc_in), dgate(ref.getLattice().dgate), mimg(ref.getLattice().mimg)
{
  if (ref.getLattice().ndim != 2)
    throw std::runtime_error("2D potential requires 2D Lattice");
  if (dgate < 0)
    throw std::runtime_error("screened potential requires distance_to_gate");
  LR_rc = ref.getLattice().LR_rc; // CoulombPBC needs get_rc() to createSpline4RbyVs
  LR_kc = ref.getLattice().LR_kc; // get_kc() is used in QMCFiniteSize
  area = ref.getLattice().Volume/ref.getLattice().R(2,2);
  int nlat = ref.getLattice().nlat;
  mRealType reach = LR_rc + nlat*2.0*LR_rc;
  // find necessary rcut
  mRealType new_rcut = findRcut();
  LR_rc = std::max(new_rcut, LR_rc);
  // report
  app_log() << "    dgate = " << dgate << std::endl;
  app_log() << "    mimg  = " << mimg << std::endl;
  app_log() << "     rcut = " << LR_rc << std::endl;
  app_log() << "     nlat = " << nlat << std::endl;
  app_log() << "    reach = " << reach << std::endl;
  bool ltoo_small = (new_rcut > reach);
  if (ltoo_small)
  {
    const int mlat = 10;
    while (nlat < mlat)
    {
      nlat++;
      reach = LR_rc + nlat*2.0*LR_rc;
      if (reach >= new_rcut) break;
    }
    std::ostringstream msg;
    msg << "nlat = " << ref.getLattice().nlat << " is too small"
        << "  try increase to " << nlat << "\n";
    throw std::runtime_error(msg.str());
  }
  fillFk(ref.getSimulationCell().getKLists());
}

void Screen2D::fillFk(const KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  MaxKshell = KList.kshell.size() - 1;
  Fk_symm.resize(MaxKshell);

  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    Fk_symm[ks] = 0.0;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
      Fk[ki++] = 0.0;
  }
}

Screen2D::mRealType Screen2D::evaluate(mRealType r, mRealType rinv) const
{
  return screen2DKernel<mRealType>(r, dgate, mimg);
}

Screen2D::mRealType Screen2D::findRcut()
{
  const int miter = 1000;
  const mRealType edge = 0.9;
  const mRealType tol = 1e-10;
  mRealType r = LR_rc;
  mRealType rmax = 10*LR_rc;
  mRealType rmin = 0.0;
  int iter = 0;
  mRealType vtail = evaluate(r, 1./r);
  mRealType vedge = evaluate(edge*r, 1./r/edge);
  app_log() << "  Screen2D::findRcut\n";
  while ((iter < miter) and ((vedge<tol) or (vtail>tol)))
  {
    app_log() << iter << " " << rmin << " " << r << " " << rmax << "\n";
    if ((rmax-rmin) < 2*tol) break;
    if (vtail>tol)
      r = (r+rmax)/2.0;
    else if (vedge<tol)
      r = (r+rmin)/2.0;
    else
      break;
    vtail = evaluate(r, 1./r);
    vedge = evaluate(edge*r, 1./r/edge);
    if (vtail>tol)
      rmin = r;
    else if (vedge<tol)
      rmax = r;
    else
      break;
    iter++;
  }
  return r;
}

} // qmcplusplus
