//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_QMCDRIFTOPERATORS_H
#define QMCPLUSPLUS_QMCDRIFTOPERATORS_H
#include "type_traits/scalar_traits.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/RandomSeqGenerator.h"
namespace qmcplusplus
{

template<class T, class TG, unsigned D>
inline T getDriftScale(T tau, const ParticleAttrib<TinyVector<TG,D> >& ga)
{
  T vsq=Dot(ga,ga);
  return (vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
}

template<class T, class TG, unsigned D>
inline T getDriftScale(T tau, const TinyVector<TG,D>& qf)
{
  T vsq=dot(qf,qf);
  return (vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
}

/** evaluate a drift with a real force
 * @param tau timestep
 * @param qf quantum force
 * @param drift
 */
template<class Tt, class TG, class T, unsigned D>
inline void getScaledDrift(Tt tau, const TinyVector<TG,D>& qf, TinyVector<T,D>& drift)
{
  T vsq=dot(qf,qf);
  vsq= (vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
  drift=vsq*qf;
}

template<class Tt, class TG, class T, unsigned D>
inline void getScaledDrift(Tt tau, const TinyVector<std::complex<TG>,D>& qf, TinyVector<T,D>& drift)
{
  convert(qf,drift);
  T vsq = OTCDot<TG,TG,D>::apply(drift,drift);
  T sc  = (-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq;
  // restore to naive drift if denominator is too small
  if (vsq < std::numeric_limits<T>::epsilon()) sc = tau; 
  drift *= sc;
}

/** scale drift
 * @param tau_au timestep au
 * @param qf quantum forces
 * @param drift scaled quantum forces
 * @param return correction term
 *
 * Assume, mass=1
 */
template<class T, class T1, unsigned D>
inline T setScaledDriftPbyPandNodeCorr(T tau,
                                       const ParticleAttrib<TinyVector<T1,D> >& qf,
                                       ParticleAttrib<TinyVector<T,D> >& drift)
{
  T norm=0.0, norm_scaled=0.0, tau2=tau*tau, vsq;
  for(int iat=0; iat<qf.size(); ++iat)
  {
    convert(dot(qf[iat],qf[iat]),vsq);
    //T vsq=dot(qf[iat],qf[iat]);
    T sc=(vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
    norm_scaled+=vsq*sc*sc;
    norm+=vsq*tau2;
    drift[iat]=qf[iat]*T1(sc);
  }
  return std::sqrt(norm_scaled/norm);
}

/** scale drift
 * @param tau_au timestep au
 * @param massinv 1/m per particle
 * @param qf quantum forces
 * @param drift scaled quantum forces
 * @param return correction term
 *
 * Fill the drift vector one particle at a time (pbyp).
 * The norm of the drift vector is limited by Umrigar's re-scaled drift.
 *
 * T is likely RealType for mass and drift
 * T1 may be ComplexType for wavefunction
 * D is likely int for the number of spatial dimensions
 */
template<class T, class T1, unsigned D>
inline T setScaledDriftPbyPandNodeCorr(T tau_au, const std::vector<T>& massinv,
                                       const ParticleAttrib<TinyVector<T1,D> >& qf,
                                       ParticleAttrib<TinyVector<T,D> >& drift)
{
  T vsq, tau, sc; // temp. variables to be assigned
  T norm2=0.0, norm2_scaled=0.0; // variables to be accumulated
  for(int iat=0; iat<massinv.size(); ++iat)
  {
    tau = tau_au*massinv[iat]; // !!!! assume timestep is scaled by mass

    // drift multiplier of Umrigar, JCP 99, 2865 (1993); eq. (33) * tau
    // keep real part of grad_psi_over_psi
    convert(qf[iat],drift[iat]);
    vsq = dot(drift[iat],drift[iat]);

    // use naive drift if vsq may cause numerical instability in the denominator
    vsq = (vsq < std::numeric_limits<T>::epsilon()) ? tau : vsq; // YY: does this ever happen?

    // drift multiplier of Umrigar, JCP 99, 2865 (1993); eq. (34) * tau
    sc = ((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
    drift[iat] *= sc;

    norm2_scaled += vsq*sc*sc; // accumulate scaled drift norm^2
    norm2 += vsq*tau*tau;      // accumulate naive drift norm^2
  }
  T node_corr = std::sqrt(norm2_scaled/norm2);
  return node_corr;
}

/** get scaled drift vector, allow a different timestep for each particle
 * @param a_vec a list of Umrigar "a" parameter, one for each particle
 * @param grad_vec a list of wavefunction gradients, one for each particle
 * @param drift_vec scaled importance-sampling drift vectors, one for each particle
 *
 * The naive drift is drift = tau*grad, where grad is grad_psi_over_psi, 
 *  namely the log derivative of the guiding wavefunction; tau is timestep.
 * The naive drift diverges at a node, causing persistent configurations.
 * Fill the drift vector one particle at a time.
 * The norm of the drift vector should be limited in two ways:
 *  1. Umrigar: suppress drift divergence to mimic wf divergence
 *  2. Ceperley: limit max drift rate to diffusion rate -> set Umrigar "a" parameter
 * The choice of drift vector does not affect VMC correctness
 *  so long as the proposal probabilities are correctly calculated.
 * The choice of drift vector changes the DMC time-step error v.s. time step, but
 *  does not affect the tau->0 DMC result.
 *
 * fdtype should be either float or double
 * rctype should be either real or complex
 * ndim should be the number of dimensions
 */
template<class fdtype, class rctype, unsigned ndim>
void fill_umrigar_drift(
 const std::vector<fdtype>& tau_vec, // time step
 const std::vector<fdtype>& a_vec,   // Umrigar "a"
 const ParticleAttrib<TinyVector<rctype,ndim> >& grad_vec, // dlogpsi/dr
 ParticleAttrib<TinyVector<fdtype,ndim> >& drift_vec)
{
  fdtype vsq, sc;  // temp. variables to be assigned
  for(int iat=0; iat<tau_vec.size(); ++iat)
  {
    // save real part of grad_psi_over_psi in drift
    convert(grad_vec[iat],drift_vec[iat]); 

    // calculate real gradient magnitude^2
    vsq = dot(drift_vec[iat],drift_vec[iat]); 
    if (vsq<0) APP_ABORT("negative square"); // this should never happen

    // use naive drift if vsq may cause numerical instability in the denominator
    // YY: does this ever happen?
    vsq = (vsq < std::numeric_limits<fdtype>::epsilon()) ? tau_vec[iat] : vsq; 

    // drift multiplier of Umrigar, JCP 99, 2865 (1993); eq. (34) * tau
    sc = (-1.0+std::sqrt(1.0+2.0*a_vec[iat]*vsq*tau_vec[iat]))/(a_vec[iat]*vsq);

    // scale drift
    drift_vec[iat] *= sc;
  }
}

/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf real quantum forces
 * @param drift drift
 */
template<class T, class TG, unsigned D>
inline void setScaledDrift(T tau,
                           const ParticleAttrib<TinyVector<TG,D> >& qf,
                           ParticleAttrib<TinyVector<T,D> >& drift)
{
  T s = getDriftScale(tau,qf);
  PAOps<T,D,TG>::scale(s,qf,drift);
}

/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf real quantum forces
 * @param drift drift
 */
template<class T, unsigned D>
inline void setScaledDrift(T tau,
                           ParticleAttrib<TinyVector<T,D> >& drift)
{
  T s = getDriftScale(tau,drift);
  drift *=s;
}


/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf complex quantum forces
 * @param drift drift
 */
template<class T, class TG, unsigned D>
inline void setScaledDrift(T tau,
                           const ParticleAttrib<TinyVector<std::complex<TG>,D> >& qf,
                           ParticleAttrib<TinyVector<T,D> >& drift)
{
  T s = getDriftScale(tau,qf);
  PAOps<T,D,TG>::scale(s,qf,drift);
}

template<class T, class TG, unsigned D>
inline void assignDrift(T s,
                        const ParticleAttrib<TinyVector<TG,D> >& ga,
                        ParticleAttrib<TinyVector<T,D> >& da)
{
  PAOps<T,D,TG>::scale(s,ga,da);
}

template<class T, class TG, unsigned D>
inline void assignDrift(T s,
                        const ParticleAttrib<TinyVector<std::complex<TG>,D> >& ga,
                        ParticleAttrib<TinyVector<T,D> >& da)
{
  PAOps<T,D,TG>::scale(s,ga,da);
}

template<class T, class T1, unsigned D>
inline void assignDrift(T tau_au, const std::vector<T>& massinv,
                                       const ParticleAttrib<TinyVector<T1,D> >& qf,
                                       ParticleAttrib<TinyVector<T,D> >& drift)
{
  for(int iat=0; iat<massinv.size(); ++iat)
  {
    T tau=tau_au*massinv[iat];
    // naive drift "tau/mass*qf" can diverge
    getScaledDrift(tau,qf[iat],drift[iat]);
  }
}

}
#endif
