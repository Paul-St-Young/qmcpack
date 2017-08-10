//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


/** @file HybridCplxAdoptor.h
 *
 * Adoptor classes to handle complex hybrid orbitals with arbitrary precision
 */
#ifndef QMCPLUSPLUS_HYBRID_CPLX_SOA_ADOPTOR_H
#define QMCPLUSPLUS_HYBRID_CPLX_SOA_ADOPTOR_H

#include <QMCWaveFunctions/BsplineFactory/HybridAdoptorBase.h>
namespace qmcplusplus
{

/** adoptor class to match
 *
 */
template<typename BaseAdoptor>
struct HybridCplxSoA: public BaseAdoptor, public HybridAdoptorBase<typename BaseAdoptor::DataType>
{
  using HybridBase       = HybridAdoptorBase<typename BaseAdoptor::DataType>;
  using ST               = typename BaseAdoptor::DataType;
  using PointType        = typename BaseAdoptor::PointType;
  using SingleSplineType = typename BaseAdoptor::SingleSplineType;

  typename OrbitalSetTraits<ST>::ValueVector_t psi_AO, d2psi_AO;
  typename OrbitalSetTraits<ST>::GradVector_t dpsi_AO;

  using BaseAdoptor::myV;
  using BaseAdoptor::myG;
  using BaseAdoptor::myL;
  using BaseAdoptor::myH;
  using HybridBase::dist_r;
  using HybridBase::dist_dr;
  using HybridBase::df_dr;
  using HybridBase::d2f_dr2;

  HybridCplxSoA(): BaseAdoptor()
  {
    this->AdoptorName="Hybrid"+this->AdoptorName;
    this->KeyWord="Hybrid"+this->KeyWord;
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseAdoptor::resizeStorage(n,nvals);
    HybridBase::resizeStorage(myV.size());
  }

  void bcast_tables(Communicate* comm)
  {
    BaseAdoptor::bcast_tables(comm);
    HybridBase::bcast_tables(comm);
  }

  void reduce_tables(Communicate* comm)
  {
    BaseAdoptor::reduce_tables(comm);
    HybridBase::reduce_atomic_tables(comm);
  }

  bool read_splines(hdf_archive& h5f)
  {
    return BaseAdoptor::read_splines(h5f) && HybridBase::read_splines(h5f);
  }

  bool write_splines(hdf_archive& h5f)
  {
    return BaseAdoptor::write_splines(h5f) && HybridBase::write_splines(h5f);
  }

  inline void flush_zero()
  {
    BaseAdoptor::flush_zero();
    HybridBase::flush_zero();
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const ST smooth_factor=HybridBase::evaluate_v(P,iat,myV);
    if(smooth_factor<0)
    {
      BaseAdoptor::evaluate_v(P,iat,psi);
    }
    else if (smooth_factor==ST(1))
    {
      const PointType& r=P.R[iat];
      BaseAdoptor::assign_v(r,psi);
    }
    else
    {
      const PointType& r=P.R[iat];
      const ST cone(1);
      psi_AO.resize(psi.size());
      BaseAdoptor::assign_v(r,psi_AO);
      BaseAdoptor::evaluate_v(P,iat,psi);
      for(size_t i=0; i<psi.size(); i++)
        psi[i] = psi_AO[i]*smooth_factor + psi[i]*(cone-smooth_factor);
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const ST smooth_factor=HybridBase::evaluate_vgl(P,iat,myV,myG,myL);
    if(smooth_factor<0)
    {
      BaseAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);
    }
    else if(smooth_factor==ST(1))
    {
      const PointType& r=P.R[iat];
      BaseAdoptor::assign_vgl_from_l(r,psi,dpsi,d2psi);
    }
    else
    {
      const PointType& r=P.R[iat];
      const ST cone(1), ctwo(2);
      const ST rinv(1.0/dist_r);
      psi_AO.resize(psi.size());
      dpsi_AO.resize(psi.size());
      d2psi_AO.resize(psi.size());
      BaseAdoptor::assign_vgl_from_l(r,psi_AO,dpsi_AO,d2psi_AO);
      BaseAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);
      for(size_t i=0; i<psi.size(); i++)
      {
        d2psi[i] = d2psi_AO[i]*smooth_factor + d2psi[i]*(cone-smooth_factor)
                 + dot(dpsi[i]-dpsi_AO[i], dist_dr) * df_dr * rinv * ctwo
                 + (psi_AO[i]-psi[i]) * (d2f_dr2 + ctwo * rinv *df_dr);
         dpsi[i] =  dpsi_AO[i]*smooth_factor +  dpsi[i]*(cone-smooth_factor)
                 + (psi[i]-psi_AO[i]) * df_dr * rinv * dist_dr;
          psi[i] =   psi_AO[i]*smooth_factor +   psi[i]*(cone-smooth_factor);
      }
    }
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& vgl)
  {
    if(HybridBase::evaluate_vgh(P,iat,myV,myG,myH))
    {
      const PointType& r=P.R[iat];
      BaseAdoptor::assign_vgl_soa(r,vgl);
    }
    else
      BaseAdoptor::evaluate_vgl_combo(P,iat,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  inline void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    if(HybridBase::evaluate_vgh(P,iat,myV,myG,myH))
    {
      const PointType& r=P.R[iat];
      BaseAdoptor::assign_vgh(r,psi,dpsi,grad_grad_psi);
    }
    else
      BaseAdoptor::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
  }
};

}
#endif
