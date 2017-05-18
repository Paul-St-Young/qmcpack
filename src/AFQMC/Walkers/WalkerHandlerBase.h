
// -*- C++ -*-
// /**@file WalkerHandler 
//  * @brief Virtual Class for walker handlers. 
//   */

#ifndef QMCPLUSPLUS_AFQMC_WALKERHANDLERBASE_H
#define QMCPLUSPLUS_AFQMC_WALKERHANDLERBASE_H

#include "OhmmsData/libxmldefs.h"
#include "io/hdf_archive.h"

#include"AFQMC/config.h"
#include <Message/MPIObjectBase.h>
//#include "Message/CommOperators.h"

namespace qmcplusplus
{

/*
 * Class that contains and handles walkers.
 * Implements communication, load balancing, and I/O operations.   
 * Walkers are always accessed through the handler.
 */
class WalkerHandlerBase: public MPIObjectBase, public AFQMCInfo
{

  public:

  /// constructor
  WalkerHandlerBase(Communicate* c, std::string type=""): MPIObjectBase(c),name("")
           ,load_balance_alg("all"),core_rank(0),ncores_per_TG(1),walkerType(type)
  { }

  /// destructor
  ~WalkerHandlerBase() {}


  virtual int size()=0; 

  virtual bool restartFromHDF5(int n, hdf_archive&, const std::string&, bool set_to_target)=0; 
  virtual bool dumpToHDF5(hdf_archive&, const std::string&)=0; 

  virtual bool dumpSamplesHDF5(hdf_archive& dump, int nW)=0; 

  // reads xml and performs setup
  virtual bool parse(xmlNodePtr cur)=0; 

  // performs setup
  virtual bool setup(int a,int b,int c,MPI_Comm, MPI_Comm comm, MPI_Comm,myTimer*)=0;

  // cleans state of object. 
  //   -erases allocated memory 
  virtual bool clean()=0; 

  // called at the beginning of each executable section
  virtual void resetNumberWalkers(int n, bool a=true, ComplexMatrix* S=NULL)=0; 

  virtual  bool initWalkers(int n)=0; 

  virtual int numWalkers(bool dummy=false)=0; 

  virtual int numWalkers2() { return 0; } 

  virtual int GlobalPopulation()=0; 
  virtual RealType GlobalWeight()=0; 

  // load balancing algorithm
  virtual void loadBalance()=0; 

  // population control algorithm
  virtual void popControl()=0; 

  virtual void setHF(const ComplexMatrix& HF)=0;

  virtual void Orthogonalize(int i)=0; 

  virtual void Orthogonalize()=0; 

  virtual ComplexType* getSM(int n)=0; 
  virtual ComplexType getWeight(int n)=0; 
  virtual ComplexType getEloc(int n)=0;
  virtual ComplexType getEloc2(int n)=0; 
  virtual ComplexType getOldEloc(int n)=0; 
  virtual ComplexType getOldEloc2(int n)=0; 
  virtual ComplexType getOvlpAlpha(int n)=0; 
  virtual ComplexType getOvlpAlpha2(int n)=0; 
  virtual ComplexType getOvlpBeta(int n)=0; 
  virtual ComplexType getOvlpBeta2(int n)=0; 
  virtual ComplexType getOldOvlpAlpha(int n)=0; 
  virtual ComplexType getOldOvlpAlpha2(int n)=0; 
  virtual ComplexType getOldOvlpBeta(int n)=0; 
  virtual ComplexType getOldOvlpBeta2(int n)=0; 
  virtual void setWeight(int n, ComplexType Q)=0; 
  virtual void setEloc(int n, ComplexType Q)=0; 
  virtual void setEloc2(int n, ComplexType Q)=0; 
  virtual void setOldEloc(int n, ComplexType Q)=0; 
  virtual void setOldEloc2(int n, ComplexType Q)=0; 
  virtual void setOvlp(int n, ComplexType Q1, ComplexType Q2)=0;
  virtual void setOvlp2(int n, ComplexType Q1, ComplexType Q2)=0; 
  virtual void setOldOvlp(int n, ComplexType Q1, ComplexType Q2)=0; 
  virtual void setOldOvlp2(int n, ComplexType Q1, ComplexType Q2)=0; 
  virtual void setCurrToOld(int n)=0; 
  virtual void setCurrToOld2(int n)=0; 
  virtual bool isAlive(int n)=0; 
  virtual void scaleWeight(RealType w0)=0; 
  virtual void setWalker(int,ComplexType,ComplexType)=0; 
  virtual void setWalker(int,ComplexType,ComplexType,ComplexType)=0; 
  virtual void setWalker(int,ComplexType,ComplexType,ComplexType,ComplexType)=0; 
  virtual void setWalker2(int,ComplexType,ComplexType,ComplexType)=0; 
  virtual ComplexType* getWalker(int,ComplexType&,ComplexType&,ComplexType&)=0; 
  virtual ComplexType* getWalker2(int,ComplexType&,ComplexType&,ComplexType&)=0; 
  virtual ComplexType* getWalker(int,ComplexType&,ComplexType&,ComplexType&,ComplexType&)=0; 
  virtual void getOldWalker(int,ComplexType&,ComplexType&,ComplexType&)=0; 

  // name of the object
  std::string name;

  // type of walker
  std::string walkerType;

  int nwalk_global, nwalk_min, nwalk_max;
  std::vector<int> nwalk_counts_old, nwalk_counts_new;

  // type of load balancing
  std::string load_balance_alg;

  int core_rank;
  int ncores_per_TG; 

};
}

#endif
