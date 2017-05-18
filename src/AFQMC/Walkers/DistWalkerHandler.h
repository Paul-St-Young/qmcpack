
// -*- C++ -*-
// /**@file WalkerHandler 
//  * @brief Virtual Class for walker handlers. 
//   */

#ifndef QMCPLUSPLUS_AFQMC_DISTWALKERHANDLER_H
#define QMCPLUSPLUS_AFQMC_DISTWALKERHANDLER_H

#include<random>

#include "OhmmsData/libxmldefs.h"
#include "io/hdf_archive.h"

#include"AFQMC/config.h"
#include"AFQMC/Walkers/WalkerHandlerBase.h"
#include<Message/MPIObjectBase.h>
#include "Message/CommOperators.h"
#include "Utilities/NewTimer.h"

#include "AFQMC/Numerics/DenseMatrixOperations.h"


namespace qmcplusplus
{

/*
 * Class that contains and handles walkers.
 * Implements communication, load balancing, and I/O operations.   
 * Walkers are always accessed through the handler.
 */
class DistWalkerHandler: public WalkerHandlerBase 
{

  public:

  /// constructor
  DistWalkerHandler(Communicate* c): WalkerHandlerBase(c,std::string("DistWalkerHandler")),
      min_weight(0.05),max_weight(4.0)
      ,reset_weight(1.0),extra_empty_spaces(10),distribution(0.0,1.0)
      ,head(false),walker_memory_usage(0),tot_num_walkers(0),maximum_num_walkers(0)
  { }

  /// destructor
  ~DistWalkerHandler() {}

  inline int size() { 
    assert( maximum_num_walkers == walkers.size()/walker_size );
    return maximum_num_walkers; 
  } 

  bool restartFromXML(); 
  bool dumpToXML();

  bool restartFromHDF5(int n, hdf_archive&, const std::string&, bool set_to_target); 
  bool dumpToHDF5(hdf_archive&, const std::string&); 

  bool dumpSamplesHDF5(hdf_archive& dump, int nW) {
    std::cerr<<" DistlWalkerHandler:dumpSamplesHDF5() not implemented. Not writing anything. \n";
    return true;
  }

  // reads xml and performs setup
  bool parse(xmlNodePtr cur); 

  // performs setup
  bool setup(int,int,int,MPI_Comm,MPI_Comm,MPI_Comm,myTimer*);

  // cleans state of object. 
  //   -erases allocated memory 
  bool clean(); 

  // called at the beginning of each executable section
  void resetNumberWalkers(int n, bool a=true, ComplexMatrix* S=NULL); 

  inline bool initWalkers(int n) {
   resetNumberWalkers(n,true,&HFMat);
   return true;   
  } 

  inline int numWalkers(bool dummy=false) {
// checking
    int cnt=0;
    for(ComplexSMVector::iterator it=walkers.begin()+data_displ[INFO]; 
       it<walkers.end(); it+=walker_size)
      if(it->real() > 0) cnt++; 
    if(cnt != tot_num_walkers)
      APP_ABORT(" Error in DistWalkerHandler::numWalkers(): Incorrect number of walkers. \n"); 

    if(dummy)
      return size();
    else
      return tot_num_walkers;
  } 

  inline int GlobalPopulation() {
    std::vector<int> res(1);
    res[0]=0;
    if(head)
      res[0] += tot_num_walkers;
    myComm->gsum(res);
    return res[0];
  }

  inline RealType GlobalWeight() {
    std::vector<RealType> res(1);
    res[0]=0;
    if(head)
      for(int i=0; i<maximum_num_walkers; i++)
       if(isAlive(i)) 
        res[0] += std::abs(getWeight(i)); 
    myComm->gsum(res);
    return res[0];
  }

  // load balancing algorithm
  void loadBalance(); 

  // population control algorithm
  void popControl(); 

  void setHF(const ComplexMatrix& HF);

  inline void Orthogonalize(int i) {

    if(walkerType == "closed_shell")  {
      DenseMatrixOperators::GeneralizedGramSchmidt(&(walkers[walker_size*i+data_displ[SM]]),NAEA,NMO,NAEA);
    } else if(walkerType == "collinear") {
      DenseMatrixOperators::GeneralizedGramSchmidt(&(walkers[walker_size*i+data_displ[SM]]),NAEA,NMO,NAEA);
      DenseMatrixOperators::GeneralizedGramSchmidt(&((walkers)[walker_size*i+data_displ[SM]])+NAEA*NMO,NAEA,NMO,NAEB);
    } else if(walkerType == "non-collinear") {
      APP_ABORT("ERROR: non-collinear not implemented in Orthogonalize. \n\n\n");
    }
  }

  inline void Orthogonalize() {
    for(int i=0; i<maximum_num_walkers; i++)
      if( i%ncores_per_TG == core_rank && isAlive(i))
        Orthogonalize(i);
  } 

  //private:
 
  // using std::random for simplicity now
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution;

  enum walker_data { INFO=0, SM=1, WEIGHT=2, ELOC=3, ELOC_OLD=4, OVLP_A=5, OVLP_B=6, OLD_OVLP_A=7, OLD_OVLP_B=8};

  // n is zero-based
  ComplexType* getSM(int n) { if(n>=maximum_num_walkers) {return NULL;} else {return &((walkers)[walker_size*n+data_displ[SM]]);} }
  ComplexType getWeight(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[WEIGHT]];} }
  ComplexType getEloc(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC]];} }
  ComplexType getOldEloc(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC_OLD]];} }
  ComplexType getOvlpAlpha(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_A]];} }
  ComplexType getOvlpBeta(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_B]];} }
  ComplexType getOldOvlpAlpha(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_A]];} }
  ComplexType getOldOvlpBeta(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_B]];} }
  void setWeight(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else {(walkers)[walker_size*n+data_displ[WEIGHT]]=Q;} }
  void setEloc(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC]] = Q;} }
  void setOldEloc(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC_OLD]] = Q;} }
  void setOvlp(int n, ComplexType Q1, ComplexType Q2) { 
    if(n>=maximum_num_walkers) {return;} 
    else { 
      (walkers)[walker_size*n+data_displ[OVLP_A]] = Q1; 
      (walkers)[walker_size*n+data_displ[OVLP_B]] = Q2;
    } 
  }
  void setOldOvlp(int n, ComplexType Q1, ComplexType Q2) { 
    if(n>=maximum_num_walkers) {return;} 
    else { 
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]] = Q1; 
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]] = Q2;
    } 
  }
  void setCurrToOld(int n) { 
    if(n>=maximum_num_walkers) {return;} 
    else { 
      (walkers)[walker_size*n+data_displ[ELOC_OLD]] = (walkers)[walker_size*n+data_displ[ELOC]];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]] = (walkers)[walker_size*n+data_displ[OVLP_A]];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]] = (walkers)[walker_size*n+data_displ[OVLP_B]];
    } 
  }
  ComplexType getEloc2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC]+1];} }
  ComplexType getOldEloc2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC_OLD]+1];} }
  ComplexType getOvlpAlpha2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_A]+1];} }
  ComplexType getOvlpBeta2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_B]+1];} }
  ComplexType getOldOvlpAlpha2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_A]+1];} }
  ComplexType getOldOvlpBeta2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_B]+1];} }
  void setEloc2(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC]+1] = Q;} }
  void setOldEloc2(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC_OLD]+1] = Q;} }
  void setOvlp2(int n, ComplexType Q1, ComplexType Q2) {
    if(n>=maximum_num_walkers) {return;}
    else {
      (walkers)[walker_size*n+data_displ[OVLP_A]+1] = Q1;
      (walkers)[walker_size*n+data_displ[OVLP_B]+1] = Q2;
    }
  }
  void setOldOvlp2(int n, ComplexType Q1, ComplexType Q2) {
    if(n>=maximum_num_walkers) {return;}
    else {
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]+1] = Q1;
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]+1] = Q2;
    }
  }
  void setCurrToOld2(int n) {
    if(n>=maximum_num_walkers) {return;}
    else {
      (walkers)[walker_size*n+data_displ[ELOC_OLD]+1] = (walkers)[walker_size*n+data_displ[ELOC]+1];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]+1] = (walkers)[walker_size*n+data_displ[OVLP_A]+1];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]+1] = (walkers)[walker_size*n+data_displ[OVLP_B]+1];
    }
  }
  void setWalker(int n, ComplexType eloc, ComplexType oa, ComplexType ob) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[ELOC]] = eloc; 
    (walkers)[walker_size*n+data_displ[OVLP_A]] = oa;
    (walkers)[walker_size*n+data_displ[OVLP_B]] = ob;
  }
  void setWalker(int n, ComplexType w0 ,ComplexType eloc) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[WEIGHT]] = w0; 
    (walkers)[walker_size*n+data_displ[ELOC]] = eloc; 
  }
  void setWalker(int n, ComplexType w0 ,ComplexType eloc, ComplexType oa, ComplexType ob) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[WEIGHT]] = w0; 
    (walkers)[walker_size*n+data_displ[ELOC]] = eloc; 
    (walkers)[walker_size*n+data_displ[OVLP_A]] = oa;
    (walkers)[walker_size*n+data_displ[OVLP_B]] = ob;
  }
  void setWalker2(int n, ComplexType eloc, ComplexType oa, ComplexType ob) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[ELOC]+1] = eloc; 
    (walkers)[walker_size*n+data_displ[OVLP_A]+1] = oa;
    (walkers)[walker_size*n+data_displ[OVLP_B]+1] = ob;
  }
  void getOldWalker(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return ;} 
    eloc = (walkers)[walker_size*n+data_displ[ELOC_OLD]]; 
    oa = (walkers)[walker_size*n+data_displ[OLD_OVLP_A]];
    ob = (walkers)[walker_size*n+data_displ[OLD_OVLP_B]];
  } 
  ComplexType* getWalker(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return NULL;} 
    eloc = (walkers)[walker_size*n+data_displ[ELOC]]; 
    oa = (walkers)[walker_size*n+data_displ[OVLP_A]];
    ob = (walkers)[walker_size*n+data_displ[OVLP_B]];
    return &((walkers)[walker_size*n+data_displ[SM]]); 
  }
  ComplexType* getWalker2(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return NULL;} 
    eloc = (walkers)[walker_size*n+data_displ[ELOC]+1]; 
    oa = (walkers)[walker_size*n+data_displ[OVLP_A]+1];
    ob = (walkers)[walker_size*n+data_displ[OVLP_B]+1];
    return &((walkers)[walker_size*n+data_displ[SM]+1]); 
  }
  ComplexType* getWalker(int n, ComplexType& w, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return NULL;}
    w = (walkers)[walker_size*n+data_displ[WEIGHT]];
    eloc = (walkers)[walker_size*n+data_displ[ELOC]];
    oa = (walkers)[walker_size*n+data_displ[OVLP_A]];
    ob = (walkers)[walker_size*n+data_displ[OVLP_B]];
    return &((walkers)[walker_size*n+data_displ[SM]]);
  }

  bool isAlive(int n) {
    return ((n>=maximum_num_walkers)?false:(walkers[walker_size*n+data_displ[INFO]].real()>0)); // for now
  }
  
  void scaleWeight(RealType w0) {
    if(!head) return;
    for(int i=0; i<maximum_num_walkers; i++) 
      (walkers)[walker_size*i+data_displ[WEIGHT]] *= w0;
  } 

  void push_walkers_to_front();

  // Stored data [all assumed std::complex numbers]:
  //   - INFO:                 1  (e.g. alive, init, etc)  
  //   - SlaterMatrix:         NCOL*NROW 
  //        type = 1 for closed_shell
  //        type = 2 for non-closed shell collinear
  //        type = 4 for non-collinear
  //   - weight:               1 
  //   - eloc:                 2 
  //   - eloc_old:             2
  //   - overlap_alpha:        2
  //   - overlap_beta:         2
  //   - old_overlap_alpha:    2
  //   - old_overlap_beta:     2 
  //   Total: 14+NROW*NCOL
  int type, nrow, ncol; 
  int walker_size, data_displ[9], walker_memory_usage; 

  int tot_num_walkers;
  int maximum_num_walkers;

  ComplexType zero = ComplexType(0.0,0.0);
  RealType min_weight, max_weight, reset_weight;

  int extra_empty_spaces; 
  
  ComplexMatrix HFMat; 

  std::vector<int> empty_spots;

  bool head;

  // container with walker data 
  ComplexSMVector walkers; 

  MPI_Comm MPI_COMM_TG_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL_HEADS; 
  int nproc_heads, rank_heads;
  
  std::vector<std::tuple<int,int>> outgoing, incoming; 
  std::vector<int> counts,displ;
  std::vector<char> bufferall;
  std::vector<ComplexType> commBuff;

  myTimer* LocalTimer;
  TimerList_t Timers;

};
}

#endif
