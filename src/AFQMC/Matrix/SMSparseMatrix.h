
#ifndef QMCPLUSPLUS_AFQMC_SMSPARSEMATRIX_H
#define QMCPLUSPLUS_AFQMC_SMSPARSEMATRIX_H

#include<tuple>
#include<algorithm>
#include"AFQMC/Utilities/tuple_iterator.hpp"
#include<iostream>
#include<vector>
#include<assert.h>
#include"AFQMC/config.0.h"
#include"Utilities/UtilityFunctions.h"
#include<mpi.h>
#include <sys/time.h>
#include <ctime>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/exceptions.hpp>

#define ASSERT_SPARSEMATRIX 

#if defined(USE_EIGEN)
namespace qmcplusplus
{
}

#else  // In this case, use OhhmsPETE and your sparse matrix class

namespace qmcplusplus
{

// class that implements a sparse matrix in CSR format
template<class T>
class SMSparseMatrix
{
  public:

  template<typename spT> using ShmemAllocator = boost::interprocess::allocator<spT, boost::interprocess::managed_shared_memory::segment_manager>;
  template<typename spT> using SMVector = boost::interprocess::vector<spT, ShmemAllocator<spT>>;


  typedef T            Type_t;
  typedef T            value_type;
  typedef T*           pointer;
  typedef const T*     const_pointer;
  typedef const int*   const_indxPtr;
  typedef int          indxType;
  typedef int*           indxPtr;
  typedef typename SMVector<T>::iterator iterator;
  typedef typename SMVector<T>::const_iterator const_iterator;
  typedef typename SMVector<int>::iterator int_iterator;
  typedef typename SMVector<int>::const_iterator const_int_iterator;
  typedef SMSparseMatrix<T>  This_t;

  SMSparseMatrix<T>():nr(0),nc(0),compressed(false),zero_based(true),head(false),ID(""),SMallocated(false),storage_format(0),vals(NULL),rowIndex(NULL),myrows(NULL),colms(NULL),mutex(NULL),share_buff(NULL),segment(NULL) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  SMSparseMatrix<T>(int n):nr(n),nc(n),compressed(false),zero_based(true),head(false),ID(""),SMallocated(false),storage_format(0),vals(NULL),rowIndex(NULL),myrows(NULL),colms(NULL),mutex(NULL),share_buff(NULL),segment(NULL)
  {
    remover.ID="NULL";
    remover.head=false;
  }

  SMSparseMatrix<T>(int n,int m):nr(n),nc(m),compressed(false),zero_based(true),head(false),ID(""),SMallocated(false),storage_format(0),vals(NULL),rowIndex(NULL),myrows(NULL),colms(NULL),mutex(NULL),share_buff(NULL),segment(NULL) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  ~SMSparseMatrix<T>()
  {
    if(segment != NULL) {
     delete segment;
     boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  }

  // this should probably be disallowed
  SMSparseMatrix(const SMSparseMatrix<T> &rhs)
  {
    APP_ABORT(" Error: SMSparseMatrix copy constructor has been disabled.");
  }

  inline void setup(bool hd, std::string ii, MPI_Comm comm_) {
    head=hd;
    ID=ii;
    remover.ID=ii;
    remover.head=hd;
    comm=comm_;
  }

  inline void barrier() {
    MPI_Barrier(comm);
  }

  inline void reserve(int n, bool allow_reduce = false)
  {
    if(vals==NULL || (vals!=NULL && vals->capacity() < n) || (vals!=NULL && vals->capacity() > n && allow_reduce)) 
      allocate(n,allow_reduce);
    if(head) {
      vals->reserve(n);
      myrows->reserve(n);
      colms->reserve(n); 
      rowIndex->reserve(std::max(nr,nc)+1);
      share_buff->resize(1000*sizeof(std::complex<double>));
    }
    barrier();
  }

  // does not allow grow/shrink
  inline bool allocate_serial(long n)
  {
    if(!head) { SMallocated=true; return true; }
    if(vals!=NULL && vals->capacity() >= n) { SMallocated=true; return true; }

    memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+(2*n+std::max(nr,nc)+1)*sizeof(int)+1000*sizeof(std::complex<double>)+8000;

    try {
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory); 
    } catch(boost::interprocess::interprocess_exception &ex) {
      std::cout<<"\n Found managed_shared_memory segment, removing. CAREFUL WITH PERSISTENT SHM MEMORY. !!! \n\n";
      boost::interprocess::shared_memory_object::remove(ID.c_str());
      segment=NULL;
    }

    if(segment==NULL) {
      try {
        segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory); 
      } catch(boost::interprocess::interprocess_exception &ex) {
        std::cerr<<"Problems setting up managed_shared_memory in SMSparseMatrix." <<std::endl;
        return false;
      }
    }

    try {
          
      alloc_int = new ShmemAllocator<int>(segment->get_segment_manager());
      alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());
      alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());
          
      mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();
          
      share_buff = segment->construct<SMVector<unsigned char>>("share_buff")(*alloc_uchar);
      share_buff->resize(1000*sizeof(std::complex<double>));
          
      rowIndex = segment->construct<SMVector<int>>("rowIndex")(*alloc_int);
      rowIndex->reserve(std::max(nr,nc)+1);
          
      myrows = segment->construct<SMVector<int>>("myrows")(*alloc_int);
      myrows->reserve(n);
          
      colms = segment->construct<SMVector<int>>("colms")(*alloc_int);
      colms->reserve(n);
          
      vals = segment->construct<SMVector<T>>("vals")(*alloc_T);
      vals->reserve(n);
        
    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMSparseMatrix." <<std::endl;
      return false;
    }
    SMallocated=true;
    return true;
  }

  // all processes must call this routine
  inline bool allocate(long n, bool allow_reduce=false)
  {
    bool grow = false;
    uint64_t old_sz = (segment==NULL)?0:(segment->get_size());
    if(SMallocated) {
      if(vals!=NULL && vals->capacity() >= n && !allow_reduce) return true;
      grow = true; 
      if(!head) {  // delay delete call on head in case you need to shrink vector
        delete segment;
        segment=NULL;
      }
    }
    barrier();    
    if(head) { 
      memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+(2*n+std::max(nr,nc)+1)*sizeof(int)+1000*sizeof(std::complex<double>)+8000;

      if(grow) {
        if(memory > old_sz) {
          uint64_t extra = memory - old_sz;
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::grow(ID.c_str(), extra)) {
            std::cerr<<" Error growing shared memory in  SMSparseMatrix::allocate(). \n";
            return false;
          }
        } else {
          segment->destroy<SMVector<T>>("vals");
          segment->destroy<SMVector<int>>("myrows");
          segment->destroy<SMVector<int>>("colms");
          myrows = segment->construct<SMVector<int>>("myrows")(*alloc_T);
          myrows->reserve(n);
          colms = segment->construct<SMVector<int>>("colms")(*alloc_T);
          colms->reserve(n);
          vals = segment->construct<SMVector<T>>("vals")(*alloc_T);
          vals->reserve(n);
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::shrink_to_fit(ID.c_str())) {
            std::cerr<<" Error in shrink_to_fit shared memory in SMSparseMatrix::allocate(). \n";
            return false;
          }
        }

        try {
          segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str());
          vals = segment->find<SMVector<T>>("vals").first;
          colms = segment->find<SMVector<int>>("colms").first;
          myrows = segment->find<SMVector<int>>("myrows").first;
          rowIndex = segment->find<SMVector<int>>("rowIndex").first;
          share_buff = segment->find<SMVector<unsigned char>>("share_buff").first;
          mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
          assert(share_buff != 0);
          assert(mutex != 0);
          assert(vals != 0);
          assert(myrows != 0);
          assert(colms != 0);
          assert(rowIndex != 0);
          share_buff->resize(1000*sizeof(std::complex<double>)); 
          rowIndex->reserve(std::max(nr,nc)+1);
          myrows->reserve(n);
          colms->reserve(n);
          vals->reserve(n);  
        } catch(std::bad_alloc&) {
          std::cerr<<"Problems opening shared memory in SMDenseVector::allocate() ." <<std::endl;
          return false;
        }

      } else { // grow

        try {
          segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
        } catch(boost::interprocess::interprocess_exception &ex) {
          std::cout<<"\n Found managed_shared_memory segment, removing. CAREFUL WITH PERSISTENT SHM MEMORY. !!! \n\n";
          boost::interprocess::shared_memory_object::remove(ID.c_str());
          segment=NULL;
        }
    
        if(segment==NULL) {
          try {
            segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
          } catch(boost::interprocess::interprocess_exception &ex) {
            std::cerr<<"Problems setting up managed_shared_memory in SMSparseMatrix." <<std::endl;
            return false;
          }
        }

        try {

          alloc_int = new ShmemAllocator<int>(segment->get_segment_manager());
          alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());
          alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());

          mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();

          share_buff = segment->construct<SMVector<unsigned char>>("share_buff")(*alloc_uchar);
          share_buff->resize(1000*sizeof(std::complex<double>));

          rowIndex = segment->construct<SMVector<int>>("rowIndex")(*alloc_int);
          rowIndex->reserve(std::max(nr,nc)+1);

          myrows = segment->construct<SMVector<int>>("myrows")(*alloc_int);
          myrows->reserve(n);

          colms = segment->construct<SMVector<int>>("colms")(*alloc_int);
          colms->reserve(n);

          vals = segment->construct<SMVector<T>>("vals")(*alloc_T);
          vals->reserve(n);

        } catch(std::bad_alloc&) {
          std::cerr<<"Problems allocating shared memory in SMSparseMatrix." <<std::endl;
          return false;
        }
      }
    }
    barrier();
    SMallocated=true;
    initializeChildren();
    return true;
  }

  // only call this when all arrays have been allocated and modified  
  inline bool initializeChildren()
  { 
    if(head) return true;
    // delete segment in case this routine is called multiple times.
    // SHM is not removed, just the mapping of the local process. 
    if(segment!=NULL) {
      delete segment;
      segment=NULL;
    }    
    try {
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str()); 
      vals = segment->find<SMVector<T>>("vals").first;
      colms = segment->find<SMVector<int>>("colms").first;
      myrows = segment->find<SMVector<int>>("myrows").first;
      rowIndex = segment->find<SMVector<int>>("rowIndex").first;
      share_buff = segment->find<SMVector<unsigned char>>("share_buff").first;
      mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
      assert(share_buff != 0);
      assert(mutex != 0);
      assert(vals != 0);
      assert(myrows != 0);
      assert(colms != 0);
      assert(rowIndex != 0);
    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMSparseMatrix: initializeChildren() ." <<std::endl;
      return false;
    }
    return true;
  }

  // does not allow grow/shrink, aborts if resizing beyond capacity
  inline void resize_serial(int nnz)
  {
    if(!head) return;
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz))
      APP_ABORT(" Error: Call to SMSparseMatrix::resize_serial without enough capacity. \n");
    if(head) {
      vals->resize(nnz);
      myrows->resize(nnz);
      colms->resize(nnz);
      rowIndex->resize(nr+1);
    }
  } 

  // this routine does not preserve information when allow_reduce=true  
  inline void resize(int nnz, bool allow_reduce=false)
  {
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz) ) {
      allocate(nnz,allow_reduce);
    } else if(vals!=NULL && vals->capacity() > nnz && allow_reduce) {
      // not keeping information if reduction occurs 
      allocate(nnz,allow_reduce);
    }
    if(head) {
      vals->resize(nnz);
      myrows->resize(nnz);
      colms->resize(nnz);
      rowIndex->resize(nr+1);
    }
    barrier();
  }

  inline void clear() { 
    compressed=false;
    zero_based=true;
    if(!SMallocated) return; 
    if(!head) return;
    vals->clear();
    colms->clear();
    myrows->clear();
    rowIndex->clear();
  }

  inline void setDims(int n, int m)
  {
    nr=n;
    nc=m;
    compressed=false;
    zero_based=true;
    clear();
  }

  template<class T1>
  inline bool copyToBuffer(std::vector<T1>& v)
  {
    if( v.size()*sizeof(T1) > share_buff->capacity() )
      return false; 
    {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      std::memcpy( share_buff->data(), v.data(), v.size()*sizeof(T1) );
    }
  }
 
  template<class T1>
  inline bool copyFromBuffer(std::vector<T1>& v)
  {
    if( v.size()*sizeof(T1) > share_buff->capacity() )
      return false; 
    std::memcpy( v.data(), share_buff->data(), v.size()*sizeof(T1) );
  }  

  template<typename T1>
  inline void share(T1* x, int n, bool sender) {
    if(!SMallocated)
      APP_ABORT("Error: Call to SMDenseVector::share with unallocated object. \n");
    assert( sizeof(T1)*n < sizeof(unsigned char)*share_buff->size() );
    if(sender) {
      std::memcpy(&((*share_buff)[0]),x,sizeof(T1)*n);
      barrier();
    } else {
      barrier();
      std::memcpy(x,&((*share_buff)[0]),sizeof(T1)*n);
    }
    barrier();
  }

  template<typename T1>
  inline void share(T1& x, bool sender) {
    if(!SMallocated)
      APP_ABORT("Error: Call to SMDenseVector::share with unallocated object. \n");
    if(sender) {
      std::memcpy(&((*share_buff)[0]),&x,sizeof(T1));
      barrier();
    } else {
      barrier();
      std::memcpy(&x,&((*share_buff)[0]),sizeof(T1));
    }
    barrier();
  }

  template<typename T1>
  inline void share(std::vector<T1>& x, int n, bool sender) {
    if(!SMallocated)
      APP_ABORT("Error: Call to SMDenseVector::share with unallocated object. \n");
    assert( sizeof(T1)*n < sizeof(unsigned char)*share_buff->size() );
    assert( x.size() >= n);
    if(sender) {
      std::memcpy(&((*share_buff)[0]),x.data(),sizeof(T1)*n);
      barrier();
    } else {
      barrier();
      std::memcpy(x.data(),&((*share_buff)[0]),sizeof(T1)*n);
    }
    barrier();
  }

  inline void setCompressed() 
  {
    compressed=true;
  }

  inline bool isCompressed() const
  {
    return compressed;
  }

  inline uint64_t memoryUsage() { return memory; }

  inline int capacity() const
  {
    return (vals!=NULL)?(vals->capacity()):0;
  }
  inline unsigned long size() const
  {
    return (vals!=NULL)?(vals->size()):0;
  }
  inline int rows() const
  {
    return nr;
  }
  inline int cols() const
  {
    return nc;
  }

  inline const_pointer values() const 
  {
    return (vals!=NULL)?(&((*vals)[0])):NULL;
  }

  inline pointer values() 
  {
    return (vals!=NULL)?(&((*vals)[0])):NULL;
  }

  inline const_indxPtr column_data() const 
  {
    return (colms!=NULL)?(&((*colms)[0])):NULL;
  }
  inline indxPtr column_data() 
  {
    return (colms!=NULL)?(&((*colms)[0])):NULL;
  }

  inline const_indxPtr row_data() const 
  {
    return (myrows!=NULL)?(&((*myrows)[0])):NULL;
  }
  inline indxPtr row_data() 
  {
    return (myrows!=NULL)?(&((*myrows)[0])):NULL;
  }

  inline const_indxPtr row_index() const 
  {
    return (rowIndex!=NULL)?(&((*rowIndex)[0])):NULL;
  }
  inline indxPtr row_index() 
  {
    return  (rowIndex!=NULL)?(&((*rowIndex)[0])):NULL;
  }

  inline This_t& operator=(const SMSparseMatrix<T> &rhs) 
  { 
    compressed=rhs.compressed;
    zero_based=rhs.zero_based;
    nr=rhs.nr;
    nc=rhs.nc;
    if(!head) return *this;
    (*vals)=*(rhs.vals);
    (*myrows)=*(rhs.myrows);
    (*colms)=*(rhs.colms);
    (*rowIndex)=*(rhs.rowIndex);
    return *this;
  }  

  inline int find_element(int i, int j) {
    return 0;
  }

  inline Type_t& operator()(int i, int j)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc && compressed); 
#endif
    return (*vals)[find_element(i,j)]; 
  }

  inline Type_t operator()( int i, int j) const
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc && compressed); 
#endif
    return (*vals)[find_element(i,j)]; 
  }

  inline void add(const int i, const int j, const T& v, bool needs_locks=false) 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc);
#endif
    compressed=false;
//    if(needs_locks) {
      {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
        myrows->push_back(i);
        colms->push_back(j);
        vals->push_back(v);
      }
//    } else {
//      if(!head) return;
//      myrows->push_back(i);
//      colms->push_back(j);
//      vals->push_back(v);
//    }
  }

  inline bool remove_repeated_and_compress() 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
#endif

    compressed=true;
    if(!head) return true;
    if(myrows->size() <= 1) return true;

    // first order arrays
    // order along myrows
    int n=myrows->size(); 
    sort_rows(0,n-1);
    if(!std::is_sorted(myrows->begin(),myrows->end())) 
      std::cout<<"ERROR: list is not sorted. \n" <<std::endl;

    // define rowIndex
    rowIndex->resize(nr+1);
    int curr=-1;
    for(int n=0; n<myrows->size(); n++) {
      if( (*myrows)[n] != curr ) {
        int old = curr;
        curr = (*myrows)[n];
        for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
      }
    }
    for(int i=myrows->back()+1; i<rowIndex->size(); i++)
      (*rowIndex)[i] = vals->size();
   
    // order within each rowIndex block
    for(int k=0; k<nr; k++) {
      if((*rowIndex)[k] == (*rowIndex)[k+1]) continue;       
      sort_colms((*rowIndex)[k],(*rowIndex)[k+1]-1);
    } 
   
    int_iterator first_r=myrows->begin(), last_r=myrows->end(); 
    int_iterator first_c=colms->begin(), last_c=colms->end(); 
    iterator first_v=vals->begin(), last_v = vals->end(); 
    int_iterator result_r = first_r;
    int_iterator result_c = first_c;
    iterator result_v = first_v;

    while ( ( (first_r+1) != last_r)  )
    {
      ++first_r; 
      ++first_c; 
      ++first_v; 
      if (!( (*result_r == *first_r) && (*result_c == *first_c) ) ) { 
        *(++result_r)=*first_r;
        *(++result_c)=*first_c;
        *(++result_v)=*first_v;
      } else {
        if( std::abs(*result_v - *first_v) > 1e-8 ) { //*result_v != *first_v) {
          app_error()<<" Error in call to SMSparseMatrix::remove_repeate_and_compressd. Same indexes with different values. \n";
          app_error()<<"ri, ci, vi: "
                     <<*result_r <<" "
                     <<*result_c <<" "
                     <<*result_v <<" "
                     <<"rj, cj, vj: "
                     <<*first_r <<" "
                     <<*first_c <<" "
                     <<*first_v <<std::endl;
          return false;
        }
      }
    }
    ++result_r;
    ++result_c;
    ++result_v;

    int sz1 = std::distance(myrows->begin(),result_r); 
    int sz2 = std::distance(colms->begin(),result_c); 
    int sz3 = std::distance(vals->begin(),result_v); 
    if(sz1 != sz2 || sz1 != sz2) {
      std::cerr<<"Error: Different number of erased elements in SMSparseMatrix::remove_repeate_and_compressed. \n" <<std::endl;  
      return false;
    }
  
    myrows->resize(sz1);
    colms->resize(sz1);
    vals->resize(sz1);

    // define rowIndex
    curr=-1;
    for(int n=0; n<myrows->size(); n++) {
      if( (*myrows)[n] != curr ) {
        int old = curr;
        curr = (*myrows)[n];
        for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
      }
    }
    for(int i=myrows->back()+1; i<rowIndex->size(); i++)
      (*rowIndex)[i] = vals->size();

    return true;
  }

  // implement parallel compression by using serial sorts on segments of the vector
  // and then either 1. using a parallel in_place_merge of the sorted segments
  // or using serial iterative (every 2 adjacent segments are joined) in_place_merge where every 
  // iteration half of the tasks do nothing (compared to the last iteration)  
  //
  // Possible algorithm:
  // 00. sort rows in parallel, generate rowIndex in serial, round robin for columns
  // 0. only processor numbers in powers of 2 work, so the largest power of 2 in the TG do work, 
  //    everybody else (from 2^lmax to nproc) do nothing 
  // 1. divide the array in Nseg = 2^lmax equalish segments. each working core performs serial in_place_sort on each segment. You might want to implement std::sort with boost::iterator_facade, which should be faster than your pure quick_sort
  // 2. now you have Nseg consecutive sorted subsequences. Now iteratively inplace_merge neighboring subsequences in parallel. Iterate from i=0:lmax and at each iteration 2^(i+1) processors fo work on a given inplace_merge. After each iteration, we adjust the boundaries of the new places to be merged and the number of processors involved.     
  //   - the parallel inplace_merge proceeds by performing iteratively block swap calls to break up the 2 sorted subsets into 2^(i+1) sorted consecutive subsets which can be inplace_merged serially by each working core. 
  //   While it is not ideal, you can use MPI_Barrier() to synchronize the work until you find something better. 
  //
  inline void compress_parallel(MPI_Comm local_comm)
  {

    double t1,t2,t3,t4;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t1 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;

    int npr,rank;
    MPI_Comm_rank(local_comm,&rank); 
    MPI_Comm_size(local_comm,&npr); 

    MPI_Barrier(local_comm);

    gettimeofday(&tv, NULL);
    t2 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time waiting: " <<t2-t1 <<std::endl;

    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
    if(vals->size() == 0) return;

    // sort equal segments in parallel
    std::vector<int> pos(npr+1);     
    // a core processes elements from pos[rank]-pos[rank+1]
    FairDivide(vals->size(),npr,pos); 
    
    std::sort(make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank],colms->begin()+pos[rank],vals->begin()+pos[rank]),
              make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank+1],colms->begin()+pos[rank+1],vals->begin()+pos[rank+1]),
              [](std::tuple<indxType, indxType, value_type> const& a, std::tuple<indxType, indxType, value_type> const& b){return std::get<0>(a) < std::get<0>(b) || (!(std::get<0>(b) < std::get<0>(a)) && std::get<1>(a) < std::get<1>(b));} 
             );
    MPI_Barrier(local_comm);

    gettimeofday(&tv, NULL);
    t1 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time sorting: " <<t1-t2 <<std::endl;


    // merging in serial for now 
    // only head merges and creates index table
    if(head) {

      gettimeofday(&tv, NULL);
      t4 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;

      int mid=1;
      int end=2;
      for(int i=0; i<(npr-1); i++) {
        std::inplace_merge( make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[0],colms->begin()+pos[0],vals->begin()+pos[0]),
                            make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[mid],colms->begin()+pos[mid],vals->begin()+pos[mid]),
                            make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[end],colms->begin()+pos[end],vals->begin()+pos[end]),
                            [](std::tuple<indxType, indxType, value_type> const& a, std::tuple<indxType, indxType, value_type> const& b){return std::get<0>(a) < std::get<0>(b) || (!(std::get<0>(b) < std::get<0>(a)) && std::get<1>(a) < std::get<1>(b));}
            );
        mid++;
        end++;

        gettimeofday(&tv, NULL);
        t3 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
        app_log()<<" Time merging: " <<i <<" " <<t3-t4 <<std::endl;
        t4=t3;
      }


      //if(!std::is_sorted(myrows->begin(),myrows->end())) 
      //  std::cout<<"ERROR: list is not sorted. \n" <<std::endl;

      // define rowIndex
      rowIndex->resize(nr+1);
      int curr=-1; 
      for(int n=0; n<myrows->size(); n++) {
        if( (*myrows)[n] != curr ) {
          int old = curr;
          curr = (*myrows)[n];
          for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
        }
      }
      for(int i=myrows->back()+1; i<rowIndex->size(); i++)
        (*rowIndex)[i] = vals->size();
    
    }
    MPI_Barrier(local_comm);
    compressed=true;   

    gettimeofday(&tv, NULL);
    t2 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time merging + indexing: " <<t2-t1 <<std::endl;
  }

  inline void compress() 
  {
    if(!(myrows->size() == colms->size() && myrows->size() == vals->size()))
      std::cout<<"problem: " <<myrows->size() <<" " <<colms->size() <<" " <<vals->size() <<std::endl; 
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
#endif
    if(!head) { compressed=true; return; }
    if(vals->size() == 0) return;

    // order along myrows
    int n=myrows->size(); 
    sort_rows(0,n-1);
    if(!std::is_sorted(myrows->begin(),myrows->end())) 
      std::cout<<"ERROR: list is not sorted. \n" <<std::endl;

    // define rowIndex
    rowIndex->resize(nr+1);
    int curr=-1;
    for(int n=0; n<myrows->size(); n++) {
      if( (*myrows)[n] != curr ) {
        int old = curr;
        curr = (*myrows)[n];
        for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
      }
    }
    for(int i=myrows->back()+1; i<rowIndex->size(); i++)
      (*rowIndex)[i] = vals->size();
   
    // order within each rowIndex block
    for(int k=0; k<nr; k++) {
      if((*rowIndex)[k] == (*rowIndex)[k+1]) continue;       
      sort_colms((*rowIndex)[k],(*rowIndex)[k+1]-1);
    } 
   
    compressed=true;
  }

  void sort_rows(int left, int right) {
    int i = left, j = right;
    auto pivot = (*myrows)[(left + right) / 2];

    /* partition */
    while (i <= j) {
      while ((*myrows)[i] < pivot)
        i++;
      while ((*myrows)[j] > pivot)
        j--;
      if (i <= j) {
        std::swap((*myrows)[i],(*myrows)[j]);
        std::swap((*colms)[i],(*colms)[j]);
        std::swap((*vals)[i++],(*vals)[j--]);
      }
    };

    /* recursion */
    if (left < j)
      sort_rows(left, j);
    if (i < right)
      sort_rows(i, right);
  }

  void sort_colms(int left, int right) {
    int i = left, j = right;
    auto pivot = (*colms)[(left + right) / 2];

    /* partition */
    while (i <= j) {
      while ((*colms)[i] < pivot)
        i++;
      while ((*colms)[j] > pivot)
        j--;
      if (i <= j) {
        std::swap((*colms)[i],(*colms)[j]);
        std::swap((*vals)[i++],(*vals)[j--]);
      }
    };

    /* recursion */
    if (left < j)
      sort_colms(left, j);
    if (i < right)
      sort_colms(i, right);
  }
  
  inline void transpose(bool paral=false) {
    if(!head) return;
    if(!SMallocated) return;
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
    for(int_iterator itR=myrows->begin(),itC=colms->begin(); itR!=myrows->end(); ++itR,++itC)
      std::swap(*itR,*itC);
    std::swap(nr,nc); 
    compress();
  }

  inline void transpose_parallel(MPI_Comm local_comm)   
  {
    if(!SMallocated) return;
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
    if(head) {
      // can parallelize this if you want
      for(int_iterator itR=myrows->begin(),itC=colms->begin(); itR!=myrows->end(); ++itR,++itC)
        std::swap(*itR,*itC);
    }
    std::swap(nr,nc);
    compress_parallel(local_comm);
  }

  inline void initFroms1D(std::vector<std::tuple<IndexType,RealType> >& V, bool sorted)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(nr==1);
#endif
    if(!head) { compressed=true; return; }
    if(!sorted) 
      //std::sort(V.begin(),V.end(),my_sort);
      std::sort(V.begin(),V.end(), [](const std::tuple<IndexType,RealType>& lhs, const std::tuple<IndexType,RealType>& rhs){return (bool)(std::get<0>(lhs) < std::get<0>(rhs));} );
    myrows->clear();
    rowIndex->clear();
    vals->clear();
    colms->clear();
  
    int nnz=V.size();
    myrows->resize(nnz);
    vals->resize(nnz);
    colms->resize(nnz);
    rowIndex->resize(nr+1);
    (*rowIndex)[0]=0;
    for(int i=0; i<V.size(); i++) {
      if( std::is_same<T,std::complex<double> >::value  ) {
        (*vals)[i] = std::complex<double>(std::get<1>(V[i]),0.0);
      } else {
       (*vals)[i] = static_cast<T>(std::get<1>(V[i]));
      }
      (*myrows)[i] = 0; 
      (*colms)[i] = std::get<0>(V[i]); 
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nc);
#endif
    }
    (*rowIndex)[1]=V.size();
    compressed=true;
  }

  inline void initFroms1D(std::vector<s1D<std::complex<RealType> > >& V, bool sorted)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(nr==1);
#endif
    if(!head) { compressed=true; return; }
    if(!sorted) 
      //std::sort(V.begin(),V.end(),my_sort);
      std::sort(V.begin(),V.end(), [](const std::tuple<IndexType,RealType>& lhs, const std::tuple<IndexType,RealType>& rhs){return (bool)(std::get<0>(lhs) < std::get<0>(rhs));} );
    myrows->clear();
    rowIndex->clear();
    vals->clear();
    colms->clear();
  
    int nnz=V.size();
    myrows->resize(nnz);
    vals->resize(nnz);
    colms->resize(nnz);
    rowIndex->resize(nr+1);
    (*rowIndex)[0]=0;
    for(int i=0; i<V.size(); i++) {
      if( std::is_same<T,std::complex<double> >::value  ) {
        (*vals)[i] = std::get<1>(V[i]);
      } else {
       assert(false);
      }
      (*myrows)[i] = 0; 
      (*colms)[i] = std::get<0>(V[i]); 
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nc);
#endif
    }
    (*rowIndex)[1]=V.size();
    compressed=true;
  }

  inline void initFroms2D(std::vector<s2D<std::complex<RealType> > >&  V, bool sorted)
  {
    if(!head) { compressed=true; return; }
    if(!sorted) 
      //std::sort(V.begin(),V.end(),my_sort);
      std::sort(V.begin(),V.end(), [](const std::tuple<IndexType,RealType>& lhs, const std::tuple<IndexType,RealType>& rhs){return (bool)(std::get<0>(lhs) < std::get<0>(rhs));} );
    myrows->clear();
    rowIndex->clear();
    vals->clear();
    colms->clear();

    int nnz=V.size();
    myrows->resize(nnz);
    vals->resize(nnz);
    colms->resize(nnz);
    rowIndex->resize(nr+1);
    for(int i=0; i<V.size(); i++) {
      if( std::is_same<T,std::complex<double> >::value  ) {
        (*vals)[i] = std::get<2>(V[i]);
      } else {
       assert(false);
      }
      (*myrows)[i] = std::get<0>(V[i]);
      (*colms)[i] = std::get<1>(V[i]);
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nr);
      assert(std::get<1>(V[i]) >= 0 && std::get<1>(V[i]) < nc);
#endif
    }
    int curr=-1;
    for(int n=0; n<myrows->size(); n++) {
      if( (*myrows)[n] != curr ) {
        int old = curr;
        curr = (*myrows)[n];
        for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
      }
    }
    for(int i=myrows->back()+1; i<rowIndex->size(); i++)
      (*rowIndex)[i] = vals->size();
    compressed=true;
  }

  inline void initFroms2D(std::vector<s2D<RealType> >&  V, bool sorted)
  {
    if(!head) { compressed=true; return; }
    if(!sorted) 
      //std::sort(V.begin(),V.end(),my_sort);
      //std::sort(V.begin(),V.end(), [](const std::tuple<IndexType,RealType>& lhs, const std::tuple<IndexType,RealType>& rhs){return (bool)(std::get<0>(lhs) < std::get<0>(rhs));} );
      std::sort(V.begin(),V.end(), [](const s2D<RealType>& lhs, const s2D<RealType>& rhs){return (bool)(std::get<0>(lhs) < std::get<0>(rhs));} );
    myrows->clear();
    rowIndex->clear();
    vals->clear();
    colms->clear();

    int nnz=V.size();
    myrows->resize(nnz);
    vals->resize(nnz);
    colms->resize(nnz);
    rowIndex->resize(nr+1);
    for(int i=0; i<V.size(); i++) {
      if( std::is_same<T,std::complex<double> >::value  ) {
        (*vals)[i] = std::complex<double>(std::get<2>(V[i]),0.0);
      } else {
       (*vals)[i] = static_cast<T>(std::get<2>(V[i]));
      }
      (*myrows)[i] = std::get<0>(V[i]);
      (*colms)[i] = std::get<1>(V[i]);
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nr);
      assert(std::get<1>(V[i]) >= 0 && std::get<1>(V[i]) < nc);
#endif
    }
    int curr=-1;
    for(int n=0; n<myrows->size(); n++) {
      if( (*myrows)[n] != curr ) {
        int old = curr;
        curr = (*myrows)[n];
        for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
      }
    }
    for(int i=myrows->back()+1; i<rowIndex->size(); i++)
      (*rowIndex)[i] = vals->size();
    compressed=true;
  }

  inline void check()
  {
    if(!head) return; 
    for(int i=0; i<rowIndex->size()-1; i++)
    {
      if((*rowIndex)[i+1] < (*rowIndex)[i]) std::cout<<"Error: SMSparseMatrix::check(): rowIndex-> \n" <<std::endl; 
  
    }

  }

  inline SMSparseMatrix<T>& operator*=(const RealType rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline SMSparseMatrix<T>& operator*=(const std::complex<RealType> rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline void toZeroBase() {
    if(!head) return; 
    if(zero_based) return;
    zero_based=true;
    for (int& i : colms ) i--; 
    for (int& i : myrows ) i--; 
    for (int& i : rowIndex ) i--; 
  }

  inline void toOneBase() {
    if(!head) return; 
    if(!zero_based) return;
    zero_based=false;
    for (int& i : colms ) i++; 
    for (int& i : myrows ) i++; 
    for (int& i : rowIndex ) i++; 
  }

  friend std::ostream& operator<<(std::ostream& out, const SMSparseMatrix<T>& rhs)
  {
    for(int i=0; i<rhs.vals->size(); i++)
      out<<"(" <<(*(rhs.myrows))[i] <<"," <<(*(rhs.colms))[i] <<":" <<(*(rhs.vals))[i] <<")\n"; 
    return out;
  }

  friend std::istream& operator>>(std::istream& in, SMSparseMatrix<T>& rhs)
  {
    if(!rhs.head) return in;
    T v;
    int c,r;
    in>>r >>c >>v;
    (rhs.vals)->push_back(v); 
    (rhs.myrows)->push_back(r);
    (rhs.colms)->push_back(c);
    return in;
  }

  // this is ugly, but I need to code quickly 
  // so I'm doing this to avoid adding hdf5 support here 
  inline SMVector<T>* getVals() const { return vals; } 
  inline SMVector<int>* getRows() const { return myrows; }
  inline SMVector<int>* getCols() const { return colms; }
  inline SMVector<int>* getRowIndex() const { return rowIndex; }

  inline iterator vals_begin() { assert(vals!=NULL); return vals->begin(); } 
  inline int_iterator rows_begin() { assert(myrows!=NULL); return myrows->begin(); }
  inline int_iterator cols_begin() { assert(colms!=NULL); return colms->begin(); }
  inline int_iterator rowIndex_begin() { assert(rowIndex!=NULL); return rowIndex->begin(); }
  inline const_iterator vals_begin() const { return vals->begin(); } 
  inline const_int_iterator cols_begin() const { assert(colms!=NULL); return colms->begin(); }
  inline const_int_iterator rowIndex_begin() const { assert(rowIndex!=NULL); return rowIndex->begin(); }
  inline const_iterator vals_end() const { assert(vals!=NULL); return vals->end(); } 
  inline const_int_iterator rows_end() const { assert(myrows!=NULL); return myrows->end(); }
  inline const_int_iterator cols_end() const { assert(colms!=NULL); return colms->end(); }
  inline const_int_iterator rowIndex_end() const { assert(rowIndex!=NULL); return rowIndex->end(); }
  inline iterator vals_end() { assert(vals!=NULL); return vals->end(); } 
  inline int_iterator rows_end() { assert(myrows!=NULL); return myrows->end(); }
  inline int_iterator cols_end() { assert(colms!=NULL); return colms->end(); }
  inline int_iterator rowIndex_end() { assert(rowIndex!=NULL); return rowIndex->end(); }

  inline bool isAllocated() {
    return (SMallocated)&&(vals!=NULL);
  }

  void setRowsFromRowIndex()
  {
    if(!head) return;
    int shift = zero_based?0:1;
    myrows->resize(vals->size());
    for(int i=0; i<nr; i++)
     for(int j=(*rowIndex)[i]; j<(*rowIndex)[i+1]; j++)
      (*myrows)[j]=i+shift;
  }

  bool zero_base() const { return zero_based; }
  int row_max() const { return max_in_row; }
  int format() const { return storage_format; } 

  private:

  boost::interprocess::interprocess_mutex *mutex;
  bool compressed;
  int nr,nc;
  SMVector<T> *vals;
  SMVector<int> *colms,*myrows,*rowIndex;
  SMVector<unsigned char> *share_buff;
  bool head;
  std::string ID; 
  bool SMallocated;
  bool zero_based;
  int storage_format; // 0: CSR, 1: Compressed Matrix (ESSL) 
  int max_in_row; 
  uint64_t memory=0;

  //_mySort_snD_ my_sort;

  boost::interprocess::managed_shared_memory *segment;
  ShmemAllocator<T> *alloc_T;
  ShmemAllocator<boost::interprocess::interprocess_mutex> *alloc_mutex;
  ShmemAllocator<int> *alloc_int;
  ShmemAllocator<unsigned char> *alloc_uchar;

  // using MPI for barrier calls until I find solution
  MPI_Comm comm;
 
  struct shm_remove
  {
    bool head;
    std::string ID; 
    shm_remove() {
      if(head) boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
    ~shm_remove(){
      if(head) boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  } remover;

};


}

#endif
#endif
