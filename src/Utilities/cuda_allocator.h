#ifndef CUDA_ALLOCATOR_H
#define CUDA_ALLOCATOR_H

#include <cuda_runtime_api.h>
#include <malloc.h>

template<typename T> class cuda_allocator;

template<>
class cuda_allocator<void>
{
public:
  typedef void*  pointer;
  typedef const  void* const_pointer;
  typedef void value_type;

  template<class T1> struct rebind 
  { 
    typedef cuda_allocator<T1> other; 
  };
};


template<typename T> 
class cuda_allocator 
{
public:
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef T*        pointer;
  typedef const T*  const_pointer;
  typedef T&        reference;
  typedef const T&  const_reference;
  typedef T         value_type;
  template<typename U> struct rebind { typedef cuda_allocator<U> other; };
  
  cuda_allocator() throw() { }
  cuda_allocator(const cuda_allocator&) throw() { }
  template<typename U> cuda_allocator(const cuda_allocator<U>&) throw() { }
  ~cuda_allocator() throw() { };

  pointer address(reference x) const 
  { return &x; }

  const_pointer address(const_reference x) const 
  { return &x; }
  
  pointer allocate(size_type s, cuda_allocator<void>::const_pointer hint = 0)
  {

  }
  
  void deallocate(pointer p, size_type n)
  {

  }

  size_type max_size() const throw()
  { return (size_type)1 << 32 - 1; }
  
  void construct(pointer p, const T& val)
  { new(static_cast<void*>(p)) T(val);  }
  
  void destroy(pointer p) {
    p->~T();
  }
};


#endif
