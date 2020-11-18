#ifndef PED_MEM_H
#define PED_MEM_H
#include "config.h"
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace pedmod {

/** class to hold allocate and hold memory */
template<class T>
class cache_mem {
  std::unique_ptr<T[]> mem;
  size_t cur_max_threads = 0,
               cur_n_mem = 0,
               cur_size  = 0;

public:
  /// set the size of the memory to hold
  void set_n_mem(size_t n_mem, size_t const max_threads){
    constexpr size_t const mult = cacheline_size() / sizeof(T);
    n_mem  = std::max(n_mem, mult) + mult;
    n_mem  = (n_mem + mult - 1L) / mult;
    n_mem *= mult;

    cur_max_threads = std::max(cur_max_threads, max_threads);
    cur_n_mem       = std::max(cur_n_mem      , n_mem);

    size_t const req_size = cur_n_mem * cur_max_threads;
    if(cur_size < req_size){
      // needs more elements
      mem.reset(new T[req_size]);
      cur_size = req_size;
    }
  }

  /// returns a pointer to the memory to use
  inline T * get_mem(int const thread_num) noexcept {
    return mem.get() + thread_num * cur_n_mem;
  }

  inline T * get_mem() noexcept {
#ifdef _OPENMP
    return get_mem(omp_get_thread_num());
#else
    return get_mem(0);
#endif
  }
};

} // namespace pedmod

#endif
