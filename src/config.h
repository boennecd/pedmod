#ifndef PEDMOD_CONFIG_H
#define PEDMOD_CONFIG_H

#ifdef DO_CHECKS
#define PEDMOD_NOEXCEPT
#else
#define PEDMOD_NOEXCEPT noexcept
#endif

namespace pedmod {

inline constexpr size_t cacheline_size(){
  return 128L;
}

} // namespace pedmod

#endif
