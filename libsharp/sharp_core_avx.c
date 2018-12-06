#if (!defined(__AVX__)) && defined(__GNUC__) && defined (__x86_64__) && (__GNUC__>=6)
// if we arrive here, we can benefit from an additional AVX version
// #warning entering gcc and x86_64 specific code branch

#define ARCH _avx
#pragma GCC target("avx")
#include "sharp_core_inc0.c"
#undef ARCH

#endif
