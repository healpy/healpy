#ifndef DUCC0_USEFUL_MACROS_H
#define DUCC0_USEFUL_MACROS_H

#if defined(__GNUC__)
#define DUCC0_NOINLINE [[gnu::noinline]]
#define DUCC0_RESTRICT __restrict__
#define DUCC0_PREFETCH_R(addr) __builtin_prefetch(addr);
#define DUCC0_PREFETCH_W(addr) __builtin_prefetch(addr,1);
#elif defined(_MSC_VER)
#define DUCC0_NOINLINE __declspec(noinline)
#define DUCC0_RESTRICT __restrict
#define DUCC0_PREFETCH_R(addr)
#define DUCC0_PREFETCH_W(addr)
#else
#define DUCC0_NOINLINE
#define DUCC0_RESTRICT
#define DUCC0_PREFETCH_R(addr)
#define DUCC0_PREFETCH_W(addr)
#endif

#endif
