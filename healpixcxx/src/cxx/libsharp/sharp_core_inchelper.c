#define Tb CONCAT2(Tb,nvec)
#define Y(arg) CONCAT2(arg,nvec)
#include "sharp_core_inc.c"

#if (SHARP_MAXTRANS>MAXJOB_SPECIAL)
#define NJ1 , int njobs
#define NJ2 , njobs
#define Z(arg) CONCAT2(arg,nvec)
#include "sharp_core_inc2.c"
#undef Z
#undef NJ1
#undef NJ2
#endif

#define NJ1
#define NJ2

#if ((MAXJOB_SPECIAL>=1)&&(SHARP_MAXTRANS>=1))
#define njobs 1
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if ((MAXJOB_SPECIAL>=2)&&(SHARP_MAXTRANS>=2))
#define njobs 2
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if ((MAXJOB_SPECIAL>=3)&&(SHARP_MAXTRANS>=3))
#define njobs 3
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if ((MAXJOB_SPECIAL>=4)&&(SHARP_MAXTRANS>=4))
#define njobs 4
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if ((MAXJOB_SPECIAL>=5)&&(SHARP_MAXTRANS>=5))
#define njobs 5
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#if ((MAXJOB_SPECIAL>=6)&&(SHARP_MAXTRANS>=6))
#define njobs 6
#define Z(arg) CONCAT3(arg,nvec,njobs)
#include "sharp_core_inc2.c"
#undef Z
#undef njobs
#endif

#undef NJ1
#undef NJ2

#undef Y
#undef Tb
