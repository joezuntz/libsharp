#define Tb CONCAT2(Tb,nvec)
#define Y(arg) CONCAT2(arg,nvec)
#include "sharp_core_inc.c"

#define Z(arg) CONCAT2(arg,nvec)
#include "sharp_core_inc2.c"
#undef Z

#undef Y
#undef Tb
