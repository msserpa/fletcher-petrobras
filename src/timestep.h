#ifndef _TIMESTEP
#define _TIMESTEP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "map.h"
#include "source.h"
#include "derivatives.h"


void Compare(int sx, int sy, int sz, int bord,
               float * restrict pp, float * restrict pback, float * restrict out);

// Propagate: using Fletcher's equations, propagate waves one dt,
//            either forward or backward in time


void Propagate(int sx, int sy, int sz, int bord,
	       float dx, float dy, float dz, float dt, int it, 
	       float * restrict ch1dxx, float * restrict ch1dyy, float * restrict ch1dzz, 
	       float * restrict ch1dxy, float * restrict ch1dyz, float * restrict ch1dxz, 
	       float * restrict v2px, float * restrict v2pz, float * restrict v2sz, float * restrict v2pn,
	       float * restrict pp, float * restrict pc, float * restrict qp, float * restrict qc);


// TimeForward: swap array pointers on time forward array propagation


void TimeForward(float **pp, float **pc, float **qp, float **qc);

#endif
