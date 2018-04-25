#ifndef _DERIVATIVES
#define _DERIVATIVES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "map.h"
#include "source.h"


// eight order finite differences coefficients of the cross second derivative


#define L11 0.64                    // L1*L1
#define L12 -0.16                   // L1*L2
#define L13 0.03047619047619047618  // L1*L2
#define L14 -0.00285714285714285713 // L1*L4
#define L22 0.04                    // L2*L2
#define L23 -0.00761904761904761904 // L2*L3
#define L24 0.00071428571428571428  // L2*L4
#define L33 0.00145124716553287981  // L3*L3
#define L34 -0.00013605442176870748 // L3*L4
#define L44 0.00001275510204081632  // L4*L4


// eight order finite differences coefficients of the second derivative


#define K0 -2.84722222222222222222  // -205/72
#define K1  1.6                     // 8/5
#define K2 -0.2                     // -1/5
#define K3  0.02539682539682539682  // 8/315
#define K4 -0.00178571428571428571  // -1/560


// Der2: computes second derivative


float Der2(float *p, int i, int s, float d2inv);


// DerCross: computes cross derivative


float DerCross(float *p, int i, int s1, int s2, float dinv);

#endif
