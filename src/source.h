#ifndef _SOURCE
#define _SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// values used on source computation

#define FCUT        40.0
#define PICUBE      31.00627668029982017537
#define TWOSQRTPI    3.54490770181103205458
#define THREESQRTPI  5.31736155271654808184


// Source: compute source value at time it*dt


float Source(float dt, int it);


// InsertSource: compute and insert source value at index iSource of arrays p and q


void InsertSource(float dt, int it, int iSource, 
		  float * restrict p, float * restrict q);
void InsertSourceTimestep(float dt, int it, int iSource, 
			  float * restrict p, float * restrict q);

#endif
