#ifndef _BOUNDARY
#define _BOUNDARY
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "map.h"

#define FRACABS 0.03125


// CreateSquareAbsorb: Creates cubic absortion zone around the domain


void CreateSquareAbsorb(int sx, int sy, int sz, 
			int nx, int ny, int nz,
			int bord, int absorb,
			float dx, float dy, float dz,
			float *fatAbsorb);


// CreateSphereAbsorb: Creates spherical absortion zone around the domain


void CreateSphereAbsorb(int sx, int sy, int sz, 
			int nx, int ny, int nz,
			int bord, int absorb,
			float dx, float dy, float dz,
			float *fatAbsorb);


// AbsorbingBoundary: Applies absortion zone on fields p and q


void AbsorbingBoundary(int sx, int sy, int sz,
		       float *fatAbsorb, float *p, float *q);


// RandomVelocityBoundary: creates a boundary with random velocity around domain


void RandomVelocityBoundary(int sx, int sy, int sz,
			    int nx, int ny, int nz,
			    int bord, int absorb,
			    float *vpz, float *vsv);

#endif
