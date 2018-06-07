#include "boundary.h"


// CreateSquareAbsorb: Creates cubic absortion zone around the domain


void CreateSquareAbsorb(int sx, int sy, int sz, 
			int nx, int ny, int nz,
			int bord, int absorb,
			float dx, float dy, float dz,
			float *fatAbsorb) {
  float dumpFactor;
  int i, ix, iy, iz;
//PPL  int ixCenter=sx/2;
//PPL  int iyCenter=sy/2;
//PPL  int izCenter=sz/2;
  int distz, disty, distx, dist;
  
  dumpFactor=-FRACABS/(float)(absorb*absorb);
  printf("Square absortion exponent of %f implies absortion range from %8.6f to %8.6f\n",
	 -FRACABS, expf(dumpFactor), expf(dumpFactor*(float)absorb*(float)absorb));
  for (iz=0; iz<sz; iz++) {
    for (iy=0; iy<sy; iy++) {
      for (ix=0; ix<sx; ix++) {
	i=ind(ix,iy,iz);
	if ((iz>=bord+absorb && iz<sz-bord-absorb) &&
	    (iy>=bord+absorb && iy<sy-bord-absorb) &&
	    (ix>=bord+absorb && ix<sx-bord-absorb)) {
	  fatAbsorb[i]=1.0;
	}
	else if ((iz>=bord && iz<sz-bord) &&
		 (iy>=bord && iy<sy-bord) &&
		 (ix>=bord && ix<sx-bord)) {
	  distz=(iz>=bord+absorb+nz)?iz-bord-absorb-nz+1:bord+absorb-iz;
	  disty=(iy>=bord+absorb+ny)?iy-bord-absorb-ny+1:bord+absorb-iy;
	  distx=(ix>=bord+absorb+nx)?ix-bord-absorb-nx+1:bord+absorb-ix;
	  dist=(disty>distz)?disty:distz;
	  dist=(dist >distx)?dist :distx;
	  fatAbsorb[i]=expf(dumpFactor*(float)(dist*dist));
	}
	else
	  fatAbsorb[i]=1.0;
      }
    }
  }
}  


// CreateSphereAbsorb: Creates spherical absortion zone around the domain


void CreateSphereAbsorb(int sx, int sy, int sz, 
			int nx, int ny, int nz,
			int bord, int absorb,
			float dx, float dy, float dz,
			float *fatAbsorb) {
  float dumpFactor;
  int i, ix, iy, iz;
  int ixCenter=sx/2;
  int iyCenter=sy/2;
  int izCenter=sz/2;
  float distz, disty, distx;
  float radius, radiusMin, radiusMax, radius2, radiusMin2, radiusMax2;
  
  distz=dz*(izCenter-bord-absorb);
  disty=dy*(iyCenter-bord-absorb);
  distx=dx*(ixCenter-bord-absorb);
  radiusMin2=distz*distz+disty*disty+distx*distx;
  radiusMin=sqrt(radiusMin2);
  
  radiusMax=dz*(izCenter-bord);
  radiusMax2=radiusMax*radiusMax;

  printf("radiusMin=%f; radiusMax=%f\n",radiusMin, radiusMax);

  dumpFactor=-FRACABS/((radiusMax-radiusMin)*(radiusMax-radiusMin));
  for (iz=0; iz<sz; iz++) {
    distz=dz*(float)(iz-izCenter); distz=distz*distz;
    for (iy=0; iy<sy; iy++) {
      disty=dy*(float)(iy-iyCenter); disty=disty*disty;
      for (ix=0; ix<sx; ix++) {
	distx=dx*(float)(ix-ixCenter); distx=distx*distx;
	radius2=distz+disty+distx;
	i=ind(ix,iy,iz);
	if (radius2 <= radiusMin2) {
	  fatAbsorb[i]=1.0;
	}
	else if (radius2 <= radiusMax2) {
	  radius=sqrtf(radius2)-radiusMin;
 	  fatAbsorb[i]=expf(dumpFactor*radius*radius);
	}
	else
	  fatAbsorb[i]=0.0;
      }
    }
  }
}  


// AbsorbingBoundary: Applies absortion zone on fields p and q


void AbsorbingBoundary(int sx, int sy, int sz,
		       float *fatAbsorb, float *p, float *q) {
  int i;

  if (fatAbsorb == NULL) {
    printf("invoking **(AbsorbingBoundary)** with null fatAbsorb; programming error\n");
    exit(-1);
  }

#pragma omp parallel for 
  for (i=0; i<sx*sy*sz; i++) {
    p[i]*=fatAbsorb[i];
    q[i]*=fatAbsorb[i];
  }
}


// RandomVelocityBoundary: creates a boundary with random velocity around domain


void RandomVelocityBoundary(int sx, int sy, int sz,
			    int nx, int ny, int nz,
			    int bord, int absorb,
			    float *vpz, float *vsv) {

  int i, ix, iy, iz;
  int distx, disty, distz, dist;
  int ivelx, ively, ivelz;
  float bordDist;
  float frac, rfac;
  int firstIn, bordLen;
  float maxP, maxS;

  // maximum speed of P and S within bounds
  maxP=0.0; maxS=0.0;
  for (iz=bord+absorb; iz<nz+bord+absorb; iz++) {
    for (iy=bord+absorb; iy<ny+bord+absorb; iy++) {
      for (i=ind(bord+absorb,iy,iz); i<ind(nx+bord+absorb,iy,iz); i++) {
	maxP=fmaxf(vpz[i],maxP);
	maxS=fmaxf(vsv[i],maxS);
      }
    }
  }

  bordLen=bord+absorb-1;   // last index on low absortion zone
  firstIn=bordLen+1;       // first index inside input grid 
  frac=1.0/(float)(absorb);

  for (iz=0; iz<sz; iz++) {
    for (iy=0; iy<sy; iy++) {
      for (ix=0; ix<sx; ix++) {
	i=ind(ix,iy,iz);
	// do nothing inside input grid
	if ((iz>=firstIn && iz<=bordLen+nz) &&
	    (iy>=firstIn && iy<=bordLen+ny) &&
	    (ix>=firstIn && ix<=bordLen+nx)) {
	  continue;
	}
	// random speed inside absortion zone
	else if ((iz>=bord && iz<=2*bordLen+nz) &&
		 (iy>=bord && iy<=2*bordLen+ny) &&
		 (ix>=bord && ix<=2*bordLen+nx)) {
	  if (iz>bordLen+nz) {
	    distz=iz-bordLen-nz;
	    ivelz=bordLen+nz;
	  } else if (iz<firstIn) {
	    distz=firstIn-iz;
	    ivelz=firstIn;
	  } else {
	    distz=0;
	    ivelz=iz;
	  }
	  if (iy>bordLen+ny) {
	    disty=iy-bordLen-ny;
	    ively=bordLen+ny;
	  } else if (iy<firstIn) {
	    disty=firstIn-iy;
	    ively=firstIn;
	  } else {
	    disty=0;
	    ively=iy;
	  }
	  if (ix>bordLen+nx) {
	    distx=ix-bordLen-nx;
	    ivelx=bordLen+nx;
	  } else if (ix<firstIn) {
	    distx=firstIn-ix;
	    ivelx=firstIn;
	  } else {
	    distx=0;
	    ivelx=ix;
	  }
	  dist=(disty>distz)?disty:distz;
	  dist=(dist >distx)?dist :distx;
	  bordDist=(float)(dist)*frac;
	  rfac=(float)rand()/(float)RAND_MAX;
	  vpz[i]=vpz[ind(ivelx,ively,ivelz)]*(1.0-bordDist)+
	    maxP*rfac*bordDist;
	  vsv[i]=vsv[ind(ivelx,ively,ivelz)]*(1.0-bordDist)+
	    maxS*rfac*bordDist;
	}
	// null speed at border
	else
	{
	  vpz[i]=0.0;
	  vsv[i]=0.0;
	}
      }
    }
  }
}
