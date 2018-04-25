#include "timestep.h"

#ifdef _DUMP
#include <omp.h>
#endif

#pragma acc routine(Der2) seq
#pragma acc routine(DerCross) seq

// Propagate: using Fletcher's equations, propagate waves one dt,
//            either forward or backward in time


void Propagate(int sx, int sy, int sz, int bord,
	       float dx, float dy, float dz, float dt, int it, 
	       float * restrict ch1dxx, float * restrict ch1dyy, float * restrict ch1dzz, 
	       float * restrict ch1dxy, float * restrict ch1dyz, float * restrict ch1dxz, 
	       float * restrict v2px, float * restrict v2pz, float * restrict v2sz, float * restrict v2pn,
	       float * restrict pp, float * restrict pc, float * restrict qp, float * restrict qc) {
  int ix, iy, iz, i;
  int strideX=ind(1,0,0)-ind(0,0,0);
  int strideY=ind(0,1,0)-ind(0,0,0);
  int strideZ=ind(0,0,1)-ind(0,0,0);
  float pxx, pyy, pzz, pxy, pyz, pxz;
  float qxx, qyy, qzz, qxy, qyz, qxz;
  float cpxx, cpyy, cpzz, cpxy, cpyz, cpxz;
  float cqxx, cqyy, cqzz, cqxy, cqyz, cqxz;
  float h1p, h2p;
  float h1q, h2q;
  float h1pmq, h2pmq;
  float rhsp, rhsq;
  float dxxinv=1.0/(dx*dx);
  float dyyinv=1.0/(dy*dy);
  float dzzinv=1.0/(dz*dz);
  float dxyinv=1.0/(dx*dy);
  float dxzinv=1.0/(dx*dz);
  float dyzinv=1.0/(dy*dz);

#ifdef _DUMP
  double startTime, endTime;
#endif

  // solve both equations in all internal grid points, 
  // including absortion zone

#ifdef _DUMP
  startTime = omp_get_wtime();
#endif

#ifndef _OPENACC
  
#pragma omp parallel for private(ix, iy, iz, i,				\
			     pxx, pyy, pzz, pxy, pyz, pxz,		\
			     qxx, qyy, qzz, qxy, qyz, qxz,		\
			     cpxx, cpyy, cpzz, cpxy, cpyz, cpxz,	\
			     cqxx, cqyy, cqzz, cqxy, cqyz, cqxz,	\
			     h1p, h2p,					\
			     h1q, h2q,					\
			     h1pmq, h2pmq,				\
			     rhsp, rhsq)
#else

#ifndef ACC_MANAGED

#pragma acc kernels							\
  present(ch1dxx[0:sx*sy*sz], ch1dyy[0:sx*sy*sz], ch1dzz[0:sy*sy*sz],	\
	  ch1dxy[0:sx*sy*sz], ch1dyz[0:sx*sy*sz], ch1dxz[0:sx*sy*sz],	\
	  v2px[0:sx*sy*sz], v2pz[0:sx*sy*sz], v2sz[0:sx*sy*sz],		\
	  v2pn[0:sx*sy*sz], pc[0:sz*sy*sz], qc[0:sz*sy*sz],		\
	  pp[0:sx*sy*sz], qp[0:sx*sy*sz])

#else

#pragma acc kernels

#endif

#endif

#pragma acc loop independent
    for (iz=bord; iz<sz-bord; iz++) {
#pragma acc loop independent
      for (iy=bord; iy<sy-bord; iy++) {
#pragma acc loop independent
	for (ix=bord; ix<sx-bord; ix++) {
	  i=ind(ix,iy,iz);
	      
	  // p derivatives, H1(p) and H2(p)
	  
	  pxx= Der2(pc, i, strideX, dxxinv);
	  pyy= Der2(pc, i, strideY, dyyinv);
	  pzz= Der2(pc, i, strideZ, dzzinv);
	  pxy= DerCross(pc, i, strideX, strideY, dxyinv);
	  pyz= DerCross(pc, i, strideY, strideZ, dyzinv);
	  pxz= DerCross(pc, i, strideX, strideZ, dxzinv);
	  
	  cpxx=ch1dxx[i]*pxx;
	  cpyy=ch1dyy[i]*pyy;
	  cpzz=ch1dzz[i]*pzz;
	  cpxy=ch1dxy[i]*pxy;
	  cpxz=ch1dxz[i]*pxz;
	  cpyz=ch1dyz[i]*pyz;
	  h1p=cpxx+cpyy+cpzz+cpxy+cpxz+cpyz;
	  h2p=pxx+pyy+pzz-h1p;
	  
	  // q derivatives, H1(q) and H2(q)
	  
	  qxx= Der2(qc, i, strideX, dxxinv);
	  qyy= Der2(qc, i, strideY, dyyinv);
	  qzz= Der2(qc, i, strideZ, dzzinv);
	  qxy= DerCross(qc, i, strideX,  strideY, dxyinv);
	  qyz= DerCross(qc, i, strideY,  strideZ, dyzinv);
	  qxz= DerCross(qc, i, strideX,  strideZ, dxzinv);
	  
	  cqxx=ch1dxx[i]*qxx;
	  cqyy=ch1dyy[i]*qyy;
	  cqzz=ch1dzz[i]*qzz;
	  cqxy=ch1dxy[i]*qxy;
	  cqxz=ch1dxz[i]*qxz;
	  cqyz=ch1dyz[i]*qyz;
	  h1q=cqxx+cqyy+cqzz+cqxy+cqxz+cqyz;
	  h2q=qxx+qyy+qzz-h1q;
	  
	  
	  // p-q derivatives, H1(p-q) and H2(p-q)
	  
	  h1pmq=h1p-h1q;
	  h2pmq=h2p-h2q;
	  
	  // rhs of p and q equations
	  
	  rhsp=v2px[i]*h2p + v2pz[i]*h1q + v2sz[i]*h1pmq;
	  rhsq=v2pn[i]*h2p + v2pz[i]*h1q - v2sz[i]*h2pmq;
	  
	  // new p and q
	  
	  pp[i]=2.0*pc[i] - pp[i] + rhsp*dt*dt;
	  qp[i]=2.0*qc[i] - qp[i] + rhsq*dt*dt;

	}
      }
    } // end nested for loops
#ifdef _DUMP
  endTime = omp_get_wtime();
  printf("Timestep %d calculated in %8.6f secs.\n", it, endTime - startTime);
#endif
}


// TimeForward: swap array pointers on time forward array propagation


void TimeForward(float **pp, float **pc, float **qp, float **qc) {
  float *tmp;

  tmp=*pp;
  *pp=*pc;
  *pc=tmp;

  tmp=*qp;
  *qp=*qc;
  *qc=tmp;
}
