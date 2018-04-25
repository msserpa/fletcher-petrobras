#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "timestep.h"
#include "boundary.h"
#include "source.h"
#include "utils.h"
#include "map.h"

enum Form {ISO, VTI, TTI};

#define MI 0.2           // stability factor to compute dt
#define ARGS 11          // tokens in executable command
#define dtOutput 0.01    // time spacing for section output

//#define _DUMP       // execution summary dump
//#undef  _DUMP     // execution summary dump
//#define _ABSOR_SQUARE    // use square absortion zone
#undef  _ABSOR_SQUARE  // don't use absortion zone
//#define _ABSOR_SPHERE    // use sphere absortion zone
#undef  _ABSOR_SPHERE  // don't use sphereabsortion zone
#define _RANDOM_BDRY     // use random boundary
//#undef _RANDOM_BDRY     // don'tuse random boundary

//#define SIGMA 20.0     // value of sigma (infinity) on formula 7 of Fletcher's paper
//#define SIGMA  6.0     // value of sigma on formula 7 of Fletcher's paper
//#define SIGMA  1.5     // value of sigma on formula 7 of Fletcher's paper
#define SIGMA  0.75      // value of sigma on formula 7 of Fletcher's paper
#define MAX_SIGMA 10.0   // above this value, SIGMA is considered infinite; as so, vsz=0

int main(int argc, char** argv) {

  enum Form prob;        // problem formulation
  int nx;                // grid points in x
  int ny;                // grid points in y
  int nz;                // grid points in z
  int bord=4;            // border size to apply the stencil at grid extremes
  int absorb;            // absortion zone size
  int sx;                // grid dimension in x (grid points + 2*border + 2*absortion)
  int sy;                // grid dimension in y (grid points + 2*border + 2*absortion)
  int sz;                // grid dimension in z (grid points + 2*border + 2*absortion)
  int st;                // number of time steps
  float dx;              // grid step in x
  float dy;              // grid step in y
  float dz;              // grid step in z
  float dt;              // time advance at each time step
  float tmax;            // desired simulation final time
  int ixSource;          // source x index
  int iySource;          // source y index
  int izSource;          // source z index
  int iSource;           // source index (ix,iy,iz) maped into 1D array
//PPL  int i, ix, iy, iz, it; // for indices
  int i, it;             // for indices
//PPL  char fNameAbs[128];    // prefix of absortion file
  char fNameSec[128];    // prefix of sections files

  it = 0; //PPL
    
  // input problem definition
  
  if (argc<ARGS) {
    printf("program requires %d input arguments; execution halted\n",ARGS-1);
    exit(-1);
  } 
  strcpy(fNameSec,argv[1]);
  nx=atoi(argv[2]);
  ny=atoi(argv[3]);
  nz=atoi(argv[4]);
  absorb=atoi(argv[5]);
  dx=atof(argv[6]);
  dy=atof(argv[7]);
  dz=atof(argv[8]);
  dt=atof(argv[9]);
  tmax=atof(argv[10]);

  // verify problem formulation

  if (strcmp(fNameSec,"ISO")==0) {
    prob=ISO;
  }
  else if (strcmp(fNameSec,"VTI")==0) {
    prob=VTI;
  }
  else if (strcmp(fNameSec,"TTI")==0) {
    prob=TTI;
  }
  else {
    printf("Input problem formulation (%s) is unknown\n", fNameSec);
    exit(-1);
  }

#ifdef _DUMP
  printf("Problem is ");
  switch (prob) {
  case ISO:
    printf("isotropic\n");
    break;
  case VTI:
    printf("anisotropic with vertical transversely isotropy using sigma=%f\n", SIGMA);
    break;
  case TTI:
    printf("anisotropic with tilted transversely isotropy using sigma=%f\n", SIGMA);
    break;
  }
#endif

  // grid dimensions from problem size

  sx=nx+2*bord+2*absorb;
  sy=ny+2*bord+2*absorb;
  sz=nz+2*bord+2*absorb;

  // number of time iterations

  st=ceil(tmax/dt);

  // source position

  ixSource=sx/2;
  iySource=sy/2;
  izSource=sz/2;
  iSource=ind(ixSource,iySource,izSource);

  // dump problem input data

#ifdef _DUMP
  printf("Grid size is (%d,%d,%d) with spacing (%.2f,%.2f,%.2f); simulated area (%.2f,%.2f,%.2f) \n", 
	 nx, ny, nz, dx, dy, dz, (nx-1)*dx, (ny-1)*dy, (nz-1)*dz);
  printf("Grid is extended by %d absortion points and %d border points at each extreme\n", absorb, bord);
  printf("Wave is propagated at internal+absortion points of size (%d,%d,%d)\n",
	 nx+2*absorb, ny+2*absorb, nz+2*absorb);
  printf("Source at coordinates (%d,%d,%d)\n", ixSource,iySource,izSource);
  printf("Will run %d time steps of %f to reach time %f\n", st, dt, st*dt);
#ifdef _OPENMP
#pragma omp parallel
  if (omp_get_thread_num() == 0) 
    printf("Execution with %d OpenMP threads\n", omp_get_num_threads());
#elif
  printf("Sequential execution\n");
#endif

#endif

  // allocate input anisotropy arrays
  
  float *vpz=NULL;      // p wave speed normal to the simetry plane
  vpz = (float *) malloc(sx*sy*sz*sizeof(float));

  float *vsv=NULL;      // sv wave speed normal to the simetry plane
  vsv = (float *) malloc(sx*sy*sz*sizeof(float));
  
  float *epsilon=NULL;  // Thomsen isotropic parameter
  epsilon = (float *) malloc(sx*sy*sz*sizeof(float));
  
  float *delta=NULL;    // Thomsen isotropic parameter
  delta = (float *) malloc(sx*sy*sz*sizeof(float));
  
  float *phi=NULL;     // isotropy simetry azimuth angle
  phi = (float *) malloc(sx*sy*sz*sizeof(float));
  
  float *theta=NULL;  // isotropy simetry deep angle
  theta = (float *) malloc(sx*sy*sz*sizeof(float));

  // input anisotropy arrays for selected problem formulation

  switch(prob) {

  case ISO:

    for (i=0; i<sx*sy*sz; i++) {
      vpz[i]=3000.0;
      epsilon[i]=0.0;
      delta[i]=0.0;
      phi[i]=0.0;
      theta[i]=0.0;
      vsv[i]=0.0;
    }
    break;

  case VTI:

    if (SIGMA > MAX_SIGMA) {
      printf("Since sigma (%f) is greater that threshold (%f), sigma is considered infinity and vsv is set to zero\n", 
		      SIGMA, MAX_SIGMA);
    }
    for (i=0; i<sx*sy*sz; i++) {
      vpz[i]=3000.0;
      epsilon[i]=0.24;
      delta[i]=0.1;
      phi[i]=0.0;
      theta[i]=0.0;
      if (SIGMA > MAX_SIGMA) {
	vsv[i]=0.0;
      } else {
	vsv[i]=vpz[i]*sqrtf(fabsf(epsilon[i]-delta[i])/SIGMA);
      }
    }
    break;

  case TTI:

    if (SIGMA > MAX_SIGMA) {
      printf("Since sigma (%f) is greater that threshold (%f), sigma is considered infinity and vsv is set to zero\n", 
		      SIGMA, MAX_SIGMA);
    }
    for (i=0; i<sx*sy*sz; i++) {
      vpz[i]=3000.0;
      epsilon[i]=0.24;
      delta[i]=0.1;
      phi[i]=0.0;
      theta[i]=atanf(1.0);
      if (SIGMA > MAX_SIGMA) {
	vsv[i]=0.0;
      } else {
	vsv[i]=vpz[i]*sqrtf(fabsf(epsilon[i]-delta[i])/SIGMA);
      }
    }
  } // end switch


  // stability condition
  
  float maxvel;
  maxvel=vpz[0]*sqrt(1.0+2*epsilon[0]);
  for (i=1; i<sx*sy*sz; i++) {
    maxvel=fmaxf(maxvel,vpz[i]*sqrt(1.0+2*epsilon[i]));
  }
  float mindelta=dx;
  if (dy<mindelta)
    mindelta=dy;
  if (dz<mindelta)
    mindelta=dz;
  float recdt;
  recdt=(MI*mindelta)/maxvel;
#ifdef _DUMP
  printf("Recomended maximum time step is %f; used time step is %f\n", recdt, dt);
#endif

  // random boundary speed

#ifdef _RANDOM_BDRY
  RandomVelocityBoundary(sx, sy, sz,
			 nx, ny, nz,
			 bord, absorb,
			 vpz, vsv);
#endif

#ifdef _DUMP
  DumpFieldToFile(sx, sy, sz,
		  0, sx-1,
		  0, sy-1,
		  0, sz-1,
		  dx, dy, dz,
		  vpz, "Velocity");
  
#endif

  // coeficients of derivatives at H1 operator

  float *ch1dxx=NULL;  // isotropy simetry deep angle
  ch1dxx = (float *) malloc(sx*sy*sz*sizeof(float));
  float *ch1dyy=NULL;  // isotropy simetry deep angle
  ch1dyy = (float *) malloc(sx*sy*sz*sizeof(float));
  float *ch1dzz=NULL;  // isotropy simetry deep angle
  ch1dzz = (float *) malloc(sx*sy*sz*sizeof(float));
  float *ch1dxy=NULL;  // isotropy simetry deep angle
  ch1dxy = (float *) malloc(sx*sy*sz*sizeof(float));
  float *ch1dyz=NULL;  // isotropy simetry deep angle
  ch1dyz = (float *) malloc(sx*sy*sz*sizeof(float));
  float *ch1dxz=NULL;  // isotropy simetry deep angle
  ch1dxz = (float *) malloc(sx*sy*sz*sizeof(float));
  float sinTheta, cosTheta, sin2Theta, sinPhi, cosPhi, sin2Phi;
  for (i=0; i<sx*sy*sz; i++) {
    sinTheta=sin(theta[i]);
    cosTheta=cos(theta[i]);
    sin2Theta=sin(2.0*theta[i]);
    sinPhi=sin(phi[i]);
    cosPhi=cos(phi[i]);
    sin2Phi=sin(2.0*phi[i]);
    ch1dxx[i]=sinTheta*sinTheta * cosPhi*cosPhi;
    ch1dyy[i]=sinTheta*sinTheta * sinPhi*sinPhi;
    ch1dzz[i]=cosTheta*cosTheta;
    ch1dxy[i]=sinTheta*sinTheta * sin2Phi;
    ch1dyz[i]=sin2Theta         * sinPhi;
    ch1dxz[i]=sin2Theta         * cosPhi;
  }
#ifdef _DUMP
  printf("ch1dxx[0]=%f; ch1dyy[0]=%f; ch1dzz[0]=%f; ch1dxy[0]=%f; ch1dxz[0]=%f; ch1dyz[0]=%f\n",
    	 ch1dxx[0], ch1dyy[0], ch1dzz[0], ch1dxy[0], ch1dxz[0], ch1dyz[0]);
#endif

  // coeficients of H1 and H2 at PDEs

  float *v2px=NULL;  // coeficient of H2(p)
  v2px = (float *) malloc(sx*sy*sz*sizeof(float));
  float *v2pz=NULL;  // coeficient of H1(q)
  v2pz = (float *) malloc(sx*sy*sz*sizeof(float));
  float *v2sz=NULL;  // coeficient of H1(p-q) and H2(p-q)
  v2sz = (float *) malloc(sx*sy*sz*sizeof(float));
  float *v2pn=NULL;  // coeficient of H2(p)
  v2pn = (float *) malloc(sx*sy*sz*sizeof(float));
  for (i=0; i<sx*sy*sz; i++){
    v2sz[i]=vsv[i]*vsv[i];
    v2pz[i]=vpz[i]*vpz[i];
    v2px[i]=v2pz[i]*(1.0+2.0*epsilon[i]);
    v2pn[i]=v2pz[i]*(1.0+2.0*delta[i]);
  }
#ifdef _DUMP
  printf("v2sz[0]=%f; v2pz[0]=%f; v2px[0]=%f; v2pn[0]=%f\n",
  	 v2sz[0], v2pz[0], v2px[0], v2pn[0]);
#endif

  // pressure fields at previous, current and future time steps
  
  float *pback=NULL;
  pback = (float *) malloc(sx*sy*sz*sizeof(float)); 

  float *outt=NULL;
  outt = (float *) malloc(sx*sy*sz*sizeof(float)); 

  float *pp=NULL;
  pp = (float *) malloc(sx*sy*sz*sizeof(float)); 
  float *pc=NULL;
  pc = (float *) malloc(sx*sy*sz*sizeof(float)); 
  float *qp=NULL;
  qp = (float *) malloc(sx*sy*sz*sizeof(float)); 
  float *qc=NULL;
  qc = (float *) malloc(sx*sy*sz*sizeof(float)); 
  for (i=0; i<sx*sy*sz; i++) {
    pp[i]=0.0f; pc[i]=0.0f; 
    qp[i]=0.0f; qc[i]=0.0f;
    pback[i]=0.0f; outt[i]=0.0f;
  }
  InsertSource(dt,0,iSource,pc,qc);
  InsertSource(dt,0,iSource,pp,qp);

  // absortion zone

  float *fatAbsorb=NULL;
#ifdef _ABSOR_SQUARE
  fatAbsorb = (float *) malloc(sx*sy*sz*sizeof(float));
  CreateSquareAbsorb(sx, sy, sz, 
		     nx, ny, nz,
		     bord, absorb,
		     dx, dy, dz,
		     fatAbsorb);
#endif
#ifdef _ABSOR_SPHERE
  fatAbsorb = (float *) malloc(sx*sy*sz*sizeof(float));
  CreateSphereAbsorb(sx, sy, sz, 
		     nx, ny, nz,
		     bord, absorb,
		     dx, dy, dz,
		     fatAbsorb);
#endif

#ifdef _DUMP
  if (fatAbsorb != NULL) {
    DumpFieldToFile(sx, sy, sz,
		    0, sx-1,
		    0, sy-1,
		    0, sz-1,
		    dx, dy, dz,
		    fatAbsorb, "Absorb");
  }
#endif

  double f_start, f_stop, f_total = 0.0;
  double f_io_start, f_io_stop, f_io_total = 0.0;
  double b_start, b_stop, b_total = 0.0;
  double b_io_start, b_io_stop, b_io_total = 0.0;

  // slices

//PPL  char fName[10];
  int ixStart=0;
  int ixEnd=sx-1;
  int iyStart=0;
  int iyEnd=sy-1;
  int izStart=0;
  int izEnd=sz-1;

  SlicePtr sPtr;
  sPtr=OpenSliceFile(ixStart, ixEnd,
		     iyStart, iyEnd,
		     izStart, izEnd,
		     dx, dy, dz, dt,
		     fNameSec);

  float tSim=0.0;
  int nOut=1;
  float tOut=nOut*dtOutput;
  f_io_start = omp_get_wtime();
  DumpSliceFile(sx,sy,sz,pc,sPtr);
  f_io_stop = omp_get_wtime();
  f_io_total += (f_io_stop - f_io_start);
    fprintf(stderr, "step %3d writted %.2lf MB in %8.5lf sec\n", 0, (sPtr->izEnd - sPtr->izStart + 1) * (sPtr->iyEnd - sPtr->iyStart + 1) * (sPtr->ixEnd-sPtr->ixStart+1) * sizeof(float) / 1024.0 / 1024.0, f_io_stop - f_io_start);
#ifdef _DUMP
  DumpSlicePtr(sPtr);
  DumpSliceSummary(sx,sy,sz,sPtr,dt,it,pc);
#endif
  
  // time advance

#ifndef ACC_MANAGED

#pragma acc data copyin(ch1dxx[0:sx*sy*sz], ch1dyy[0:sx*sy*sz], ch1dzz[0:sy*sy*sz],	\
			ch1dxy[0:sx*sy*sz], ch1dyz[0:sx*sy*sz], ch1dxz[0:sx*sy*sz], \
			v2px[0:sx*sy*sz], v2pz[0:sx*sy*sz], v2sz[0:sx*sy*sz], \
			v2pn[0:sx*sy*sz], pc[0:sz*sy*sz], qc[0:sz*sy*sz], \
			pp[0:sx*sy*sz], qp[0:sx*sy*sz], fatAbsorb[sx*sy*sz])

#endif
  for (it=1; it<=st; it++) {
    f_start = omp_get_wtime();

    Propagate(sx, sy, sz, bord,
	      dx, dy, dz, dt, it,
	      ch1dxx, ch1dyy, ch1dzz, 
	      ch1dxy, ch1dyz, ch1dxz, 
	      v2px, v2pz, v2sz, v2pn,
	      pp, pc, qp, qc);

    TimeForward(&pp, &pc, &qp, &qc);

    InsertSourceTimestep(dt,it,iSource,pc,qc);

#if ((defined _ABSOR_SPHERE) || (defined _ABSOR_SQUARE))
    AbsorbingBoundary(sx, sy, sz, fatAbsorb, pc, qc);
#endif

    f_stop = omp_get_wtime();
    f_total += (f_stop - f_start);

    tSim=it*dt;
    if (tSim >= tOut) {
#ifndef ACC_MANAGED

#pragma acc update host(pc[0:sx*sy*sz])

#endif
      f_io_start = omp_get_wtime();
      DumpSliceFile(sx,sy,sz,pc,sPtr);
      f_io_stop = omp_get_wtime();
      f_io_total += (f_io_stop - f_io_start);

      fprintf(stderr, "step %3d writted %.2lf MB in %8.5lf sec\n", it, (sPtr->izEnd - sPtr->izStart + 1) * (sPtr->iyEnd - sPtr->iyStart + 1) * (sPtr->ixEnd-sPtr->ixStart+1) * sizeof(float) / 1024.0 / 1024.0, f_io_stop - f_io_start);

      tOut=(++nOut)*dtOutput;
#ifdef _DUMP
      DumpSliceSummary(sx,sy,sz,sPtr,dt,it,pc);
#endif
    }
  }

  fprintf(stderr, "\nForward took %8.5lf sec\n", f_total);
  fprintf(stderr, "Write took  %8.5lf sec\n\n", f_io_total);

CloseSliceFile(sPtr);

  /* benchmark */
  SlicePtr sPtr2;
  sPtr2=OpenSliceFile2(ixStart, ixEnd,
         iyStart, iyEnd,
         izStart, izEnd,
         dx, dy, dz, dt,
         fNameSec);

  tSim=0.0;
  nOut=1;
  tOut=nOut*dtOutput;

    #ifndef ACC_MANAGED

    #pragma acc data copyin(ch1dxx[0:sx*sy*sz], ch1dyy[0:sx*sy*sz], ch1dzz[0:sy*sy*sz], \
          ch1dxy[0:sx*sy*sz], ch1dyz[0:sx*sy*sz], ch1dxz[0:sx*sy*sz], \
          v2px[0:sx*sy*sz], v2pz[0:sx*sy*sz], v2sz[0:sx*sy*sz], \
          v2pn[0:sx*sy*sz], pc[0:sz*sy*sz], qc[0:sz*sy*sz], \
          pp[0:sx*sy*sz], pback[0:sx*sy*sz], qp[0:sx*sy*sz], fatAbsorb[sx*sy*sz]) copyout(outt[0:sx*sy*sz])

    #endif
    for (it=1; it<=st; it++) {
    b_start = omp_get_wtime();

    Propagate(sx, sy, sz, bord,
        dx, dy, dz, dt, it,
        ch1dxx, ch1dyy, ch1dzz, 
        ch1dxy, ch1dyz, ch1dxz, 
        v2px, v2pz, v2sz, v2pn,
        pp, pc, qp, qc);

    TimeForward(&pp, &pc, &qp, &qc);

    InsertSourceTimestep(dt,it,iSource,pc,qc);

#if ((defined _ABSOR_SPHERE) || (defined _ABSOR_SQUARE))
    AbsorbingBoundary(sx, sy, sz, fatAbsorb, pc, qc);
#endif

    b_stop = omp_get_wtime();
    b_total += (b_stop - b_start);

    tSim=it*dt;
    if (tSim >= tOut) {
#ifndef ACC_MANAGED

#endif
      b_io_start = omp_get_wtime();
      DumpSliceFile2(sx,sy,sz,pback,sPtr2);
      b_io_stop = omp_get_wtime();
      b_io_total += (b_io_stop - b_io_start);

      fprintf(stderr, "step %3d read %.2lf MB in %8.5lf sec\n", it, (sPtr2->izEnd - sPtr2->izStart + 1) * (sPtr2->iyEnd - sPtr2->iyStart + 1) * (sPtr2->ixEnd-sPtr2->ixStart+1) * sizeof(float) / 1024.0 / 1024.0, b_io_stop - b_io_start);

      #pragma acc update device(pback[0:sx*sy*sz])

      Compare(sx, sy, sz, bord,
        pp, pback, outt);

      tOut=(++nOut)*dtOutput;
#ifdef _DUMP
      DumpSliceSummary(sx,sy,sz,sPtr,dt,it,pc);
#endif
    }
  }

  fprintf(stderr, "\nBackward took %8.5lf sec\n", b_total);
  fprintf(stderr, "Read took  %8.5lf sec\n\n", b_io_total);

CloseSliceFile2(sPtr2);

  exit(0);    
}
  