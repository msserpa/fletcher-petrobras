#include "source.h"


// Source: compute source value at time it*dt

#pragma acc routine seq
float Source(float dt, int it){

  float tf, fc, fct, expo;
  tf=TWOSQRTPI/FCUT;
  fc=FCUT/THREESQRTPI;
  fct=fc*(((float)it)*dt-tf);
  expo=PICUBE*fct*fct;
  return ((1.0f-2.0f*expo)*expf(-expo));
}


// InsertSource: compute and insert source value at index iSource of arrays p and q


void InsertSource(float dt, int it, int iSource, 
		  float * restrict p, float * restrict q) {
  float src;

  {
    src=Source(dt,it);
    p[iSource]+=src;
    q[iSource]+=src;
  }
}

// InsertSourceTimestep: created only because OpenACC pragma reference two arrays
// that is not present on target when InsertSource is called outside of timestep

void InsertSourceTimestep(float dt, int it, int iSource, 
		     float * restrict p, float * restrict q) {
  float src;

#pragma acc kernels present(p, q)
    //PPL laço inserido abaixo para que a cláusula present funcione
    // descoberto na base da tentativa e erro
    for (int i = 0; i < 1; i++) {
    src=Source(dt,it);
    p[iSource]+=src;
    q[iSource]+=src;
  }
}
