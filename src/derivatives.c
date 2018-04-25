#include "derivatives.h"


// Der2: computes second derivative

#pragma acc routine seq
float Der2(float *p, int i, int s, float d2inv){
  float res;
  res=(K0*p[i]+
       K1*(p[i+s]+p[i-s])+
       K2*(p[i+2*s]+p[i-2*s]) +
       K3*(p[i+3*s]+p[i-3*s]) +
       K4*(p[i+4*s]+p[i-4*s])
       )*(d2inv);
  return(res);
}


// DerCross: computes cross derivative


#pragma acc routine seq
float DerCross(float *p, int i, int s11, int s21, float dinv){
  float res;
  int s12=2*s11;
  int s13=3*s11;
  int s14=4*s11;
  int s22=2*s21;
  int s23=3*s21;
  int s24=4*s21;
  res=(L11*(p[i+s21+s11]-p[i+s21-s11]-p[i-s21+s11]+p[i-s21-s11])+
       L12*(p[i+s21+s12]-p[i+s21-s12]-p[i-s21+s12]+p[i-s21-s12]+p[i+s22+s11]-p[i+s22-s11]-p[i-s22+s11]+p[i-s22-s11])+
       L13*(p[i+s21+s13]-p[i+s21-s13]-p[i-s21+s13]+p[i-s21-s13]+p[i+s23+s11]-p[i+s23-s11]-p[i-s23+s11]+p[i-s23-s11])+
       L14*(p[i+s21+s14]-p[i+s21-s14]-p[i-s21+s14]+p[i-s21-s14]+p[i+s24+s11]-p[i+s24-s11]-p[i-s24+s11]+p[i-s24-s11])+
       L22*(p[i+s22+s12]-p[i+s22-s12]-p[i-s22+s12]+p[i-s22-s12])+
       L23*(p[i+s22+s13]-p[i+s22-s13]-p[i-s22+s13]+p[i-s22-s13]+p[i+s23+s12]-p[i+s23-s12]-p[i-s23+s12]+p[i-s23-s12])+
       L24*(p[i+s22+s14]-p[i+s22-s14]-p[i-s22+s14]+p[i-s22-s14]+p[i+s24+s12]-p[i+s24-s12]-p[i-s24+s12]+p[i-s24-s12])+
       L33*(p[i+s23+s13]-p[i+s23-s13]-p[i-s23+s13]+p[i-s23-s13])+
       L34*(p[i+s23+s14]-p[i+s23-s14]-p[i-s23+s14]+p[i-s23-s14]+p[i+s24+s13]-p[i+s24-s13]-p[i-s24+s13]+p[i-s24-s13])+
       L44*(p[i+s24+s14]-p[i+s24-s14]-p[i-s24+s14]+p[i-s24-s14]))*(dinv);
  return(res);
}

