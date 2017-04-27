// Codes for generating knot vector, continuity vector and manipulating splines.
#include "../Settings.h"


//----------------------------------------------------------------------------//
// Script for evaluating a spline.                                            //
//----------------------------------------------------------------------------//
void IGA_1D::Evaluate(double EV[], double K[], int p, double t)
{
  double R[2*p+1];
  ZeroVector(R,2*p+1);
  R[p] = 1.0;
  for (int k = 0; k < p; k++){
    for (int i = 0; i < 2*p-k; i++){
      double LBF = R[i];
      double a = K[i+k+1]-K[i];
      double LRU = 0.0;
      if (a > 0)
        LRU = (t-K[i])/a;
      double RBF = R[i+1];
      a = K[i+k+2]-K[i+1];
      double RRD = 0;
      if (a > 0)
        RRD = (K[i+k+2]-t)/a;
      R[i] = LBF*LRU+RBF*RRD;
    }
  }
  for (int i = 0; i <= p; i++)
    EV[i] = R[i];
}


//----------------------------------------------------------------------------//
// Script for evaluating the first derivative of a spline.                    //
//----------------------------------------------------------------------------//
void IGA_1D::Derivative_1(double EV[], double K[], int p, double t)
{
  double C[p],E[p-1];
  ZeroVector(C,p);
  ZeroVector(E,p-1);
  for (int i = 0; i < p; i++){
    double r = (K[i+p+1]-K[i+1])/p;
    if (r != 0)
      C[i] = 1/r;
  }
  double K0[2*p];
  for (int i = 0; i < 2*p; i++)
    K0[i] = K[i+1];
  Evaluate(E,K0,p-1,t);
  EV[0] = -C[0]*E[0];
  for (int i = 0; i < p-1; i++)
    EV[i+1] = C[i]*E[i]-C[i+1]*E[i+1];
  EV[p] = C[p-1]*E[p-1];
}


//----------------------------------------------------------------------------//
// Script for evaluating the second derivative of a spline.                   //
//----------------------------------------------------------------------------//
void IGA_1D::Derivative_2(double EV[], double K[], int p, double t)
{
  double C[p],D[p+1],E[p-1];
  ZeroVector(C,p);
  ZeroVector(E,p+1);
  for (int i = 0; i < p+1; i++){
    if (i < p){
      double r = K[i+p+1]-K[i+1];
      if (r != 0)
        C[i] = 1/r;
    }
    double r = K[i+p]-K[i+1];
      if (r != 0)
        D[i] = 1/r;
  }
  if (p == 2){
    EV[0] = C[0]*D[1];
    EV[1] = -(C[0]+C[1])*D[1];
    EV[2] = C[1]*D[1];
  }
  else if (p >= 3){
    double K0[2*p];
    for (int i = 0; i < 2*p-2; i++)
      K0[i] = K[i+2];
    Evaluate(E,K0,p-2,t);
    EV[0] = C[0]*D[1]*E[0];
    EV[1] = -(C[0]+C[1])*D[1]*E[0]+C[1]*D[2]*E[1];
    if (p >= 4){
      for (int i = 2; i < p-2; i++)
        EV[i] = C[i]*D[i]*E[i-2];
    }
    EV[p-1] = C[p-2]*D[p-2]*E[p-3]-(C[p-2]+C[p-1])*D[p-1]*E[p-2];
    EV[p] = C[p-2]*D[p-1]*E[p-2];
  }
  double s = static_cast<double>(p*(p-1));
  for (int i = 0; i < p+1; i++)
    EV[i] = s*EV[i];
}


//----------------------------------------------------------------------------//
// Script for initializing a zero vector.                                     //
//----------------------------------------------------------------------------//
void IGA_1D::ZeroVector(double V[], int length)
{
  for (int i = 0; i < length; i++)
    V[i] = 0;
}
