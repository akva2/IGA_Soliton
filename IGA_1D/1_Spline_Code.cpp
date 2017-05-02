// Codes for generating knot vector, continuity vector and manipulating splines.
#include "../Settings.h"
#include "IGA_1D.h"


//----------------------------------------------------------------------------//
// Script for creating the knot vector.                                       //
//----------------------------------------------------------------------------//
double* IGA_1D::GenerateKnotVector(double a, double b, int n, int p, int c)
{
  double h = static_cast<double>((b-a)/n);
  int r = p-c;
  int L = 2*(p+1)+r*(n-1);
  double* K = new double[L];
  for (int i = 0; i <= p; i++){
    K[i] = a;
    K[L-1-i] = b;
  }
  for (int i = 1; i <= n-1; i++){
    double k = a+i*h;
    for (int j = p+(i-1)*r+1; j <= p+i*r; j++)
      K[j] = k;
  }
  return K;
}


//----------------------------------------------------------------------------//
// Script for creating the knot vector.                                       //
//----------------------------------------------------------------------------//
int* IGA_1D::GenerateContinuityVector(double* K, int n, int p, int c)
{
  int* C = new int[n];
  for (int i = 0; i < n; i++)
    C[i] = 0;
  int r = p-c;
  int it = 0;
  int d = 0;
  for (int i = 1; i <= n+(r-1)*(n-1); i++){
    double K0[2*p+2];
    for (int j = i; j < i+2*p+1; j++)
      K0[j-i] = K[j];
    if (K0[p] == K0[p+1])
      d++;
    else{
      it++;
      C[it] = C[it-1]+d;
      d = 0;
    }
  }
  return C;
}


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
      for (int i = 1; i < p-2; i++)
        EV[i+1] = C[i]*D[i]*E[i-1]-(C[i]+C[i+1])*D[i+1]*E[i]+C[i+1]*D[i+2]*E[i+1];
    }
    EV[p-1] = C[p-2]*D[p-2]*E[p-3]-(C[p-2]+C[p-1])*D[p-1]*E[p-2];
    EV[p] = C[p-2]*D[p-1]*E[p-2];
  }
  double s = static_cast<double>(p*(p-1));
  for (int i = 0; i < p+1; i++)
    EV[i] = s*EV[i];
}
