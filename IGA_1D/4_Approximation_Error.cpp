// Codes for calculating the Lebesgue and Sobolev norms of the final numerical
// solution, provided that the IVP/BVP has an exact analytical solution.
#include "../Settings.h"
#include "IGA_1D.h"


//----------------------------------------------------------------------------//
// Script for computing the Lebesgue L2-norm.                                 //
//----------------------------------------------------------------------------//
double IGA_1D::Lebesgue_L2(int p, double C[], double k[], double n, double w, double F)
{
  double E[p+1];
  ZeroVector(E,p+1);
  Evaluate(E,k,p,n);
  double FV = DotProduct(C,E,p+1);
  return w*pow((F-FV),2);
}


//----------------------------------------------------------------------------//
// Script for computing the Sobolev H1-seminorm.                              //
//----------------------------------------------------------------------------//
double IGA_1D::Sobolev_H1(int p, double C[], double k[], double n, double w, double F)
{
  double E[p+1];
  ZeroVector(E,p+1);
  Derivative_1(E,k,p,n);
  double FV = DotProduct(C,E,p+1);
  return w*pow((F-FV),2);
}


//----------------------------------------------------------------------------//
// Script for computing the Sobolev H2-seminorm.                              //
//----------------------------------------------------------------------------//
double IGA_1D::Sobolev_H2(int p, double C[], double k[], double n, double w, double F)
{
  double E[p+1];
  ZeroVector(E,p+1);
  Derivative_2(E,k,p,n);
  double FV = DotProduct(C,E,p+1);
  return w*pow((F-FV),2);
}


//----------------------------------------------------------------------------//
// Computing the dot product of two vectors.                                  //
//----------------------------------------------------------------------------//
double IGA_1D::DotProduct(double a[], double b[], int s)
{
  double d = 0.0;
  for (int i = 0; i < s; i++)
    d += a[i]*b[i];
}
