// Codes for generating quadrature data of any order based on the orthogonal
// Gauss-Legendre polynomials.
#include "../Settings.h"
#include "IGA_1D.h"


//----------------------------------------------------------------------------//
// Program for finding the quadrature nodes and weights based on orthogonal   //
// Gauss-Legendre polynomial of order n.                                      //
//----------------------------------------------------------------------------//
void IGA_1D::QuadratureData(double** G, int n)
{
  int s = floor(n/2);
  for (int i = 1; i <= s; i++){
    double x0 = cos((i*pi)/(static_cast<double>(n+0.5)));
    double x = Olver(n,x0);
    double w = 2/((1-x*x)*pow(Legendre_1(n,x),2));
    G[n-i][0] = x;
    G[i-1][0] = -x;
    G[n-i][1] = w;
    G[i-1][1] = w;
  }
  if (n%2 == 1){
    G[s][0] = 0;
    G[s][1] = 2/(pow(Legendre_1(n,0),2));
  }
}


//----------------------------------------------------------------------------//
// Olver's method for finding the root of an equation, with a third order     //
// convergence rate.                                                          //
//----------------------------------------------------------------------------//
double IGA_1D::Olver(int n, double x)
{
  double x0 = x;
  double tol = pow(10,-14);
  double err = 0.0;
  do {
    double d0 = Legendre_0(n,x0);
    double d1 = Legendre_1(n,x0);
    double d2 = (2*x0*d1-n*(n+1)*d0)/(1-pow(x,2));
    double x1 = x0-d0/d1-0.5*(d2/d1)*pow(d0/d1,2);
    err = fabs(x1-x0);
    x0 = x1;
  } while (err > tol);
  return x0;
}


//----------------------------------------------------------------------------//
// Evaluation of the n-th order Gauss-Legendre polynomial at a given point x, //
// using dynamic programming (memoization).                                   //
//----------------------------------------------------------------------------//
double IGA_1D::Legendre_0(int n, double x)
{
  double Y;
  if (n == 0)
    Y = 1;
  else if (n == 1)
    Y = x;
  else if (n >= 2){
    double T[3] = {0,1,x};
    for (int i = 1; i <= n-1; i++){
      T[0] = T[1];
      T[1] = T[2];
      T[2] = ((2*i+1)*x*T[1]-i*T[0])/(i+1);
    }
    Y = T[2];
  }
  return Y;
}


//----------------------------------------------------------------------------//
// Evaluation of the first derivative of the n-th order Gauss-Legendre        //
// polynomial at a given point x.                                             //
//----------------------------------------------------------------------------//
double IGA_1D::Legendre_1(int n, double x)
{
  double Y;
  if (n == 0)
    Y = 0;
  else if (n >= 1){
    Y = ((pow(n,2)+n)/(2*n+1))*(Legendre_0(n-1,x)-Legendre_0(n+1,x));
    Y = Y/(1-pow(x,2));
  }
  return Y;
}


//----------------------------------------------------------------------------//
// Transforming quadrature data from reference element to physical element.   //
//----------------------------------------------------------------------------//
void IGA_1D::QuadratureTransform(double g[][2], double** G, int s, double a, double b)
{
  for (int i = 0; i < s; i++){
    g[i][0] = a*G[i][0]+b;
    g[i][1] = a*G[i][1];
  }
}
