// Codes for computing the error of the discretization.
#include "../Settings.h"
#include "../IGA_1D/IGA_1D.h"
#include "KuramotoSivashinsky.h"


//----------------------------------------------------------------------------//
// Program for calculating the approximation error, in the Sobolev seminorm.  //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::CalculateError(double U, int p, int k, int nx, double* KX, int* CX)
{
  double Y0 = 0.0, Y1 = 0.0, Y2 = 0.0;
  int c0 = 2*p+1;
  for (int n = 0; n < nx; n++){
    // Extracting local knot vector in the x-direction.
    double kx[2*p+2];
    for (int i = 0; i < 2*p+2; i++)
      kx[i] = KX[n+CX[n]+i];

    // Transforming the quadrature data in the x-direction.
    double ax = (kx[p+1]-kx[p])/2;
    double bx = (kx[p+1]+kx[p])/2;
    double g0[p+1][2],g1[p][2],g2[p-1][2];
    QuadratureTransform(g0,G0,p+1,ax,bx);
    QuadratureTransform(g1,G1,p,ax,bx);
    QuadratureTransform(g2,G2,p-1,ax,bx);

    // Extracting local coefficient vector in the x-direction.
    double CO[p+1];
    for (int i = 0; i <= p; i++){
      int in = (p-k)*n+1;
      CO[i] = in;
    }

    // Computing local error on element.
    double y0 = 0.0, y1 = 0.0, y2 = 0.0;
    ComputeLocalError(CO,kx,g0,g1,g2,p,y0,y1,y2);
    Y0 += y0;
    Y1 += y1;
    Y2 += y2;
  }

  // Finding norm of the exact errors.
  double L2h = sqrt(Y0/Lebesgue);
  double H1h = sqrt(Y1/Sobolev_1);
  double H2h = sqrt(Y2/Sobolev_2);

  // Storing the results in an .asc-file.
  ofstream FILE;
  FILE.setf(ios::fixed);
  FILE.setf(ios::showpoint);
  FILE.precision(15);
  FILE.open("Norms.asc");
  FILE << L2h << endl;
  FILE << H1h << endl;
  FILE << H2h << endl;
  FILE.close();
}


//----------------------------------------------------------------------------//
// Script for computing local error on an element.                            //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::ComputeLocalError(double CO[], double kx[], double g0[][2],
  double g1[][2], double g2[][2], int p, double y0, double y1, double y2)
{
  for (int i = 0; i < p+1; i++){
    // Computing the Lebesgue L2-norm.
    if (i < p+1){
      double F0 = 0; //!!!!
      double r = Lebesgue_L2(p,CO,kx,g0[i][0],g0[i][1],F0);
      y0 += r;
    }
    // Computing the Sobolev H1-seminorm.
    if (i < p){
      double F1 = 0; //!!!!
      double r = Sobolev_H1(p,CO,kx,g1[i][0],g1[i][1],F1);
      y1 += r;
    }
    // Computing the Sobolev H2-seminorm.
    if (i < p-1){
      double F2 = 0; //!!!!
      double r = Sobolev_H2(p,CO,kx,g2[i][0],g2[i][1],F2);
      y2 += r;
    }
  }
}
