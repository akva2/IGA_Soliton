// Codes for generating linear element matrices used in the assembly of global
// finite element matrices arising in the isogeometric discretization.
#include "../Settings.h"
#include "IGA_1D.h"


//----------------------------------------------------------------------------//
// Creating the element contribution for the mass matrix, second order or     //
// fourth order stiffness matrices.                                           //
//----------------------------------------------------------------------------//
void IGA_1D::LinAss024(int TYPE, int p, double K[], int s, double g[][2], double** m)
{
  for (int i = 0; i <= p; i++){
    for (int j = 0; j <= p; j++)
      m[i][j] = 0.0;
  }
  for (int k = 0; k <= s; k++){
    double E[p+1];
    if (TYPE == 0)
      Evaluate(E,K,p,g[k][0]);
    else if (TYPE == 2)
      Derivative_1(E,K,p,g[k][0]);
    else if (TYPE == 4)
      Derivative_2(E,K,p,g[k][0]);
    double g0 = g[k][1];
    for (int i = 0; i <= p; i++)
      m[i][i] += g0*pow(E[i],2);
    for (int i = 0; i <= p-1; i++){
      for (int j = 1; j <= p; j++){
        double d = g0*E[i]*E[j];
        m[i][j] += d;
        m[j][i] += d;
      }
    }
  }
}


//----------------------------------------------------------------------------//
// Creating the element contribution for the first order or third order       //
// convection matrices.                                                       //
//----------------------------------------------------------------------------//
void IGA_1D::LinAss13(int TYPE, int p, double K[], int s, double g[][2], double** m)
{
  for (int i = 0; i <= p; i++){
    for (int j = 0; j <= p; j++)
      m[i][j] = 0.0;
  }
  for (int k = 0; k <= s; k++){
    double E1[p+1],E2[p+1];
    if (TYPE == 1){
      Evaluate(E1,K,p,g[k][0]);
      Derivative_1(E2,K,p,g[k][0]);
    }
    else if (TYPE == 3){
      Derivative_1(E1,K,p,g[k][0]);
      Derivative_2(E2,K,p,g[k][0]);
    }
    double g0 = g[k][1];
    for (int i = 0; i <= p; i++){
      for (int j = 0; j <= p; j++){
        m[i][j] = g0*E1[i]*E2[j];
      }
    }
  }
}


//----------------------------------------------------------------------------//
// Script for initializing a zero vector.                                     //
//----------------------------------------------------------------------------//
void IGA_1D::ZeroVector(double v[], int length)
{
  for (int i = 0; i < length; i++)
    v[i] = 0;
}


//----------------------------------------------------------------------------//
// Script for initializing a dynamic matrix.                                  //
//----------------------------------------------------------------------------//
void IGA_1D::DynamicMatrix(double*** m, int row, int col)
{
  *m = new double*[row];
  for (int i = 0; i < row; i++)
    (*m)[i] = new double[col];
}


//----------------------------------------------------------------------------//
// Script for deleting a dynamic matrix.                                      //
//----------------------------------------------------------------------------//
void IGA_1D::DeleteMatrix(double** m, int n)
{
  for (int i = 0; i < n; i++)
    delete[] m[i];
  delete[] m;
}
