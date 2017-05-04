// Codes for generating linear element matrices used in the assembly of global
// finite element matrices arising in the isogeometric discretization.
#include "../Settings.h"
#include "IGA_1D.h"


//----------------------------------------------------------------------------//
// Initializing three vectors for declaring the matrices' sparsity pattern,   //
// using the format Compressed Row Storage, in one spatial dimension.         //
//----------------------------------------------------------------------------//
void IGA_1D::CompressedRowStorage(int p, int n, double** VP, int** CP, int** RP,
  double* KX, int* CX)
{
  // Finding the number of nonzero elements and the matrix size.
  int s = n*pow((p+1),2);
  int ms = p+1;
  for (int i = 1; i < n; i++){
    int r = CX[i]-CX[i-1]+1;
    s -= pow(p+1-r,2);
    ms += r;
  }

  // Initializing the value-pointer vector.
  *VP = new double[s];
  for (int i = 0; i < s; i++)
    (*VP)[i] = 0.0;

  // Creating two auxiliary vectors for the elements' position.
  int COL[ms], ROW[ms];
  int pr = 0, pc = p, it = -1;
  for (int i = 0; i <= p; i++)
    ROW[i] = pr;
  for (int i = 1; i < n; i++){
    int r = CX[i]-CX[i-1]+1;
    pr += r;
    for (int j = 0; j < r; j++){
      it++;
      ROW[it+p+1] = pr;
      COL[it] = pc;
    }
    pc += r;
  }
  for (int i = ms-p-1; i < ms; i++)
    COL[i] = pc;

  // Initializing the column-pointer vector.
  *CP = new int[s];
  it = -1;
  for (int i = 0; i < ms; i++){
    for (int j = ROW[i]; j <= COL[i]; j++){
      it++;
      (*CP)[it] = j;
    }
  }

  // Initializing the row-pointer vector.
  *RP = new int[s];
  (*RP)[0] = 0;
  for (int i = 1; i < ms; i++){
    (*RP)[i] = (*RP)[i-1]+COL[i-1]-ROW[i-1]+1;
  }
}


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
  for (int k = 0; k < s; k++){
    double E[p+1];
    if (TYPE == 0)
      Evaluate(E,K,p,g[k][0]);
    else if (TYPE == 2)
      Derivative_1(E,K,p,g[k][0]);
    else if (TYPE == 4)
      Derivative_2(E,K,p,g[k][0]);
    double w = g[k][1];
    for (int i = 0; i <= p; i++)
      m[i][i] += w*pow(E[i],2);
    for (int i = 0; i <= p-1; i++){
      for (int j = 1; j <= p; j++){
        double d = w*E[i]*E[j];
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
  for (int k = 0; k < s; k++){
    double E1[p+1],E2[p+1];
    if (TYPE == 1){
      Evaluate(E1,K,p,g[k][0]);
      Derivative_1(E2,K,p,g[k][0]);
    }
    else if (TYPE == 3){
      Derivative_1(E1,K,p,g[k][0]);
      Derivative_2(E2,K,p,g[k][0]);
    }
    double w = g[k][1];
    for (int i = 0; i <= p; i++){
      for (int j = 0; j <= p; j++){
        m[i][j] = w*E1[i]*E2[j];
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
