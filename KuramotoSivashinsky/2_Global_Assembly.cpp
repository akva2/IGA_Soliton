// Codes for performing global assembly of the matrices and vectors obtained
// from the isogeometric discretization.
#include "../Settings.h"
#include "../IGA_1D/IGA_1D.h"
#include "KuramotoSivashinsky.h"


//----------------------------------------------------------------------------//
// Program for discretizing the Kuramoto-Sivashinsky equation.                //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::IsogeometricDiscretization()
{
  CompressedRowStorage(polynomial,elements,&RowSize,KnotVector,ContinuityVector);
  Assembly(EquationType,polynomial,elements,KnotVector,ContinuityVector);
}


//----------------------------------------------------------------------------//
// Module for assembling the matrices and vectors.                            //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::Assembly(string Type, int p, int nx, double* KX, int* CX)
{
  // Initializing the vector and relevant matrices in sparse PETSc-format.
  VecCreateSeq(PETSC_COMM_SELF,MatrixSize,&U0);
  MatCreateSeqAIJ(PETSC_COMM_SELF,MatrixSize,MatrixSize,SpanNumber,RowSize,&M0);
  MatCreateSeqAIJ(PETSC_COMM_SELF,MatrixSize,MatrixSize,SpanNumber,RowSize,&K2);
  MatCreateSeqAIJ(PETSC_COMM_SELF,MatrixSize,MatrixSize,SpanNumber,RowSize,&K4);
  if (Type == "GKS")
    MatCreateSeqAIJ(PETSC_COMM_SELF,MatrixSize,MatrixSize,SpanNumber,RowSize,&C3);
  else if (Type == "NKS")
    MatCreateSeqAIJ(PETSC_COMM_SELF,MatrixSize,MatrixSize,SpanNumber,RowSize,&H3);

  // Initializing the quadrature data.
  DynamicMatrix(&GF,p+4,2);
  DynamicMatrix(&G0,p+1,2);
  DynamicMatrix(&G1,p,2);
  DynamicMatrix(&G2,p-1,2);
  QuadratureData(GF,p+4);
  QuadratureData(G0,p+1);
  QuadratureData(G1,p);
  QuadratureData(G2,p-1);
  if (Type == "NKS"){
    DynamicMatrix(&GI,p+4,2);
    DynamicMatrix(&GO,p+5,2);
    QuadratureData(GI,p+4);
    QuadratureData(GO,p+5);
  }

  // Making the vector and matrices ready for assembly
  VecAssemblyBegin(U0);
  MatAssemblyBegin(M0,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K2,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K4,MAT_FINAL_ASSEMBLY);
  if (Type == "GKS")
    MatAssemblyBegin(C3,MAT_FINAL_ASSEMBLY);
  if (Type == "GKS")
    MatAssemblyBegin(H3,MAT_FINAL_ASSEMBLY);

  // Starting the isogeometric assembly.
  for (int n = 0; n < nx; n++){
    // Extracting local knot vector in the x-direction.
    double kx[2*p+2];
    for (int i = 0; i < 2*p+2; i++)
      kx[i] = KX[n+CX[n]+i];

    // Transforming the quadrature data in the x-direction.
    double ax = (kx[p+1]-kx[p])/2;
    double bx = (kx[p+1]+kx[p])/2;
    double gf[p+1][2];
    double g0[p+1][2],g1[p][2],g2[p-1][2];
    double gi[p+4][2],go[p+5][2];
    QuadratureTransform(gf,GF,p+4,ax,bx);
    QuadratureTransform(g0,G0,p+1,ax,bx);
    QuadratureTransform(g1,G1,p,ax,bx);
    QuadratureTransform(g2,G2,p-1,ax,bx);
    if (Type == "NKS"){
      QuadratureTransform(gi,GI,p+4,ax,bx);
      QuadratureTransform(go,GO,p+5,ax,bx);
    }

    // Computing the element quantities.
    double FE[p+1];
    double** M0E; double** K2E; double** K4E;
    double** C3E; double** H3E;
    DynamicMatrix(&M0E,p+1,p+1);
    DynamicMatrix(&K2E,p+1,p+1);
    DynamicMatrix(&K4E,p+1,p+1);
    ZeroVector(FE,p+1);
    LinAss024(0,p,kx,p+1,g0,M0E);
    LinAss024(1,p,kx,p,g1,K2E);
    LinAss024(2,p,kx,p-1,g2,K4E);
    if (Type == "GKS"){
      DynamicMatrix(&C3E,p+1,p+1);
      LinAss13(3,p,kx,p-1,g2,C3E);
    }
    if (Type == "NKS"){
      DynamicMatrix(&H3E,p+1,p+1);
      Hilbert3(p,kx,p+4,gi,p+5,go,H3E);
    }

    // Incrementing the element quantities into the global ones.
    int IX[p+1];
    ComputeElementIndex(IX,p,n+CX[n]);
    // MatSetValues(M0,p+1,const int idxm[],p+1,const int idxn[],M0E,ADD_VALUES);

    // Deleting the element matrices.
    DeleteMatrix(M0E,p+1);
    DeleteMatrix(K2E,p+1);
    DeleteMatrix(K4E,p+1);
    if (Type == "GKS")
      DeleteMatrix(C3E,p+1);
    else if (Type == "NKS")
      DeleteMatrix(H3E,p+1);
  }

  // Completing the vector and matrix assembly.
  VecAssemblyEnd(U0);
  MatAssemblyEnd(M0,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K2,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K4,MAT_FINAL_ASSEMBLY);
  if (Type == "GKS")
    MatAssemblyEnd(C3,MAT_FINAL_ASSEMBLY);
  else if (Type == "NKS")
    MatAssemblyEnd(H3,MAT_FINAL_ASSEMBLY);
}


//----------------------------------------------------------------------------//
// Creating the element contribution for the initial condition vector.        //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::InitAss(int p, double K[], int s, double g[][2], double f[])
{
  for (int k = 0; k <= s; k++){
    double E[p+1];
    Evaluate(E,K,p,g[k][0]);
    double f0 = 0; //!!!!
    double g0 = g[k][1];
    for (int i = 0; i <= p; i++){
      f[i] += g0*E[i]*f0;
    }
  }
}


//----------------------------------------------------------------------------//
// Creating the element contribution for the Hilbert transform of the third   //
// order convection matrix.                                                   //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::Hilbert3(int p, double K[], int si, double gi[][2],
  int so, double go[][2], double** m)
{
  for (int i = 0; i <= p; i++){
    for (int j = 0; j <= p; j++)
      m[i][j] = 0.0;
  }
  double R[p+1][si];
  for (int i = 0; i < si; i++){
    double EI[p+1];
    Derivative_2(EI,K,p,gi[i][0]);
    for (int j = 0; j <= p; j++){
      R[j][i] = EI[j];
    }
  }
  for (int k = 0; k < so; k++){
    double EO[p+1],I[p+1];
    Derivative_1(EO,K,p,go[k][0]);
    ZeroVector(I,p+1);
    double no = go[k][0];
    double wo = go[k][1];
    for (int l = 0; l < si; l++){
      double ni = gi[l][0];
      double wi = gi[l][1];
      for (int i = 0; i <= p; i++){
        I[i] += (wi*R[i][l])/(ni-no);
      }
    }
    for (int i = 0; i <= p; i++){
      for (int j = 0; j <= p; j++)
        m[i][j] += wo*EO[i]*I[j];
    }
  }
  for (int i = 0; i <= p; i++){
    for (int j = 0; j <= p; j++)
      m[i][j] /= pi;
  }
}
