// Parent class for all the functions which can execute algorithms related to
// isogeometric analysis in one spatial dimension.

#pragma once

class IGA_1D
{
public:
  // Functions for the spline manipulation.
  double* GenerateKnotVector(double a, double b, int n, int p, int c);
  int* GenerateContinuityVector(double* K, int n, int p, int c);
  void Evaluate(double E[], double K[], int p, double t);
  void Derivative_1(double EV[], double K[], int p, double t);
  void Derivative_2(double EV[], double K[], int p, double t);

  // Functions for the quadrature initialization.
  void QuadratureData(double** G, int n);
  double Olver(int n, double x);
  double Legendre_0(int n, double x);
  double Legendre_1(int n, double x);
  void QuadratureTransform(double g[][2], double** G, int n, double a, double b);

  // Functions for the linear assembly.
  void LinAss024(int TYPE, int p, double K[], int s, double g[][2], double** m);
  void LinAss13(int TYPE, int p, double K[], int s, double g[][2], double** m);
  void ZeroVector(double v[], int length);
  void DynamicMatrix(double*** m, int row, int col);
  void DeleteMatrix(double** m, int n);

protected:
  int polynomial;
  int continuity;
  int elements;
  double* KnotVector;
  int* ContinuityVector;
};
