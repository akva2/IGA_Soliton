// Derived class for all the functions and variables which will be used to solve
// the Kuramoto-Sivashinsky equation with isogeometric analysis.

#pragma once

class KuramotoSivashinsky : public IGA_1D
{
public:
  // Functions for handling the XML-file.
  void InitializeClass(char* XMLFILE);
  void Extract_XML(char* XMLFILE);
  void Control_XML();

  // Functions for the isogeometric assembly.
  void IsogeometricDiscretization();
  void Assembly(string Type, int p, int nx, double* KX, int* CX);
  void InitAss(int p, double K[], int s, double g[][2], double f[]);
  void Hilbert3(int p, double K[], int si, double gi[][2],int so, double go[][2], double** m);

  void CalculateError(double U, int p, int k, int nx, double* KX, int* CX);
  void ComputeLocalError(double CO[], double kx[], double g0[][2], double g1[][2],
    double g2[][2], int p, double y0, double y1, double y2);

private:
  // Variables classifying the PDE.
  string EquationType;
  string BoundaryCondition;
  double alpha, beta, gamma;

  // Parameters for the time-integration.
  double Time;
  int Order, TimeSteps;

  // Quadrature matrices.
  double** G0; double** G1; double** G2;
  double** GF; double** GC;
  double** GI; double** GO;

  // Parameters for post-processing.
  bool Error, Visualize;
  double Lebesgue,Sobolev_1,Sobolev_2;
  int nviz;
};
