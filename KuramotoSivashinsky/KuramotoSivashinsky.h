// Derived class for all the functions and variables which will be used to solve
// the Kuramoto-Sivashinsky equation with isogeometric analysis.

#pragma once

class KuramotoSivashinsky : public IGA_1D
{
  ExprEval::Expression* expr;
  ExprEval::FunctionList* f;
  ExprEval::ValueList* v;
  std::vector<double*> arg;
public:
  // Functions for handling the XML-file.
  void InitializeClass(char* XMLFILE);
  void Extract_XML(char* XMLFILE);
  void Control_XML();
  void InitializeConstantParameters();

  // Functions for the isogeometric assembly.
  void IsogeometricDiscretization();
  void Assembly(string Type, int p, int nx, double* KX, int* CX);
  void InitAss(int p, double K[], int s, double g[][2], double f[]);
  void Hilbert3(int p, double K[], int si, double gi[][2],int so, double go[][2], double** m);

  // Functions for the time-integration.
  void SBDF_integrator(int n, int p, int ms, int T, int nt, int order);

  // Functions for computing the approximation error.
  void CalculateError(double U, int p, int k, int nx, double* KX, int* CX);
  void ComputeLocalError(double CO[], double kx[], double g0[][2], double g1[][2],
    double g2[][2], int p, double y0, double y1, double y2);

private:
  // Variables classifying the PDE.
  string EquationType;
  string BoundaryCondition;
  double alpha, beta, gamma;

  // Global matrices in PETSc-format.
  Vec U0;
  Mat M0; Mat K2; Mat K4;
  Mat C3; Mat H3;

  // Quadrature matrices.
  double** G0; double** G1; double** G2;
  double** GF; double** GC;
  double** GI; double** GO;

  // Parameters for the time-integration.
  double Time;
  int Order, TimeSteps;
  double EXP[6][6], BDF[6][7];
  double** Ut;

  // Parameters for post-processing.
  bool Error, Visualize;
  double Lebesgue,Sobolev_1,Sobolev_2;
  int nviz;

  double InitialCondition(double x, double t);
};
