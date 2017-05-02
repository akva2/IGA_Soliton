// Derived class for all the functions and variables which will be used to solve
// the Kuramoto-Sivashinsky equation with isogeometric analysis.

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

private:
  string EquationType;
  string BoundaryCondition;
  double StartPoint, EndPoint;
  double alpha, beta, gamma;
  double Time;
  int Order, TimeSteps;

  // Quadrature matrices.
  double** G0; double** G1; double** G2;
  double** GF; double** GC;
  double** GI; double** GO;

  // Parameters for post-processing.
  bool Error, Visualize;
  int nviz;
};
