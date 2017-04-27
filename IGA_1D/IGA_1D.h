// Parent class for all the functions which can execute algorithms related to
// isogeometric analysis in one spatial dimension. 

class IGA_1D
{
public:
  void Evaluate(double E[], double K[], int p, double t);
  void Derivative_1(double EV[], double K[], int p, double t);
  void Derivative_2(double EV[], double K[], int p, double t);
  void ZeroVector(double V[], int length);

protected:
  int polynomial;
  int continuity;
  int elements;
};
