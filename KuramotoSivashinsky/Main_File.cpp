#include "../Settings.h"

int main(int argc, char* argv[])
{
  // Initializing the decimal precision globally.
  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(5);

  IGA_1D CCC;

  int p = atoi(argv[1]);
  double t = atof(argv[2]);
  double E[p+1];
  double K[2*p+2] = {0,0,0,0,1,1,1,1}; //

  CCC.Derivative_2(E,K,p,t);
  for (int i = 0; i <= p; i++)
    cout << E[i] << " ";
  cout << endl;
  return 0;
}
