#include "../Settings.h"
#include "../IGA_1D/IGA_1D.h"
#include "KuramotoSivashinsky.h"

int main(int argc, char* argv[])
{
  // Initializing the decimal precision globally.
  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(15);

  char* XMLFILE = argv[1];

  KuramotoSivashinsky KS;
  PetscInitialize(&argc,&argv,(char*)0,"help");
  KS.InitializeClass(XMLFILE);
  KS.IsogeometricDiscretization();
  PetscFinalize();

  return 0;
}
