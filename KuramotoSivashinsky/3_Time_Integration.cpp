// Codes for performing time-integration of the PDE.
#include "../Settings.h"
#include "../IGA_1D/IGA_1D.h"
#include "KuramotoSivashinsky.h"


//----------------------------------------------------------------------------//
// Module for time-integrating the Kuramoto-Sivashinsky equation with the     //
// semi-implicit backward differentiation formula (SBDF), up to order 6.      //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::SBDF_integrator(int n, int p, int ms, int T, int nt, int order)
{
  DynamicMatrix(&Ut,ms,nt+1);
}
