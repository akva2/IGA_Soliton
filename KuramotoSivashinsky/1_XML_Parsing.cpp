// Codes for extracting XML-data and controlling that they are correct, before
// using them to intialize the class for the PDE.
#include "../Settings.h"
#include "../IGA_1D/IGA_1D.h"
#include "KuramotoSivashinsky.h"


//----------------------------------------------------------------------------//
// Program for initializing the class KuramotoSivashinsky, based on the data  //
// obtained from the XML-file.                                                //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::InitializeClass(char* XMLFILE)
{
  Extract_XML(XMLFILE);
  Control_XML();
  GenerateKnotVector(StartPoint,EndPoint,elements,polynomial,continuity);
  GenerateContinuityVector(KnotVector,elements,polynomial,continuity);
}


//----------------------------------------------------------------------------//
// Module for extracting raw data from the XML-file.                          //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::Extract_XML(char* XMLFILE)
{

}


//----------------------------------------------------------------------------//
// Module for controlling that the data from the XML-file are correct.        //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::Control_XML()
{
  // Initializing the global control variable.
  int control = 0;

  // Test 1: Checking the equation type and corresponding constants.
  if (EquationType == "KS"){
    if ((alpha == 1)&&(beta == 0)&&(gamma > 0))
      control++;
  }
  if (EquationType == "GKS"){
    if ((alpha > 0)&&(beta > 0)&&(gamma > 0))
      control++;
  }
  if (EquationType == "NKS"){
    if ((alpha == 0)&&(beta > 0)&&(gamma == 1))
      control++;
  }

  // Test 2: Controlling the chosen type of boundary condition.
  if (BoundaryCondition == "Dirichlet")
    control++;
  else if (BoundaryCondition == "Periodic")
    control++;

  // Test 3: Controlling the spline parameters.
  if ((polynomial >= 2)&&(elements >= 2)){
    continuity = polynomial-1;
    if (StartPoint < EndPoint)
      control++;
  }

  // Test 4: Controlling the time-integrator parameters.
  if (Time > 0){
    if ((Order >= 1)&&(Order <= 6)){
      if (TimeSteps >= Order)
        control++;
    }
  }

  // Test 5: Controlling the parameters for post-processing.
  int post = 0;
  if (Error == 1)
    post++;
  if (Visualize == 1){
    if (nviz >= 2)
      post++;
  }
  if (post > 0)
    control++;

  // Stopping the program if the XML-file is incorrect.
  if (control < 5){
    cout << "Invalid XML-input" << endl;
    return;
  }
  else
    cout << "Valid XML-input" << endl;
}
