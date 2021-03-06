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
  InitializeConstantParameters();
  KnotVector = GenerateKnotVector(StartPoint,EndPoint,elements,polynomial,continuity);
  ContinuityVector = GenerateContinuityVector(KnotVector,elements,polynomial,continuity);
  cout << InitialCondition(1,2) << std::endl;
}


//----------------------------------------------------------------------------//
// Module for extracting raw data from the XML-file.                          //
//----------------------------------------------------------------------------//
void KuramotoSivashinsky::Extract_XML(char* XMLFILE)
{
  // Opening the XML-file.
  TiXmlDocument doc(XMLFILE);
  doc.LoadFile();

  // Extracting all data from the XML-file.
  TiXmlElement* EL = doc.RootElement()->FirstChildElement();
  while (EL){
    // Looping through the KnotVector-tag.
    TiXmlElement* EL1 = EL->FirstChildElement();
    while (EL1){
      if (std::string(EL1->Value()) == "StartPoint")
        StartPoint = atof(EL1->FirstChild()->Value());
      if (std::string(EL1->Value()) == "EndPoint")
        EndPoint = atof(EL1->FirstChild()->Value());
      if (std::string(EL1->Value()) == "polynomial")
        polynomial = atoi(EL1->FirstChild()->Value());
      if (std::string(EL1->Value()) == "elements")
        elements = atoi(EL1->FirstChild()->Value());
      EL1 = EL1->NextSiblingElement();
    }

    // Looping through the PDE-tag.
    TiXmlElement* EL2 = EL->FirstChildElement();
    while (EL2){
      if (std::string(EL2->Value()) == "EquationType")
        EquationType = EL2->FirstChild()->Value();
      if (std::string(EL2->Value()) == "alpha")
        alpha = atof(EL2->FirstChild()->Value());
      if (std::string(EL2->Value()) == "beta")
        beta = atof(EL2->FirstChild()->Value());
      if (std::string(EL2->Value()) == "gamma")
        gamma = atof(EL2->FirstChild()->Value());
      if (std::string(EL2->Value()) == "BoundaryCondition")
        BoundaryCondition = EL2->FirstChild()->Value();
      EL2 = EL2->NextSiblingElement();
    }

    // Looping through the TimeIntegrator-tag.
    TiXmlElement* EL3 = EL->FirstChildElement();
    while (EL3){
      if (std::string(EL3->Value()) == "Time")
        Time = atof(EL3->FirstChild()->Value());
      if (std::string(EL3->Value()) == "Order")
        Order = atoi(EL3->FirstChild()->Value());
      if (std::string(EL3->Value()) == "TimeSteps")
        TimeSteps = atoi(EL3->FirstChild()->Value());
      if (std::string(EL3->Value()) == "InitialCondition")
      {
        expr = new ExprEval::Expression;
        f = new ExprEval::FunctionList;
        v = new ExprEval::ValueList;
        f->AddDefaultFunctions();
        v->AddDefaultValues();
        v->Add("x",0,false);
        v->Add("t",0,false);
        expr->SetFunctionList(f);
        expr->SetValueList(v);
        expr->Parse(EL3->FirstChild()->Value());
        arg.push_back(v->GetAddress("x"));
        arg.push_back(v->GetAddress("t"));
      }
      EL3 = EL3->NextSiblingElement();
    }

    // Looping through the PostProcessing-tag.
    TiXmlElement* EL4 = EL->FirstChildElement();
    while (EL4){
      if (std::string(EL4->Value()) == "Error"){
        int temp = atoi(EL4->FirstChild()->Value());
        Error = static_cast<bool>(temp);
      }
      if (std::string(EL4->Value()) == "Lebesgue")
        Lebesgue = atof(EL4->FirstChild()->Value());
      if (std::string(EL4->Value()) == "Sobolev_1")
        Sobolev_1 = atof(EL4->FirstChild()->Value());
      if (std::string(EL4->Value()) == "Sobolev_2")
        Sobolev_2 = atof(EL4->FirstChild()->Value());
      if (std::string(EL4->Value()) == "Visualize"){
        int temp = atoi(EL4->FirstChild()->Value());
        Visualize = static_cast<bool>(temp);
      }
      if (std::string(EL4->Value()) == "nviz")
        nviz = atoi(EL4->FirstChild()->Value());
      EL4 = EL4->NextSiblingElement();
    }

    // Stopping the while-loop.
    EL = EL->NextSiblingElement();
  }
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
  if (Error == 1){
    if ((Lebesgue > 0)&&(Sobolev_1 > 0)&&(Sobolev_2 > 0))
      post++;
  }
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


void KuramotoSivashinsky::InitializeConstantParameters()
{
  double A[6][6] = {1.0,0.0,0.0,0.0,0.0,0.0,
          2.0,-1.0,0.0,0.0,0.0,0.0,
          3.0,-3.0,1.0,0.0,0.0,0.0,
          4.0,-6.0,4.0,-1.0,0.0,0.0,
          5.0,-10.0,10.0,-5.0,1.0,0.0,
          6.0,-15.0,20.0,-15.0,6.0,-1.0};
  double B[6][7] = {1.0,-1.0,0.0,0.0,0.0,0.0,0.0,
          3.0/2,-2.0,1.0/2,0.0,0.0,0.0,0.0,
          11.0/6,-3.0,3.0/2,-1.0/3,0.0,0.0,0.0,
          25.0/12,-4.0,3.0,-4.0/3,1.0/4,0.0,0.0,
          137.0/60,-5.0,5.0,-10.0/3,5.0/4,-1.0/5,0.0,
          147.0/60,-6.0,15.0/2,-20.0/3,15.0/4,-6.0/5,1.0/6};
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 6; j++)
      EXP[i][j] = A[i][j];
    for (int j = 0; j < 7; j++)
      BDF[i][j] = B[i][j];
  }
}
