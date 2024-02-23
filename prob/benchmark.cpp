// -------------------------------------------------------------------------
// benchmark.cpp - implementation of the benchmark class.
// -------------------------------------------------------------------------
// Copyright (c) 2021 LMCV/UFC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------
// Created:      21-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     14-Mar-2013    Evandro Parente Junior
//               Evaluate input parameter changed to int*.
//
// Modified:     18-Mar-2013    Evandro Parente Junior
//               Created Booth and Rastrigin problems.
//
// Modified:     04-Dec-2017    Marina Alves Maia
//               Created CONSTR, TNK, ZDT1, ZDT6, SCH and KUR multiobjective
//               problems.
//
// Modified:     18-Oct-2018    Leonardo Gonçalves Ribeiro
//               Created PeaksC and PeaksSAO problem.
//
// Modified:     07-Jul-2019    Marina Alves Maia
//               Created BraninSAO, GoldSteinSAO, ColvilleSAO, TridSAO,
//               RastriginSAO, Hart3SAO and Hart6SAO problems.
//
// -------------------------------------------------------------------------

#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "benchmark.h"
#include "gbldef.h"
#include "gblvar.h"
#include "input.h"
#include "utl.h"
#include "vec.h"
#include "mat.h"
#include "matvec.h"
#include "sysmat.h"
#include "group.h"
#include "optsolution.h"

using namespace std;


// -------------------------------------------------------------------------
// Static Variables:
//
int cBenchmark :: InpNumVar = 0;

// -------------------------------------------------------------------------
// Set read functions labels:
//
static const bool ReadFuncRegister[] =
{
  ProbMap( ).Insert("PROBLEM.NUMBER.VARIABLES",cBenchmark::ReadInpNumVar),
};
// -------------------------------------------------------------------------
// Register problems on the problem factory:
//
static const bool registeredProb[] =
{
  cProblemFactory :: Register("Peaks"                   , MakeProb<cPeaksC,cPeaksD>),
  cProblemFactory :: Register("Branin"                  , MakeProb<cBraninC,cBraninD>),
  cProblemFactory :: Register("Hartmann3"               , MakeProb<cHart3C,cHart3D>),
  cProblemFactory :: Register("Hartmann6"               , MakeProb<cHart6C,cHart6D>),
  cProblemFactory :: Register("Rastrigin"               , MakeProb<cRastriginC,cRastriginD>),
  cProblemFactory :: Register("ConstrainedBranin"       , MakeProb<cConstrainedBraninC,cConstrainedBraninD>),
  cProblemFactory :: Register("Kitayama5"               , MakeProb<cKit5C,cKit5D>),
  cProblemFactory :: Register("ThreeBarTruss"           , MakeProb<c3BarTrussC,c3BarTrussD>),
  cProblemFactory :: Register("ThreeBarTrussFAST"       , MakeProb<c3BarTrussCFAST>),
  cProblemFactory :: Register("ThreeBarTrussABAQUS"     , MakeProb<c3BarTrussCABAQUS>),
  cProblemFactory :: Register("ThreeBarTrussDIANA"      , MakeProb<c3BarTrussCDIANA>),
  cProblemFactory :: Register("TenBarTrussFAST"         , MakeProb<c10BarTrussFAST>),
  cProblemFactory :: Register("TenBarTrussABAQUS"       , MakeProb<c10BarTrussABAQUS>),
  cProblemFactory :: Register("TenBarTrussFrequencyFAST", MakeProb<c10BarTrussFrequencyFAST>),
  cProblemFactory :: Register("NowackiBeam"             , MakeProb<cNowackiBeamC,cNowackiBeamD>),
  cProblemFactory :: Register("Beam"                    , MakeProb<cBeamC,cBeamD>),
  cProblemFactory :: Register("CONSTR"                  , MakeProb<cCONSTRC>),
  cProblemFactory :: Register("TNK"                     , MakeProb<cTNKC>),
  cProblemFactory :: Register("ZDT6"                    , MakeProb<cZDT6C>),
  cProblemFactory :: Register("ZDT1"                    , MakeProb<cZDT1C>),
  cProblemFactory :: Register("SCH"                     , MakeProb<cSCHC>),
  cProblemFactory :: Register("KUR"                     , MakeProb<cKURC>)
};

// -------------------------------------------------------------------------
// Class cBenchmark:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cReadInpNumVar ==============================

void cBenchmark :: ReadInpNumVar(void)
{
  if (!(in >> InpNumVar))
  {
    cout << "Error in the input of the number of optimizations." << endl;
    exit(0);
  }
}

// =========================== cBenchmark ==================================

cBenchmark :: cBenchmark(void)
{
}

// ========================== cBenchmark ===================================

cBenchmark :: ~cBenchmark(void)
{
}

// -------------------------------------------------------------------------
// Class cBenchDiscrete:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cBenchDiscrete ==============================

cBenchDiscrete :: cBenchDiscrete(void)
{
}

// ========================== cBenchContinuous =============================

cBenchContinuous :: cBenchContinuous(void)
{
}

// =========================== ~cBenchDiscrete =============================

cBenchDiscrete :: ~cBenchDiscrete(void)
{
  delete []ListDim;
  delete []List;
}

// ========================== ~cBenchContinuous ============================

cBenchContinuous :: ~cBenchContinuous(void)
{
  delete[] Low;
  delete[] Upp;
}

// ============================== PrintVar =================================

void cBenchDiscrete :: PrintVar(int *algvar)
{
  for (int i = 0; i < NumVar; i++)
    cout << 1+i << " " << List[i][algvar[i]] << endl;

  cout << endl;
}

// ============================== PrintVar =================================

void cBenchDiscrete :: DecodeVar(int *algvar, cVector &xsamp)
{
  xsamp.Resize(NumVar);

  for (int i = 0; i < NumVar; i++)
      xsamp[i] = (List[i][algvar[i]]-List[i].Min())/(List[i].Max()-List[i].Min());
}

// ============================= WriteVar ==================================

void cBenchDiscrete :: WriteVar(int *algvar, ostream &out)
{
  for (int i = 0; i < NumVar; i++)
    out << List[i][algvar[i]] << "  ";

  out << endl;
}

// ============================= GetBounds =================================

void cBenchDiscrete :: GetBounds(int i, int *low, int *upp)
{
  *low = 0;
  *upp = ListDim[i]-1;
}

// ============================= GetBounds =================================

void cBenchDiscrete :: GetVarBounds(int i, double &low, double &upp)
{
  low = List[i].Min();
  upp = List[i].Max();
}

// ========================== GetBoundsDouble ==============================

void cBenchContinuous :: GetDblBounds(double *low,double *upp)
{
  for (int i = 0; i < NumVar; i++)
  {
    low[i] = Low[i];
    upp[i] = Upp[i];
  }
}

// -------------------------------------------------------------------------
// Class cPeaksC
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPeaksC ===============================

cPeaksC :: cPeaksC(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = -3;
    Upp[i] = 3;
  }
}

// ============================ Evaluate ==============================

void cPeaksC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  float a1 = 3*(1-x[0])*(1-x[0])*exp(-(x[0]*x[0])-(x[1]+1)*(x[1]+1));
  float a2 = -10*(x[0]/5 - x[0]*x[0]*x[0] - x[1]*x[1]*x[1]*x[1]*x[1])*exp(-(x[0]*x[0])-(x[1])*(x[1]));
  float a3 = -exp(-(x[0]+1)*(x[0]+1)-(x[1])*(x[1]));
  float fobj = a1+a2+a3/3;

  fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cPeaksD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPeaksD ==================================

cPeaksD :: cPeaksD(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-2;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(6.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = -3.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================= Evaluate ================================

void cPeaksD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Objective function evaluation.

  double a1 = 3*(1-x[0])*(1-x[0])*exp(-(x[0]*x[0])-(x[1]+1)*(x[1]+1));
  double a2 = -10*(x[0]/5 - x[0]*x[0]*x[0] - x[1]*x[1]*x[1]*x[1]*x[1])*exp(-(x[0]*x[0])-(x[1])*(x[1]));
  double a3 = -exp(-((x[0]+1)*(x[0]+1))-(x[1])*(x[1]));
  double fobj = a1+a2+a3/3;

  fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cBraninC
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBraninC ===============================

cBraninC :: cBraninC(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = -5;
  Upp[0] = 10;
  Low[1] = 0;
  Upp[1] = 15;
}

// ============================ Evaluate ==============================

void cBraninC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5.0/PI;
    double r = 6.0;
    double s = 10.0;
    double t = 1.0/(8.0*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cBraninD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBraninD ===============================

cBraninD :: cBraninD(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(15.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  List[0].Resize(ListDim[0]);
  List[0][0] = -5.0;
  for (int j = 1; j < ListDim[0]; j++) List[0][j] = List[0][j-1] + step;

  List[1].Resize(ListDim[1]);
  List[1][0] = 0.0;
  for (int j = 1; j < ListDim[1]; j++) List[1][j] = List[1][j-1] + step;
}

// ============================ Evaluate ==============================

void cBraninD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5/PI;
    double r = 6;
    double s = 10;
    double t = 1/(8*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    fobjs[0] = fobj;
}

// ============================= cHart3C ===============================

cHart3C :: cHart3C(void)
{
  NumVar = 3;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i < NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cHart3C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,3);
    cMatrix P(4,3);

    A[0][0] = 3.0;
    A[0][1] = 10.0;
    A[0][2] = 30.0;

    A[1][0] = 0.1;
    A[1][1] = 10.0;
    A[1][2] = 35.0;

    A[2][0] = 3.0;
    A[2][1] = 10.0;
    A[2][2] = 30.0;

    A[3][0] = 0.10;
    A[3][1] = 10.0;
    A[3][2] = 35.0;

    P[0][0] = pow(10,-4)*3689;
    P[0][1] = pow(10,-4)*1170;
    P[0][2] = pow(10,-4)*2673;

    P[1][0] =pow(10,-4)*4699;
    P[1][1] =pow(10,-4)*4387;
    P[1][2] =pow(10,-4)*7470;

    P[2][0] =pow(10,-4)*1091;
    P[2][1] =pow(10,-4)*8732;
    P[2][2] =pow(10,-4)*5547;

    P[3][0] =pow(10,-4)*381;
    P[3][1] =pow(10,-4)*5743;
    P[3][2] =pow(10,-4)*8828;

    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 3; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = - outer;
}

// -------------------------------------------------------------------------
// Class cHart3D
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cHart3D ===============================

cHart3D :: cHart3D(void)
{
  NumVar = 3;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(1.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = 0.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================ Evaluate ==============================

void cHart3D :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,3);
    cMatrix P(4,3);

    A[0][0] = 3.0;
    A[0][1] = 10.0;
    A[0][2] = 30.0;

    A[1][0] = 0.1;
    A[1][1] = 10.0;
    A[1][2] = 35.0;

    A[2][0] = 3.0;
    A[2][1] = 10.0;
    A[2][2] = 30.0;

    A[3][0] = 0.10;
    A[3][1] = 10.0;
    A[3][2] = 35.0;

    P[0][0] = pow(10,-4)*3689;
    P[0][1] = pow(10,-4)*1170;
    P[0][2] = pow(10,-4)*2673;

    P[1][0] =pow(10,-4)*4699;
    P[1][1] =pow(10,-4)*4387;
    P[1][2] =pow(10,-4)*7470;

    P[2][0] =pow(10,-4)*1091;
    P[2][1] =pow(10,-4)*8732;
    P[2][2] =pow(10,-4)*5547;

    P[3][0] =pow(10,-4)*381;
    P[3][1] =pow(10,-4)*5743;
    P[3][2] =pow(10,-4)*8828;

    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 3; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = - outer;
}

// ============================= cHart6C ===============================

cHart6C :: cHart6C(void)
{
  NumVar = 6;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i < NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cHart6C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,6);
    cMatrix P(4,6);

    A[0][0] = 10;
    A[0][1] = 3;
    A[0][2] = 17;
    A[0][3] = 3.5;
    A[0][4] = 1.7;
    A[0][5] = 8;

    A[1][0] = 0.05;
    A[1][1] = 10;
    A[1][2] = 17;
    A[1][3] = 0.1;
    A[1][4] = 8;
    A[1][5] = 14;

    A[2][0] = 3;
    A[2][1] = 3.5;
    A[2][2] = 1.7;
    A[2][3] = 10;
    A[2][4] = 17;
    A[2][5] = 8;

    A[3][0] = 17;
    A[3][1] = 8;
    A[3][2] = 0.05;
    A[3][3] = 10;
    A[3][4] = 0.1;
    A[3][5] = 14;

    P[0][0] = pow(10,-4)*1312;
    P[0][1] = pow(10,-4)*1696;
    P[0][2] = pow(10,-4)*5569;
    P[0][3] = pow(10,-4)*124;
    P[0][4] = pow(10,-4)*8283;
    P[0][5] = pow(10,-4)*5886;

    P[1][0] =pow(10,-4)*2329;
    P[1][1] =pow(10,-4)*4135;
    P[1][2] =pow(10,-4)*8307;
    P[1][3] =pow(10,-4)*3736;
    P[1][4] =pow(10,-4)*1004;
    P[1][5] =pow(10,-4)*9991;

    P[2][0] =pow(10,-4)*2348;
    P[2][1] =pow(10,-4)*1451;
    P[2][2] =pow(10,-4)*3522;
    P[2][3] =pow(10,-4)*2883;
    P[2][4] =pow(10,-4)*3047;
    P[2][5] =pow(10,-4)*6650;

    P[3][0] =pow(10,-4)*4047;
    P[3][1] =pow(10,-4)*8828;
    P[3][2] =pow(10,-4)*8732;
    P[3][3] =pow(10,-4)*5743;
    P[3][4] =pow(10,-4)*1091;
    P[3][5] =pow(10,-4)*381;


    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 6; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = -outer;
}

// -------------------------------------------------------------------------
// Class cHart6D
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cHart6D ===============================

cHart6D :: cHart6D(void)
{
  NumVar = 6;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(1.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = 0.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================ Evaluate ==============================

void cHart6D :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,6);
    cMatrix P(4,6);

    A[0][0] = 10;
    A[0][1] = 3;
    A[0][2] = 17;
    A[0][3] = 3.5;
    A[0][4] = 1.7;
    A[0][5] = 8;

    A[1][0] = 0.05;
    A[1][1] = 10;
    A[1][2] = 17;
    A[1][3] = 0.1;
    A[1][4] = 8;
    A[1][5] = 14;

    A[2][0] = 3;
    A[2][1] = 3.5;
    A[2][2] = 1.7;
    A[2][3] = 10;
    A[2][4] = 17;
    A[2][5] = 8;

    A[3][0] = 17;
    A[3][1] = 8;
    A[3][2] = 0.05;
    A[3][3] = 10;
    A[3][4] = 0.1;
    A[3][5] = 14;

    P[0][0] = pow(10,-4)*1312;
    P[0][1] = pow(10,-4)*1696;
    P[0][2] = pow(10,-4)*5569;
    P[0][3] = pow(10,-4)*124;
    P[0][4] = pow(10,-4)*8283;
    P[0][5] = pow(10,-4)*5886;

    P[1][0] =pow(10,-4)*2329;
    P[1][1] =pow(10,-4)*4135;
    P[1][2] =pow(10,-4)*8307;
    P[1][3] =pow(10,-4)*3736;
    P[1][4] =pow(10,-4)*1004;
    P[1][5] =pow(10,-4)*9991;

    P[2][0] =pow(10,-4)*2348;
    P[2][1] =pow(10,-4)*1451;
    P[2][2] =pow(10,-4)*3522;
    P[2][3] =pow(10,-4)*2883;
    P[2][4] =pow(10,-4)*3047;
    P[2][5] =pow(10,-4)*6650;

    P[3][0] =pow(10,-4)*4047;
    P[3][1] =pow(10,-4)*8828;
    P[3][2] =pow(10,-4)*8732;
    P[3][3] =pow(10,-4)*5743;
    P[3][4] =pow(10,-4)*1091;
    P[3][5] =pow(10,-4)*381;


    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 6; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = -outer;
}

// -------------------------------------------------------------------------
// Class cRastriginC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cRastriginC ================================

cRastriginC :: cRastriginC(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  if (InpNumVar) NumVar = InpNumVar;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for(int i = 0; i < NumVar; i++)
  {
    Low[i] = -5.12;
    Upp[i] = 5.12;
  }
}

// ============================= Evaluate ================================

void cRastriginC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double f = 0;

  for(int i=0; i < NumVar; i++)
    f += x[i]*x[i] - 10*cos(2*PI*x[i]);

  f += 10*NumVar;

  fobjs[0] = f;
}

// -------------------------------------------------------------------------
// Class cRastriginD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cRastriginD ===============================

cRastriginD :: cRastriginD(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(2*5.12/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = -5.12;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================= Evaluate ================================

void cRastriginD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Objective function evaluation.

  double f = 0;

  for(int i=0; i < NumVar; i++)
    f += x[i]*x[i] - 10*cos(2*PI*x[i]);

  f += 10*NumVar;

  fobjs[0] = f;
}

// ============================= cBraninC ===============================

cConstrainedBraninC :: cConstrainedBraninC(void)
{
  NumVar = 2;
  NumConstr = 1;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = -5;
  Upp[0] = 10;
  Low[1] = 0;
  Upp[1] = 15;
}

// ============================ Evaluate ==============================

void cConstrainedBraninC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5.0/PI;
    double r = 6.0;
    double s = 10.0;
    double t = 1.0/(8.0*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    c[0] = 1 - (x[0] + 5)*x[1]/45;
    fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cConstrainedBraninC :: EvalExactConstraint(int index, cVector& x, double &c)
{
  // Single constraint evaluation.
    if (index == 0){
        c = 1 - (x[0] + 5)*x[1]/45;
    }
    else{
        cout << "Invalid index in EvalExactConstraint!";
        exit(0);
    }
}

// -------------------------------------------------------------------------
// Class cConstrainedBraninD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================= cConstrainedBraninD ============================

cConstrainedBraninD :: cConstrainedBraninD(void)
{
  NumVar = 2;
  NumConstr = 1;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(15.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  List[0].Resize(ListDim[0]);
  List[0][0] = -5.0;
  for (int j = 1; j < ListDim[0]; j++) List[0][j] = List[0][j-1] + step;

  List[1].Resize(ListDim[1]);
  List[1][0] = 0.0;
  for (int j = 1; j < ListDim[1]; j++) List[1][j] = List[1][j-1] + step;
}

// ============================ Evaluate ==============================

void cConstrainedBraninD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5/PI;
    double r = 6;
    double s = 10;
    double t = 1/(8*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    c[0] = 1 - (x[0] + 5)*x[1]/45;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cKit5C
// -------------------------------------------------------------------------

// ============================= cKit5C ===============================

cKit5C :: cKit5C(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cKit5C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
  c[0] = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  c[1] = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  c[2] = (10*x[0] + x[1])/7 - 1;

  fobjs[0] = fobj;
}


// ============================ Evaluate ==============================

void cKit5C :: EvalExactFobj(cVector &x, double &fobj)
{
  fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
}

// ============================ Evaluate ==============================

void cKit5C :: EvalExactConstraint(int index, cVector& x, double &c)
{
  // Single constraint evaluation.
  if (index == 0){
      c = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  }
  else if (index == 1){
      c = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  }
  else if (index == 2){
      c = (10*x[0] + x[1])/7 - 1;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cKit5C :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
  approxc[1] = 1;
  approxc[2] = 1;
}

// -------------------------------------------------------------------------
// Class cKit5D
// -------------------------------------------------------------------------

// ============================= cKit5D ===============================

cKit5D :: cKit5D(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-2;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(1.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = 0.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================ Evaluate ==============================

void cKit5D :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];
  // Objective function evaluation.

  double fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
  c[0] = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  c[1] = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  c[2] = (10*x[0] + x[1])/7 - 1;

  fobjs[0] = fobj;
}


// ============================ Evaluate ==============================

void cKit5D :: EvalExactFobj(int *algvar, double &fobj)
{
    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
}

// ============================ Evaluate ==============================

void cKit5D :: EvalExactConstraint(int index, int *algvar, double &c)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Single constraint evaluation.
  if (index == 0){
      c = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  }
  else if (index == 1){
      c = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  }
  else if (index == 2){
      c = (10*x[0] + x[1])/7 - 1;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cKit5D :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
  approxc[1] = 1;
  approxc[2] = 1;
}

// -------------------------------------------------------------------------
// Class c3BarTrussC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ c3BarTrussC ===============================
c3BarTrussC :: c3BarTrussC(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = Low[1] = 0.1;
  Upp[0] = Upp[1] = 5.0;
}

// ============================ Evaluate ===============================

void c3BarTrussC :: Evaluate(cVector & x, cVector & c, cVector & fobjs)
{
  // Renaming design variable
  cVector & A = x;

  // Maximum permissible stresses in compression and tension
  double sigma_c {15}; 
  double sigma_t {20};

  // Total volume of the truss
  double H = 1.0;
  double V = H * (2 * sqrt(2) * A[0] + A[1]);

  // Running Analysis
  double sigma [3];
  Analysis(A, sigma);

  // Constraints evaluation
  c[0] = sigma[0] / sigma_t - 1;
  c[1] = sigma[1] / sigma_t - 1;
  c[2] = abs(sigma[2]) / sigma_c - 1;

  // Objetive function evaluation
  fobjs[0] = V;
}

// ============================ FindPosition ===============================

int c3BarTrussC :: FindPosition(fstream & stream, string line)
{
  string entry;
  while (getline(stream, entry))
  {
    if (entry == line) return stream.tellp();
  }
  return -1;
}

// ============================ Analysis ===============================

void c3BarTrussC :: Analysis(cVector & A, double * sigma)
{ 
  // Load Parameters
  double P = 20;
  
  // Analytical Stresses
  sigma[0] = P * (A[1] + sqrt(2) * A[0]) / (sqrt(2) * A[0] * A[0] + 2 * A[0] * A[1]);
  sigma[1] = P * (1 / (A[0] + sqrt(2) * A[1]));
  sigma[2] = -P * (A[1] / (sqrt(2) * A[0] * A[0] + 2 * A[0] * A[1]));
}

// ============================ c3BarTrussCFAST :: Analysis ===============================

void c3BarTrussCFAST :: Analysis(cVector & A, double * sigma)
{
  // Files Base Name
  string base_name = "ThreeBarTrussFAST";

  // Opening Input File and Setting Float-point Format
  fstream dat_file (base_name + ".dat");
  dat_file << scientific << setprecision(5);

  // Searching for Keyword Position
  int start_position = FindPosition(dat_file, "%SECTION.BAR.GENERAL\r");

  // Replacing Area Values in Input File
  if (start_position != -1)
  {
    dat_file.seekp(start_position + 11);
    dat_file << A[0];
    dat_file.seekp(start_position + 80);
    dat_file << A[1];
  }
  else
  {
    cout << "Warning: It were not possible to replace Area Values.\n";
  }
  dat_file.close();

  // Calculating Stresses by Numerical Analysis
  system(("fast " + base_name + " -silent").c_str());

  // Reading Results
  double N [3];
  fstream pos_file (base_name + ".pos");
  string entry, trash;
  while (pos_file >> entry)
  {
    if (entry == "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR.DATA")
    {
      pos_file >> trash >> trash >> N[0];
      pos_file >> trash >> trash >> N[1];
      pos_file >> trash >> trash >> N[2];
      break;
    }
  }
  pos_file.close();
  sigma[0] = N[0] / A[0];
  sigma[1] = N[1] / A[1];
  sigma[2] = N[2] / A[0];
}

// ============================ c3BarTrussCABAQUS :: Analysis ===============================

void c3BarTrussCABAQUS :: Analysis(cVector & A, double * sigma)
{
  // Files Base Name
  string base_name = "ThreeBarTrussABAQUS";

  // Opening Input File and Setting Float-point Format
  fstream inp_file (base_name + ".inp");
  inp_file << scientific << setprecision(5);

  // Searching for Keyword Position
  int start_position_1 = FindPosition(inp_file, "*Solid Section, elset=A1-A3, material=Metal\r");
  int start_position_2 = FindPosition(inp_file, "*Solid Section, elset=A2, material=Metal\r");

  // Replacing Area Values in Input File
  if (start_position_1 != -1 && start_position_2 != -1)
  {
    inp_file.seekp(start_position_1);
    inp_file << A[0];
    inp_file.seekp(start_position_2);
    inp_file << A[1];
  } 
  else
  {
    cout << "Warning: It were not possible to replace Area Values.\n";
  }
  inp_file.close();

  // Calculating Stresses by Numerical Analysis and Extracting Results
  system(("cmd.exe /c abaqus cae noGUI=" + base_name + ".py").c_str());

  // Reading Results
  fstream txt_file (base_name + ".txt");
  txt_file >> sigma[0] >> sigma[1] >> sigma[2];
  txt_file.close();
}

// ============================ c3BarTrussCDIANA :: Analysis ===============================

void c3BarTrussCDIANA :: Analysis(cVector & A, double * sigma)
{
  // Files Base Name
  string base_name = "ThreeBarTrussDIANA";

  // Opening Input File and Setting Float-point Format
  fstream dat_file (base_name + ".dat");
  dat_file << scientific << setprecision(5);

  // Searching for Keyword Position
  int start_position = FindPosition(dat_file, "'GEOMET'\r");

  // Replacing Area Values in Input File
  if (start_position != -1)
  {
    dat_file.seekp(start_position + 67);
    dat_file << A[0];
    dat_file.seekp(start_position + 147);
    dat_file << A[1];
  }
  else
  {
    cout << "Warning: It were not possible to replace Area Values.\n";
  }
  dat_file.close();
  
  // Calculating Stresses by Numerical Analysis
  system(("cmd.exe /c diana " + base_name).c_str());

  // Reading Results
  fstream tb_file (base_name + ".tb");
  string entry, trash;
  while (tb_file >> entry)
  {
    if (entry == "Elmnr")
    {
      tb_file >> trash >> trash >> trash >> trash >> sigma[0];
      tb_file >> trash >> trash >> trash >> trash >> sigma[1];
      tb_file >> trash >> trash >> trash >> trash >> sigma[2];
      break;
    }
  }
  tb_file.close();
}

// -------------------------------------------------------------------------
// Class c3BarTrussD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ c3BarTrussD ===============================
c3BarTrussD :: c3BarTrussD(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  double step = 1.0e-3;

  ListDim = new int[NumVar];
  int i,j;
  for (i = 0; i < NumVar; i++) ListDim[i] = round((5.0 - 0.1)/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (i = 0; i < NumVar; i++) List[i].Resize(ListDim[i]);

  for (i = 0; i < NumVar; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 0.1 + step*j;
}

// ============================ Evaluate ===============================

void c3BarTrussD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Truss parameters

    double P = 20; // Nodal load

    double syc = 15; // Maximum permissible stress in compression
    double syt = 20;  // Maximum permissible stress in tension

    double A1 = x[0]; // Area of truss 1, which is the same as truss 3
    double A2 = x[1]; // Area of truss 2

    double w = 2*A1*sqrt(2) + A2; // Total weigth of the truss (divided by the constant \rho)

    double s1 = P*(A2 + sqrt(2)*A1)/(sqrt(2)*A1*A1 + 2*A1*A2);
    double s2 = P*(1/(A1 + sqrt(2)*A2));
    double s3 = -P*(A2/(sqrt(2)*A1*A1 + 2*A1*A2));

    // Constraints evaluation.

    c[0] = s1/syt - 1;
    c[1] = s2/syt - 1;
    c[2] = s3/syc - 1;

    // Objetive function evaluation

    double fobj = w;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class c10BarTruss:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ c10BarTruss ===============================
c10BarTruss :: c10BarTruss(void)
{
  NumVar = 10;
  NumConstr = 14;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
    Low[i] = 64.52e-6;
    Upp[i] = 1.0e-1;
  }
}

// ============================ Evaluate ===============================

void c10BarTruss :: Evaluate(cVector & x, cVector & c, cVector & fobjs)
{
  // Renaming design variable
  cVector & A = x;

  // Maximum displacement and stress
  double v_a {50.8e-3}; // m
  double sigma_a {172.37e6}; // Pa

  // Total volume of the truss
  double L {9.144}; // m
  double rho {2770}; // kg/m³
  double volume {0}; // m³
  for (int i = 0; i < 6; i++)
  {
    volume += A[i] * L;
  }
  for (int i = 0; i < 4; i++)
  {
    volume += A[6 + i] * L * sqrt(2);
  }
  double mass {rho * volume}; // kg

  // Running Analysis
  double v [4], sigma [10];
  Analysis(A, v, sigma);

  // Constraints evaluation
  for (int i = 0; i < 10; i++)
  {
    c[i] = abs(sigma[i]) / sigma_a - 1;
  }
  for (int i = 0; i < 4; i++)
  {
    c[10 + i] = abs(v[i]) / v_a - 1;
  }

  // Objetive function evaluation
  fobjs[0] = mass;
}

// ============================ FindPosition ===============================

int c10BarTruss :: FindPosition(fstream & stream, string line)
{
  string entry;
  while (getline(stream, entry))
  {
    if (entry == line) return stream.tellp();
  }
  return -1;
}

// ============================ FindPosition ===============================

void c10BarTruss :: ReplaceAreas(cVector & A, string base_name, int offset, int jump)
{
  fstream dat_file (base_name + ".dat");
  dat_file << scientific << setprecision(5);

  // Searching for Keyword Position
  int start_position = FindPosition(dat_file, "%SECTION.BAR.GENERAL\r");

  // Replacing Area Values in Input File
  if (start_position != -1)
  {
    for (int i = 0; i < NumVar; i++)
    {
      dat_file.seekp(start_position + offset + i * jump);
      dat_file << A[i];
    }
  }
  else
  {
    cout << "Warning: It were not possible to replace Area Values.\n";
  }
  dat_file.close();
}


// ============================ Analysis ===============================

void c10BarTruss :: Analysis(cVector & A, double * v, double * sigma)
{}

// ============================ c10BarTrussFAST :: Analysis ===============================

void c10BarTrussFAST :: Analysis(cVector & A, double * v, double * sigma)
{
  // Files Base Name
  string base_name = "TenBarTrussFAST";

  // Replacing Area Values in Input File
  ReplaceAreas(A, base_name, 13, 70);

  // Calculating Stresses by Numerical Analysis
  system(("fast " + base_name + " -silent").c_str());

  // Opening Output File
  fstream pos_file (base_name + ".pos");
  string entry, trash;

  // Reading Displacements
  while (pos_file >> entry)
  {
    if (entry == "%RESULT.CASE.STEP.NODAL.DISPLACEMENT")
    {
      pos_file >> trash >> trash >> trash >> trash;
      for (int i = 0; i < 4; i++)
      {
        pos_file >> v[i] >> trash >> trash >> trash >> trash >> trash >> trash;
      }
      break;
    }
  }

  // Reading Stresses
  double N [10];
  while (pos_file >> entry)
  {
    if (entry == "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR.DATA")
    {
      for (int i = 0; i < 10; i++)
      {
        pos_file >> trash >> trash >> N[i];
        sigma[i] = N[i] / A[i];
      }
      break;
    }
  }

  // Closing Input File
  pos_file.close();
}

// ============================ c10BarTrussABAQUS :: Analysis ===============================

void c10BarTrussABAQUS :: Analysis(cVector & A, double * v, double * sigma)
{
  // Files Base Name
  string base_name = "TenBarTrussABAQUS";

  // Opening Input File and Setting Float-point Format
  fstream inp_file (base_name + ".inp");
  inp_file << scientific << setprecision(5);

  // Searching for Keyword Position
  double start_position[10];
  start_position[0] = FindPosition(inp_file, "*Solid Section, elset=A1, material=Material-1\r");
  start_position[1] = FindPosition(inp_file, "*Solid Section, elset=A2, material=Material-1\r");
  start_position[2] = FindPosition(inp_file, "*Solid Section, elset=A3, material=Material-1\r");
  start_position[3] = FindPosition(inp_file, "*Solid Section, elset=A4, material=Material-1\r");
  start_position[4] = FindPosition(inp_file, "*Solid Section, elset=A5, material=Material-1\r");
  start_position[5] = FindPosition(inp_file, "*Solid Section, elset=A6, material=Material-1\r");
  start_position[6] = FindPosition(inp_file, "*Solid Section, elset=A7, material=Material-1\r");
  start_position[7] = FindPosition(inp_file, "*Solid Section, elset=A8, material=Material-1\r");
  start_position[8] = FindPosition(inp_file, "*Solid Section, elset=A9, material=Material-1\r");
  start_position[9] = FindPosition(inp_file, "*Solid Section, elset=A10, material=Material-1\r");

  // Replacing Area Values in Input File
  for (int i = 0; i < 10; i++){
  if (start_position[i] != -1)
  {
    inp_file.seekp(start_position[i]);
    inp_file << A[i];
  } 
  else
  {
    cout << "Warning: It were not possible to replace Area Values.\n";
  }}
  inp_file.close();

  // Calculating Stresses by Numerical Analysis and Extracting Results
  system(("cmd.exe /c abaqus cae noGUI=" + base_name + ".py").c_str());

  // Reading Results
  fstream txt_file (base_name + ".txt");
  // Reading Displacements
  for (int i = 0; i < 4; i++)
  {
  txt_file >> v[i];
  }
  // Reading Stresses
  for (int i = 0; i < 10; i++)
  {
  txt_file >> sigma[i];
  }
  txt_file.close();
}

// ============================ c10BarTrussDIANA :: Analysis ===============================

void c10BarTrussDIANA :: Analysis(cVector & A, double * v, double * sigma)
{
  // Files Base Name
  string base_name = "TenBarTrussDIANA";

  // Replacing Area Values in Input File
  ReplaceAreas(A, base_name, 67, 80);
  
  // Calculating Stresses by Numerical Analysis
  system(("cmd.exe /c diana " + base_name).c_str());

  // Reading Results
  fstream tb_file (base_name + ".tb");
  string entry, trash;
  while (tb_file >> entry)
  {
    if (entry == "Elmnr")
    {
      tb_file >> trash >> trash >> trash >> trash >> sigma[0];
      tb_file >> trash >> trash >> trash >> trash >> sigma[1];
      tb_file >> trash >> trash >> trash >> trash >> sigma[2];
      break;
    }
  }
  tb_file.close();
}

// -------------------------------------------------------------------------
// Class c10BarTrussFrequency:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ c10BarTrussFrequency ===============================
c10BarTrussFrequency :: c10BarTrussFrequency(void)
{
  NumConstr = 1;
}

// ============================ Evaluate ===============================

void c10BarTrussFrequency :: Evaluate(cVector & x, cVector & c, cVector & fobjs)
{
  // Renaming design variable
  cVector & A = x;  

  // Maximum frequencies
  double omega_a [1] {14}; // Hz

  // Total volume of the truss
  double L {9.144}; // m
  double rho {2770}; // kg/m³
  double volume {0}; // m³
  for (int i = 0; i < 6; i++)
  {
    volume += A[i] * L;
  }
  for (int i = 0; i < 4; i++)
  {
    volume += A[6 + i] * L * sqrt(2);
  }
  double mass {rho * volume}; // kg

  // Running Analysis
  double omega [1];
  Analysis(A, omega);

  // Constraints evaluation
  for (int i = 0; i < 1; i++)
  {
    c[i] = omega_a[i] / omega[i] - 1;
  }

  // Objetive function evaluation
  fobjs[0] = mass;
}

void c10BarTrussFrequency :: Analysis(cVector & A, double * omega)
{}

// ============================ c10BarTrussFrequencyFAST :: Analysis ===============================

void c10BarTrussFrequencyFAST :: Analysis(cVector & A, double * omega)
{
  // Files Base Name
  string base_name = "TenBarTrussFrequencyFAST";

  // Opening Input File and Setting Float-point Format
  fstream dat_file (base_name + ".dat");
  dat_file << scientific << setprecision(5);

  // Searching for Keyword Position
  int start_position = FindPosition(dat_file, "%SECTION.BAR.GENERAL\r");

  // Replacing Area Values in Input File
  int offsets {13};
  if (start_position != -1)
  {
    for (int i = 0; i < NumVar; i++)
    {
      dat_file.seekp(start_position + offsets + i * 70);
      dat_file << A[i];
    }
  }
  else
  {
    cout << "Warning: It were not possible to replace Area Values.\n";
  }
  dat_file.close();

  // Calculating Stresses by Numerical Analysis
  system(("fast " + base_name + " -silent").c_str());

  // Opening Output File
  fstream pos_file (base_name + ".pos");
  string entry, trash;

  // Reading Frequencies
  while (pos_file >> entry)
  {
    if (entry == "%RESULT.CASE.STEP.NATURAL.FREQUENCY")
    {
      pos_file >> trash >> omega[0];
      break;
    }
  }

  // Converting Natural Frenquency
   omega[0] /= 2 * PI;

  // Closing Input File
  pos_file.close();
}


// -------------------------------------------------------------------------
// Class cNowackiBeamC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cNowackiBeamC ===============================
cNowackiBeamC :: cNowackiBeamC(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = 5;
  Upp[0] = 50;
  Low[1] = 50;
  Upp[1] = 250;
}

// ============================ Evaluate ===============================

void cNowackiBeamC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Beam parameters

    double l = 1500;   // Length of the beam, mm
    double b = x[0];   // Beam width, mm
    double h = x[1];   // Beam height, mm
    double F = 5000;   // Tip load, N

    double E  = 216620; // Young's modulus, MPa

    double sy   = 240;  // Yield stress, MPa
    double dmax = 5;    // Maximum displacement, mm

    // Auxiliary variables

    double Iy   = b*h*h*h/12; // Cross-sectional inertia
    double Area = b*h;        // Cross-sectional area

    // Constraints evaluation.

    c[0] = F*l*l*l/(3*E*Iy*dmax) - 1; // Maximum displacement constraint
    c[1] = 6*F*l/(b*h*h*sy) - 1;      // Maximum bending stress constraint

    // Objetive function evaluation

    double fobj = Area;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cNowackiBeamD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cNowackiBeamD ===============================
cNowackiBeamD :: cNowackiBeamD(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  double step = 1.0e-3;

  ListDim = new int[NumVar];
  int i,j;
  for (i = 0; i < 1; i++) ListDim[i] = round((50.0 - 5.0)/step) + 1;
  for (i = 1; i < 2; i++) ListDim[i] = round((250.0 - 50.0)/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (i = 0; i < NumVar; i++) List[i].Resize(ListDim[i]);

  for (i = 0; i < 1; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 5.0 + step*j;
  for (i = 1; i < 2; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 50.0 + step*j;
}

// ============================ Evaluate ===============================

void cNowackiBeamD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Beam parameters

    double l = 1500;   // Length of the beam, mm
    double b = x[0];   // Beam width, mm
    double h = x[1];   // Beam height, mm
    double F = 5000;   // Tip load, N

    double E  = 216620; // Young's modulus, MPa

    double sy   = 240;  // Yield stress, MPa
    double dmax = 5;    // Maximum displacement, mm

    // Auxiliary variables

    double Iy   = b*h*h*h/12; // Cross-sectional inertia
    double Area = b*h;        // Cross-sectional area

    // Constraints evaluation.

    c[0] = F*l*l*l/(3*E*Iy*dmax) - 1; // Maximum displacement constraint
    c[1] = 6*F*l/(b*h*h*sy) - 1;      // Maximum bending stress constraint

    // Objetive function evaluation

    double fobj = Area;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cBeamC
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBeamC ===============================

cBeamC :: cBeamC(void)
{
  NumVar = 10;
  NumConstr = 11;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<5; i++)
  {
    Low[i] = 35;
    Upp[i] = 65;
  }

  for (int i = 5; i<10; i++)
  {
    Low[i] = 1;
    Upp[i] = 4;
  }
}

// ============================ Evaluate ==============================

void cBeamC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Beam rotation.

  double teta[5];
  teta[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    teta[i] = ((P*l)/(E*I[i]))*(L + (l/2.0) - aux[i]) + teta[i-1];
  }

  // Beam displacement.

  double displ[5];
  displ[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    displ[i] = ((P*l*l)/(2.0*E*I[i])*(L-aux[i]+(2.0*l/3.0))) + (teta[i-1]*l) + displ[i-1];
  }

  // Stress constraints.

  int nc = 0;
  for (int i = 0; i < 5; i++)
  {
    c[nc++] = sig[i]/Sall - 1.0;
  }

  // Geometry constraints.

  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    c[nc++] = h/(20.0*b) - 1.0;
  }

  // Displacement constraint.

  c[nc++] = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;

  // Objetive function evaluation.

  fobjs[0] = (x[0]*x[5] + x[1]*x[6] + x[2]*x[7] + x[3]*x[8] + x[4]*x[9])*l;
}

// ============================ Evaluate ==============================

void cBeamC :: EvalExactConstraint(int index, cVector &x, double &c)
{
  // Decodification of problem variables.

  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Single constraint evaluation.
  if (index < 5){
      c = sig[index]/Sall - 1.0;
  }
  else if (index >=5 && index < 10){
      double h = x[index - 5];
      double b = x[index];
      c = h/(20.0*b) - 1.0;
  }
  else if (index == 10){
      c = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cBeamC :: GetApproxConstr(bool* approxc)
{
  approxc[0]  = 0;
  approxc[1]  = 0;
  approxc[2]  = 0;
  approxc[3]  = 0;
  approxc[4]  = 0;
  approxc[5]  = 0;
  approxc[6]  = 0;
  approxc[7]  = 0;
  approxc[8]  = 0;
  approxc[9]  = 0;
  approxc[10] = 0;
  approxc[11] = 0;
}

// -------------------------------------------------------------------------
// Class cBeamD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBeamD ===============================
cBeamD :: cBeamD(void)
{
  NumVar = 10;
  NumConstr = 11;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  double step = 1.0e-3;

  ListDim = new int[NumVar];
  int i,j;
  for (i = 0; i < 5; i++) ListDim[i] = round((65 - 35)/step) + 1;
  for (i = 5; i < 10; i++) ListDim[i] = round((4 - 1)/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (i = 0; i < NumVar; i++) List[i].Resize(ListDim[i]);

  for (i = 0; i < 5; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 35 + step*j;

  for (i = 5; i < 10; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 1 + step*j;

}
// ========================= Evaluate ================================

void cBeamD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Beam rotation.

  double teta[5];
  teta[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    teta[i] = ((P*l)/(E*I[i]))*(L + (l/2.0) - aux[i]) + teta[i-1];
  }

  // Beam displacement.

  double displ[5];
  displ[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    displ[i] = ((P*l*l)/(2.0*E*I[i])*(L-aux[i]+(2.0*l/3.0))) + (teta[i-1]*l) + displ[i-1];
  }

  // Stress constraints.

  int nc = 0;
  for (int i = 0; i < 5; i++)
  {
    c[nc++] = sig[i]/Sall - 1.0;
  }

  // Geometry constraints.

  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    c[nc++] = h/(20.0*b) - 1.0;
  }

  // Displacement constraint.

  c[nc++] = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;

  // Objetive function evaluation.

  double fobj = (x[0]*x[5] + x[1]*x[6] + x[2]*x[7] + x[3]*x[8] + x[4]*x[9])*l;

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cBeamD :: EvalExactConstraint(int index, int *algvar, double &c)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Single constraint evaluation.
  if (index < 5){
      c = sig[index]/Sall - 1.0;
  }
  else if (index >=5 && index < 10){
      double h = x[index - 5];
      double b = x[index];
      c = h/(20.0*b) - 1.0;
  }
  else if (index == 10){
      c = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cBeamD :: GetApproxConstr(bool* approxc)
{
  approxc[0]  = 0;
  approxc[1]  = 0;
  approxc[2]  = 0;
  approxc[3]  = 0;
  approxc[4]  = 0;
  approxc[5]  = 0;
  approxc[6]  = 0;
  approxc[7]  = 0;
  approxc[8]  = 0;
  approxc[9]  = 0;
  approxc[10] = 0;
  approxc[11] = 0;
}

// -------------------------------------------------------------------------
// Class cCONSTRC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cCONSTRC ===============================

cCONSTRC :: cCONSTRC(void)
{
  NumVar = 2;
  NumObj = 2;
  NumConstr = 2;

  Low = new double[NumVar];
  Upp = new double[NumVar];
  Low[0] = 0.1;
  Upp[0] = 1.0;
  Low[1] = 0.0;
  Upp[1] = 5.0;
 }

// ============================= Evaluate ================================

void cCONSTRC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Constraint evaluation.

  c[0] =  6-x[1]-(9*x[0]);
  c[1] =  1+x[1]-9*x[0];

  // Objective functions evaluation.

  fobjs[0] = x[0];
  fobjs[1] = (1+ x[1])/(x[0]);

 // cout << "Fobjs " << fobjs[0] << "  " << fobjs[1] << endl;

}

// -------------------------------------------------------------------------
// Class cTNKC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cTNKC ===============================

cTNKC :: cTNKC(void)
{
  NumVar = 2;
  NumObj = 2;
  NumConstr = 2;

  Low = new double[NumVar];
  Upp = new double[NumVar];
  Low[0] = 0;
  Upp[0] = PI;
  Low[1] = 0;
  Upp[1] = PI;

 }

// ============================= Evaluate ================================

void cTNKC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Constraint evaluation.

  c[0] =  -pow(x[0],2) -pow(x[1],2) + 1 + 0.1*cos(16*atan(x[0]/x[1]));
  c[1] =  pow(x[0] - 0.5, 2)/0.5 + pow(x[1]-0.5, 2)/0.5 - 1;

  // Objective functions evaluation.

  fobjs[0] = x[0];
  fobjs[1] = x[1];

//  cout << "Fobjs " << fobjs[0] << "  " << fobjs[1] << endl;

}

// -------------------------------------------------------------------------
// Class cZDT6C:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cZDT6C ===============================

cZDT6C :: cZDT6C(void)
{
  NumVar = 10;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
      Low[i] = 0.0;
      Upp[i] = 1.0;
  }
}

// ============================= Evaluate ================================

void cZDT6C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  fobjs[0] = x[0];

  double aux = 0.0;

  for (int i = 1; i < NumVar; i++)
  {
      aux += x[i];
  }

  aux = 1 + (aux*9)/(NumVar-1);

  fobjs[1] = aux*(1-pow((fobjs[0]/aux),2));

}

// -------------------------------------------------------------------------
// Class cSCHC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cSCHC ===============================

cSCHC :: cSCHC(void)
{
  NumVar = 1;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = -1000;
  Upp[0] = 1000;

}

// ============================= Evaluate ================================

void cSCHC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  fobjs[0] = pow(x[0], 2);

  fobjs[1] = pow((x[0]-2),2);
}

// -------------------------------------------------------------------------
// Class cZDT1C:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cZDT1C ===============================

cZDT1C :: cZDT1C(void)
{
  NumVar = 30;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
      Low[i] = 0.0;
      Upp[i] = 1.0;
  }
}

// ============================= Evaluate ================================

void cZDT1C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  fobjs[0] = x[0];

  double aux = 0.0;

  for (int i = 1; i < NumVar; i++)
  {
      aux += x[i];
  }

  double aux2 = 0.0;

  aux2 = 1 + (9/(NumVar-1))*aux;

  fobjs[1] = aux2*(1-sqrt((fobjs[0]/aux2)));
}



// -------------------------------------------------------------------------
// Class cKURC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cKURC ===============================

cKURC :: cKURC(void)
{
  NumVar = 3;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
      Low[i] = -5;
      Upp[i] = 5;
  }
}

// ============================= Evaluate ================================

void cKURC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  double aux1 = 0.0;
  double aux2 = 0.0;

  aux1 = -0.2*sqrt((pow(x[0],2)+pow(x[1],2)));
  aux2 = -0.2*sqrt((pow(x[1],2)+pow(x[2],2)));

  fobjs[0] = -10*(exp(aux1)+exp(aux2));

  double aux3 = 0.0;

  for (int i = 0; i < NumVar; i++)
  {
      aux3 += pow(fabs(x[i]),0.8)+5*sin(pow(x[i],3));
  }

  fobjs[1] = aux3;
}

// =========================== End of file =================================

