// -------------------------------------------------------------------------
// saokrg.cpp - implementation of cSAOKRG class.
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
// Created:      30-Jul-2019    Marina Alves Maia
//
// Modified:     17-Dec-2019    Marina Alves Maia
//                              Code refactoring.
// -------------------------------------------------------------------------

#include <string>
#include <vector>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "prob/problem.h"
#include "ctrl/sel.h"
#include "sol/group.h"
#include "sol/individual.h"
#include "pen/penalty.h"
#include "main/utl.h"
#include "main/input.h"
#include "main/gblvar.h"
#include "main/gbldef.h"
#include "ctrl/saokrg.h"
#include "sur/krg.h"
#include "mlib/vec.h"
#include "sur/surr.h"
#include "ctrl/stdpso.h"
#include "ctrl/stdga.h"
#include "prob/probsurr.h"
#include "ctrl/sao.h"
#include "sur/samp.h"

// -------------------------------------------------------------------------
// Class SAOKRG:
//

// -------------------------------------------------------------------------
// Public methods:

// =============================== cSAOKRG =================================

cSAOKRG :: cSAOKRG(void)  : cSAO( )
{
  Type          = SAOKRG;
  CorrType      = GAUSS;
  HyperParamLow = -1.0;
  HyperParamUpp = 2.0;
}

// ============================= ~cSAOKRG ==================================

cSAOKRG :: ~cSAOKRG(void)
{
}

// ============================ ReadCorrType ===============================

void cSAOKRG :: ReadCorrType(istream &in)
{
  char signame[100];

  if (!Utl::ReadString(in, signame))
  {
    cout << "Error in the input of the sigma definition method." << endl;
    exit(0);
  }

  if(string(signame) == "Gauss" || string(signame) == "GAUSS"
          || string(signame) == "gauss")
  {
    CorrType = GAUSS;
  }
  else if(string(signame) == "Matern" || string(signame) == "MATERN"
          || string(signame) == "matern")
  {
    CorrType = MATERN52;
  }
  else
  {
    cout << "Unknown correlation function definition: " << signame << endl;
    exit(0);
  }
}

// ======================== ReadHyperParam =================================

void cSAOKRG :: ReadHyperParam(istream &in)
{
  if (!(in >> HyperParamLow) || !(in >> HyperParamUpp))
  {
    cout << "Error in the input of the hyperparameters (lower and upper bounds)" << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cSAOKRG :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cSAO :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("CORRELATION.TYPE"         ,  makeReadObj(cSAOKRG, ReadCorrType));
  im.Insert("KRG.HYPERPARAMETERS" ,  makeReadObj(cSAOKRG, ReadHyperParam));
}

// =========================== CreateSurrogate =============================

cSURR* cSAOKRG ::  CreateSurrogate(sSampData &sdata) 
{
  krg = new cKRG;
  krg->CreateModel(sdata,CorrType,HyperParamLow,HyperParamUpp);

  return krg;
}

// =========================== UpdateSurrogate =============================

void cSAOKRG ::  UpdateSurrogate(cVectorVec &nx, cVectorVec &ny) 
{
  krg->UpdateModel(CorrType,nx,ny);
}

// ======================================================= End of file =====
