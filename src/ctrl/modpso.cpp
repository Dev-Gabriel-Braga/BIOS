// -------------------------------------------------------------------------
// modpso.cpp - Implementation of the modified PSO class.
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
// Created:      22-Apr-2014    Elias Saraiva Barroso
//
// Modified:     
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "ctrl/modpso.h"
#include "prob/problem.h"
#include "ctrl/sel.h"
#include "sol/group.h"
#include "sol/particle.h"
#include "pen/penalty.h"
#include "main/utl.h"
#include "main/gblvar.h"
#include "main/gbldef.h"

// -------------------------------------------------------------------------
// cLaminatePSO class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cLaminatePSO ==============================

cLaminatePSO :: cLaminatePSO(void) : cStandardPSO( )
{
  Type = LAMINATE_PSO;
}

// ============================= ~cLaminatePSO ==============================

cLaminatePSO :: ~cLaminatePSO(void)
{
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Mutation =================================

void cLaminatePSO :: Mutation(cSwarm &son)
{
  // Perform the mutation operation.
  
  if (MutProb)
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++) 
      son[i]->Mutate(MutProb);
  }

  // Perform the special laminate mutation operation.
 
  if (LamMutProb[0] || LamMutProb[1] || LamMutProb[2])
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++) 
      son[i]->LamMutate(LamMutProb);
  }

  // Perform the layer swap operation.

  if (SwapProb)
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++) son[i]->Swap(SwapProb);
  }

  // Perform the layer addition operation.

  if (AddProb)
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++) son[i]->Add(AddProb);
  }

  // Perform the layer deletion operation.

  if (DelProb)
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++) son[i]->Delete(DelProb);
  }
}

// ============================= RandomRates ==============================

void cLaminatePSO :: RandomRates(void)
{
    if (!MutProb && MaxMut)
      MutProb = Utl::RandDouble(MinMut, MaxMut);

    if (!SwapProb && MaxSwap)
      SwapProb = Utl::RandDouble(MinSwap, MaxSwap);

    if (!AddProb && MaxAdd)
      AddProb = Utl::RandDouble(MinAdd, MaxAdd);

    if (!DelProb && MaxDel)
      DelProb = Utl::RandDouble(MinDel, MaxDel);
}

// ======================================================= End of file =====
