// -------------------------------------------------------------------------
// lamnsgaIIm.cpp - implementation of cmodLaminateNSGAIIm class.
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
// Created:      30-Aug-2017    Marina Alves Maia
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

#include "ctrl/modlamnsgaII.h"
#include "prob/problem.h"
#include "ctrl/sel.h"
#include "sol/group.h"
#include "sol/individual.h"
#include "pen/penalty.h"
#include "main/utl.h"
#include "main/gblvar.h"
#include "main/gbldef.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cmodLaminateNSGA ==============================

cmodLaminateNSGAII :: cmodLaminateNSGAII(void) : cmodNSGAII( ), cLamAlg( )
{
  Type = modLAMINATE_NSGAII;
}

// ============================= ~cmodLaminateNSGA ==============================

cmodLaminateNSGAII :: ~cmodLaminateNSGAII(void)
{
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Mutation =================================

void cmodLaminateNSGAII :: Mutation(cPopulation &son)
{
  // Perform the mutation operation.

  if (MutProb)
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++)
      son[i]->Mutate(MutProb);
  }

  // Perform the laminate special mutation operation.

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

void cmodLaminateNSGAII :: RandomRates(void)
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
