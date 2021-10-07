// -------------------------------------------------------------------------
// samp.cpp - implementation of cSamp class.
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
// Created:      27-Mai-2018    Leonardo Gonçalves Ribeiro
//
// Modified:     16-Nov-2020    Elias Saraiva Barroso
//               Added iostream reading method to eSmpType. 
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
//#include <bits/stdc++.h>

#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "matvec.h"
#include "samp.h"

using namespace std;

// -------------------------------------------------------------------------
// Auxiliary methods:
//

// ================================= operator>> ============================

istream& operator>>(istream &in, eSampType &t)
{
  char sampname[100];

  if (!Utl::ReadString(in, sampname))
  {
    in.setstate(ios::badbit);
    return in;
  }

  if(string(sampname) == "RANDOM" || string(sampname) == "Random"
          || string(sampname) == "SRS")
      t = RANDOM;
  else if(string(sampname) == "LATINHYPERSQUARE" || string(sampname) == "LatinHyperSquare"
          || string(sampname) == "LHS")
      t = LHS;
  else if(string(sampname) == "NLATINHYPERSQUARE" || string(sampname) == "nLatinHyperSquare"
          || string(sampname) == "NLHS")
      t = NLHS;
  else if(string(sampname) == "HAMMERSLEY"
          || string(sampname) == "Hammersley" || string(sampname) == "HSS")
      t = HAMMERSLEY;
  else if(string(sampname) == "SOBOL" || string(sampname) == "SobolSequence"
          || string(sampname) == "SSS")
      t = SOBOL;
  else if(string(sampname) == "CENTROIDALVORONOITESSELATION" || string(sampname) == "CentroidalVoronoiTesselation"
          || string(sampname) == "CVT")
      t = CVT;
  else if(string(sampname) == "LATINIZEDCENTROIDALVORONOITESSELATION" || string(sampname) == "LatinizedCentroidalVoronoiTesselation"
          || string(sampname) == "LCVT")
      t = LCVT;
  else
  {
    in.setstate(ios::failbit);
    cout << "Unknown sampling method: " << sampname << endl;
  }

  return in;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ================================= cSamp ==================================

cSamp :: cSamp(void)
{
}

// ================================= ~cSamp =================================

cSamp :: ~cSamp()
{
}

// ================================= InitSample ==========t========================

void cSamp :: InitSample(eSampType samptype, int nv, int ns, vector<cVector> &sx)
{
  // Read sample data.

  NumSample = ns;  // Number of sampling points
  NumVar = nv;     // Number of variables (i.e. dimension)

  cVector low(NumVar);
  cVector upp(NumVar);

  // Read sampling points and normalize to range [0, 1].

  cVector xn(NumVar);

  CreateSample(samptype, sx);
}

// ================================= CreateSample ==================================

void cSamp :: CreateSample(eSampType samptype, vector<cVector> &sx)
{
  int seed;
  seed = get_seed( );
  if (samptype == HAMMERSLEY)
    CalcSampHammersley(NumVar, NumSample, sx);
  else if (samptype == LHS)
  {
    CalcSampLHS(NumVar, NumSample, seed, sx);
  }
  else if (samptype == RANDOM)
  {
    cout << "seed = " << seed << endl;
    CalcSampRandom(NumVar, NumSample, seed, sx);
  }
  else if (samptype == SOBOL)
  {
    CalcSampSobol(NumVar, NumSample, 0, sx);
  }
  else if (samptype == NLHS)
  {
      CalcSampNLHS(NumVar, NumSample, sx);
  }
  else if (samptype == CVT)
  {
      double r[NumVar*NumSample];
      int it_num;
      double it_diff;
      double energy;
      CalcSampCVT (NumVar, NumSample, 1000, 0, 0, 10000, 100, 10, &seed, r, &it_num, &it_diff, &energy, sx );
  }
  else if (samptype == LCVT)
  {
      double r[NumVar*NumSample];
      int it_num;
      double it_diff;
      double energy;
      CalcSampLCVT (NumVar, NumSample, 1000, 0, 0, 10000, 100, 10, &seed, r, &it_num, &it_diff, &energy, sx );
  }
}

// ================================= CalcSampHammersley ==================================

void cSamp :: CalcSampHammersley ( int m, int n , vector<cVector>  &sx)

//****************************************************************************80
//
//  Purpose:
//
//    HAMMERSLEY computes an element of a Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2016
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, 1960, pages 844-874.
//
//  Parameters:
//
//    Input, int I, the index of the element of the sequence.
//    0 <= I.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the "base" for the first component.
//    1 <= N.
//
//    Output, double HAMMERSLEY[M], the element of the sequence with index I.
//
{
  int d;
  int i1;
  int j;
  cVector prime_inv(m);
  cVector r(m);
  cVector t(m);
  cVector xs(m);

  t[0] = 0;
  for (int i = 0; i < NumSample; i++)
  {
      for ( j = 1; j < m; j++ )
      {
        t[j] = i;
      }
      //
      //  Carry out the computation.
      //

      prime_inv[0] = 1.0;
      for ( j = 1; j < m; j++ )
      {
        prime_inv[j] = 1.0 / ( double ) ( prime ( j ) );
      }

      r[0] = ( double ) ( i % ( n + 1 ) ) / ( double ) ( n );
      for ( j = 1; j < m; j++ )
      {
        r[j] = 0.0;
      }

      while ( 0 < i4vec_sum ( m, t ) )
      {
        for ( j = 1; j < m; j++ )
        {
          int tj = t[j];
          d = ( tj % prime ( j ) );
          r[j] = r[j] + ( double ) ( d ) * prime_inv[j];
          prime_inv[j] = prime_inv[j] / ( double ) ( prime ( j ) );
          t[j] = ( t[j] / prime ( j ) );
        }
      }

      for (j = 0; j < m; j++)
      {
        xs[j] = r[j];
      }

      xs.Print();
      sx.push_back(xs);
   }
}


// ================================= prime ==================================

int cSamp :: prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, pages 870-873.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641,
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739,
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829,
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923,
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007,
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109,
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187,
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309,
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411,
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cerr << "\n";
    cerr << "PRIME - Fatal error!\n";
    cerr << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}

// ================================= i4vec_sum ==================================

int cSamp :: i4vec_sum ( int n, cVector a)

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:get_seed();
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}

// ================================= CalcSampLHS ==================================

void cSamp :: CalcSampLHS(int dim_num, int point_num, int& seed, vector<cVector>& sx)

//****************************************************************************80
//
//  Purpose:
//
//    LATIN_RANDOM_NEW returns points in a Latin Random square.
//
//  Discussion:
//
//    In each spatial dimension, there will be exactly one
//    point whose coordinate value lies between consecutive
//    values in the list:
//
//      ( 0, 1, 2, ..., point_num ) / point_num
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input/output, int &SEED, a seed for UNIFORM.
//
//    Output, double LATIN_RANDOM_NEW[DIM_NUM,POINT_NUM], the points.
//
//{
//  int i;
//  int *perm;
//  double r;
//  double *x;
//  cVector xs(NumVar);

//  for (int j = 0; j < point_num; j++)
//  {

//  x = r8mat_uniform_01_new ( dim_num, point_num, seed );
////
////  For spatial dimension I,
////    pick a random permutation of 1 to POINT_NUM,
////    force the corresponding I-th components of X to lie in the
////    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
////

//  for ( i = 0; i < dim_num; i++ )
//  {
//    perm = perm_uniform_new ( point_num, seed );

//    xs[i] = ( ( ( double ) perm[j] ) + x[i+j*dim_num] )/( ( double ) point_num );
//    delete [] perm;
//  }
//  sx.push_back(xs);
//  }
//}

{
  int i;
  int j;
  int *perm;
  double *x;
  cVector xs(NumVar);

  x = r8mat_uniform_01_new ( dim_num, point_num, seed );
//
//  For spatial dimension I,
//    pick a random permutation of 1 to POINT_NUM,
//    force the corresponding I-th components of X to lie in the
//    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
//
  for ( i = 0; i < dim_num; i++ )
  {
    perm = perm_uniform_new ( point_num, seed );

    for ( j = 0; j < point_num; j++ )
    {
      x[i+j*dim_num] = ( ( ( double ) perm[j] ) + (double) x[i+j*dim_num] ) / ( ( double ) point_num );
    }
    delete [] perm;
  }

  for ( int j = 0; j < point_num; j++)
  {
  for ( i = 0; i < dim_num; i++ )
  {
    xs[i] = x[i+j*dim_num];
  }
  sx.push_back(xs);
  }
}

// ================================= r8mat_uniform_01_new ==================================

double* cSamp :: r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}

// ================================= perm_uniform_new ==================================

int* cSamp :: perm_uniform_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM_NEW selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM_NEW[N], a permutation of
//    (0, 1, ..., N-1).
//
{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
}

// ================================= i4_uniform_ab ==================================

int cSamp :: i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 )
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}



// ================================= CalcSampLHS ==================================

void cSamp :: CalcSampRandom(int dim_num, int point_num, int& seed, vector<cVector>& sx)

{
  int i;
  double *x;
  cVector xs(NumVar);

  x = r8mat_uniform_01_new ( dim_num, point_num, seed );

  for ( int j = 0; j < point_num; j++)
  {
  for ( i = 0; i < dim_num; i++ )
  {
    xs[i] = x[i+j*dim_num];
  }
  sx.push_back(xs);
  }
}

// ================================= CalcSampSobol ==================================

void cSamp :: CalcSampSobol( int m, int n, int skip , vector <cVector>&sx)

//****************************************************************************80
//
//  Purpose:
//
//    I8_SOBOL_GENERATE generates a Sobol dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points to generate.
//
//    Input, int SKIP, the number of initial points to skip.
//
//    Output, double I8_SOBOL_GENERATE[M*N], the points.
//
{
  int j;
  double *r;
  long long int seed;
  cVector xs(NumVar);

  r = new double[m*n];

  seed = ( long long ) skip;

  for ( j = 0; j < n; j++ )
  {
    i8_sobol ( m, &seed, r+m*j );
  }

  for ( int j = 0; j < n; j++)
  {
  for ( int i = 0; i < m; i++ )
  {
    xs[i] = r[i+j*m];
  }
  sx.push_back(xs);
  }
}

// ================================= i8_sobol ==================================

void cSamp :: i8_sobol ( int dim_num, long long int *seed, double quasi[ ] )

//****************************************************************************80
//
//  Purpose:
//
//    I8_SOBOL generates a new quasirandom Sobol vector with each call.
//
//  Discussion:
//
//    The routine adapts the ideas of Antonov and Saleev.
//
//    This routine uses LONG LONG INT for integers and DOUBLE for real values.
//
//    Thanks to Steffan Berridge for supplying (twice) the properly
//    formatted V data needed to extend the original routine's dimension
//    limit from 40 to 1111, 05 June 2007.
//
//    Thanks to Francis Dalaudier for pointing out that the range of allowed
//    values of DIM_NUM should start at 1, not 2!  17 February 2009.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2009
//
//  Author:
//
//    FORTRAN77 original version by Bennett Fox.
//    C++ version by John Burkardt
//
//  Reference:
//
//    IA Antonov, VM Saleev,
//    An Economic Method of Computing LP Tau-Sequences,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 19, 1980, pages 252 - 256.
//
//    Paul Bratley, Bennett Fox,
//    Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 1, pages 88-100, 1988.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Stephen Joe, Frances Kuo
//    Remark on Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 29, Number 1, pages 49-57, March 2003.
//
//    Ilya Sobol,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 16, pages 236-242, 1977.
//
//    Ilya Sobol, YL Levitan,
//    The Production of Points Uniformly Distributed in a Multidimensional
//    Cube (in Russian),
//    Preprint IPM Akad. Nauk SSSR,
//    Number 40, Moscow 1976.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 1 <= DIM_NUM <= 1111.
//
//    Input/output, long long int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, double QUASI[DIM_NUM], the next quasirandom vector.
//
{
# define DIM_MAX 40
# define DIM_MAX2 1111
# define LOG_MAX 62
//
//  Here, we have commented out the definition of ATMOST, because
//  in some cases, a compiler was complaining that the value of ATMOST could not
//  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
//  so as long as we set MAXCOL (below) to what we expect it should be, we
//  may be able to get around this difficulty.
//  JVB, 24 January 2006.
//
//static long long int atmost = 4611686018427387903;
//
  static int dim_num_save = 0;
  long long int i;
  bool includ[LOG_MAX];
  static bool initialized = false;
  long long int j;
  long long int j2;
  long long int k;
  long long int l;
  static long long int lastq[DIM_MAX2];
  long long int m;
  static long long int maxcol;
  long long int newv;
  static long long int poly[DIM_MAX2] =
  {
        1,    3,    7,   11,   13,   19,   25,   37,   59,   47,
       61,   55,   41,   67,   97,   91,  109,  103,  115,  131,
      193,  137,  145,  143,  241,  157,  185,  167,  229,  171,
      213,  191,  253,  203,  211,  239,  247,  285,  369,  299,
      301,  333,  351,  355,  357,  361,  391,  397,  425,  451,
      463,  487,  501,  529,  539,  545,  557,  563,  601,  607,
      617,  623,  631,  637,  647,  661,  675,  677,  687,  695,
      701,  719,  721,  731,  757,  761,  787,  789,  799,  803,
      817,  827,  847,  859,  865,  875,  877,  883,  895,  901,
      911,  949,  953,  967,  971,  973,  981,  985,  995, 1001,
     1019, 1033, 1051, 1063, 1069, 1125, 1135, 1153, 1163, 1221,
     1239, 1255, 1267, 1279, 1293, 1305, 1315, 1329, 1341, 1347,
     1367, 1387, 1413, 1423, 1431, 1441, 1479, 1509, 1527, 1531,
     1555, 1557, 1573, 1591, 1603, 1615, 1627, 1657, 1663, 1673,
     1717, 1729, 1747, 1759, 1789, 1815, 1821, 1825, 1849, 1863,
     1869, 1877, 1881, 1891, 1917, 1933, 1939, 1969, 2011, 2035,
     2041, 2053, 2071, 2091, 2093, 2119, 2147, 2149, 2161, 2171,
     2189, 2197, 2207, 2217, 2225, 2255, 2257, 2273, 2279, 2283,
     2293, 2317, 2323, 2341, 2345, 2363, 2365, 2373, 2377, 2385,
     2395, 2419, 2421, 2431, 2435, 2447, 2475, 2477, 2489, 2503,
     2521, 2533, 2551, 2561, 2567, 2579, 2581, 2601, 2633, 2657,
     2669, 2681, 2687, 2693, 2705, 2717, 2727, 2731, 2739, 2741,
     2773, 2783, 2793, 2799, 2801, 2811, 2819, 2825, 2833, 2867,
     2879, 2881, 2891, 2905, 2911, 2917, 2927, 2941, 2951, 2955,
     2963, 2965, 2991, 2999, 3005, 3017, 3035, 3037, 3047, 3053,
     3083, 3085, 3097, 3103, 3159, 3169, 3179, 3187, 3205, 3209,
     3223, 3227, 3229, 3251, 3263, 3271, 3277, 3283, 3285, 3299,
     3305, 3319, 3331, 3343, 3357, 3367, 3373, 3393, 3399, 3413,
     3417, 3427, 3439, 3441, 3475, 3487, 3497, 3515, 3517, 3529,
     3543, 3547, 3553, 3559, 3573, 3589, 3613, 3617, 3623, 3627,
     3635, 3641, 3655, 3659, 3669, 3679, 3697, 3707, 3709, 3713,
     3731, 3743, 3747, 3771, 3791, 3805, 3827, 3833, 3851, 3865,
     3889, 3895, 3933, 3947, 3949, 3957, 3971, 3985, 3991, 3995,
     4007, 4013, 4021, 4045, 4051, 4069, 4073, 4179, 4201, 4219,
     4221, 4249, 4305, 4331, 4359, 4383, 4387, 4411, 4431, 4439,
     4449, 4459, 4485, 4531, 4569, 4575, 4621, 4663, 4669, 4711,
     4723, 4735, 4793, 4801, 4811, 4879, 4893, 4897, 4921, 4927,
     4941, 4977, 5017, 5027, 5033, 5127, 5169, 5175, 5199, 5213,
     5223, 5237, 5287, 5293, 5331, 5391, 5405, 5453, 5523, 5573,
     5591, 5597, 5611, 5641, 5703, 5717, 5721, 5797, 5821, 5909,
     5913, 5955, 5957, 6005, 6025, 6061, 6067, 6079, 6081, 6231,
     6237, 6289, 6295, 6329, 6383, 6427, 6453, 6465, 6501, 6523,
     6539, 6577, 6589, 6601, 6607, 6631, 6683, 6699, 6707, 6761,
     6795, 6865, 6881, 6901, 6923, 6931, 6943, 6999, 7057, 7079,
     7103, 7105, 7123, 7173, 7185, 7191, 7207, 7245, 7303, 7327,
     7333, 7355, 7365, 7369, 7375, 7411, 7431, 7459, 7491, 7505,
     7515, 7541, 7557, 7561, 7701, 7705, 7727, 7749, 7761, 7783,
     7795, 7823, 7907, 7953, 7963, 7975, 8049, 8089, 8123, 8125,
     8137, 8219, 8231, 8245, 8275, 8293, 8303, 8331, 8333, 8351,
     8357, 8367, 8379, 8381, 8387, 8393, 8417, 8435, 8461, 8469,
     8489, 8495, 8507, 8515, 8551, 8555, 8569, 8585, 8599, 8605,
     8639, 8641, 8647, 8653, 8671, 8675, 8689, 8699, 8729, 8741,
     8759, 8765, 8771, 8795, 8797, 8825, 8831, 8841, 8855, 8859,
     8883, 8895, 8909, 8943, 8951, 8955, 8965, 8999, 9003, 9031,
     9045, 9049, 9071, 9073, 9085, 9095, 9101, 9109, 9123, 9129,
     9137, 9143, 9147, 9185, 9197, 9209, 9227, 9235, 9247, 9253,
     9257, 9277, 9297, 9303, 9313, 9325, 9343, 9347, 9371, 9373,
     9397, 9407, 9409, 9415, 9419, 9443, 9481, 9495, 9501, 9505,
     9517, 9529, 9555, 9557, 9571, 9585, 9591, 9607, 9611, 9621,
     9625, 9631, 9647, 9661, 9669, 9679, 9687, 9707, 9731, 9733,
     9745, 9773, 9791, 9803, 9811, 9817, 9833, 9847, 9851, 9863,
     9875, 9881, 9905, 9911, 9917, 9923, 9963, 9973,10003,10025,
    10043,10063,10071,10077,10091,10099,10105,10115,10129,10145,
    10169,10183,10187,10207,10223,10225,10247,10265,10271,10275,
    10289,10299,10301,10309,10343,10357,10373,10411,10413,10431,
    10445,10453,10463,10467,10473,10491,10505,10511,10513,10523,
    10539,10549,10559,10561,10571,10581,10615,10621,10625,10643,
    10655,10671,10679,10685,10691,10711,10739,10741,10755,10767,
    10781,10785,10803,10805,10829,10857,10863,10865,10875,10877,
    10917,10921,10929,10949,10967,10971,10987,10995,11009,11029,
    11043,11045,11055,11063,11075,11081,11117,11135,11141,11159,
    11163,11181,11187,11225,11237,11261,11279,11297,11307,11309,
    11327,11329,11341,11377,11403,11405,11413,11427,11439,11453,
    11461,11473,11479,11489,11495,11499,11533,11545,11561,11567,
    11575,11579,11589,11611,11623,11637,11657,11663,11687,11691,
    11701,11747,11761,11773,11783,11795,11797,11817,11849,11855,
    11867,11869,11873,11883,11919,11921,11927,11933,11947,11955,
    11961,11999,12027,12029,12037,12041,12049,12055,12095,12097,
    12107,12109,12121,12127,12133,12137,12181,12197,12207,12209,
    12239,12253,12263,12269,12277,12287,12295,12309,12313,12335,
    12361,12367,12391,12409,12415,12433,12449,12469,12479,12481,
    12499,12505,12517,12527,12549,12559,12597,12615,12621,12639,
    12643,12657,12667,12707,12713,12727,12741,12745,12763,12769,
    12779,12781,12787,12799,12809,12815,12829,12839,12857,12875,
    12883,12889,12901,12929,12947,12953,12959,12969,12983,12987,
    12995,13015,13019,13031,13063,13077,13103,13137,13149,13173,
    13207,13211,13227,13241,13249,13255,13269,13283,13285,13303,
    13307,13321,13339,13351,13377,13389,13407,13417,13431,13435,
    13447,13459,13465,13477,13501,13513,13531,13543,13561,13581,
    13599,13605,13617,13623,13637,13647,13661,13677,13683,13695,
    13725,13729,13753,13773,13781,13785,13795,13801,13807,13825,
    13835,13855,13861,13871,13883,13897,13905,13915,13939,13941,
    13969,13979,13981,13997,14027,14035,14037,14051,14063,14085,
    14095,14107,14113,14125,14137,14145,14151,14163,14193,14199,
    14219,14229,14233,14243,14277,14287,14289,14295,14301,14305,
    14323,14339,14341,14359,14365,14375,14387,14411,14425,14441,
    14449,14499,14513,14523,14537,14543,14561,14579,14585,14593,
    14599,14603,14611,14641,14671,14695,14701,14723,14725,14743,
    14753,14759,14765,14795,14797,14803,14831,14839,14845,14855,
    14889,14895,14909,14929,14941,14945,14951,14963,14965,14985,
    15033,15039,15053,15059,15061,15071,15077,15081,15099,15121,
    15147,15149,15157,15167,15187,15193,15203,15205,15215,15217,
    15223,15243,15257,15269,15273,15287,15291,15313,15335,15347,
    15359,15373,15379,15381,15391,15395,15397,15419,15439,15453,
    15469,15491,15503,15517,15527,15531,15545,15559,15593,15611,
    15613,15619,15639,15643,15649,15661,15667,15669,15681,15693,
    15717,15721,15741,15745,15765,15793,15799,15811,15825,15835,
    15847,15851,15865,15877,15881,15887,15899,15915,15935,15937,
    15955,15973,15977,16011,16035,16061,16069,16087,16093,16097,
    16121,16141,16153,16159,16165,16183,16189,16195,16197,16201,
    16209,16215,16225,16259,16265,16273,16299,16309,16355,16375,
    16381 };
  static double recipd;
  static long long int seed_save = - 1;
  long long int seed_temp;
  static long long int v[DIM_MAX2][LOG_MAX];

  if ( !initialized || dim_num != dim_num_save )
  {
    initialized = true;
    for ( i = 0; i < DIM_MAX2; i++ )
    {
      for ( j = 0; j < LOG_MAX; j++ )
      {
        v[i][j] = 0;
      }
    }
//
//  Initialize (part of) V.
//
    v[0][0] = 1;
    v[1][0] = 1;
    v[2][0] = 1;
    v[3][0] = 1;
    v[4][0] = 1;
    v[5][0] = 1;
    v[6][0] = 1;
    v[7][0] = 1;
    v[8][0] = 1;
    v[9][0] = 1;
    v[10][0] = 1;
    v[11][0] = 1;
    v[12][0] = 1;
    v[13][0] = 1;
    v[14][0] = 1;
    v[15][0] = 1;
    v[16][0] = 1;
    v[17][0] = 1;
    v[18][0] = 1;
    v[19][0] = 1;
    v[20][0] = 1;
    v[21][0] = 1;
    v[22][0] = 1;
    v[23][0] = 1;
    v[24][0] = 1;
    v[25][0] = 1;
    v[26][0] = 1;
    v[27][0] = 1;
    v[28][0] = 1;
    v[29][0] = 1;
    v[30][0] = 1;
    v[31][0] = 1;
    v[32][0] = 1;
    v[33][0] = 1;
    v[34][0] = 1;
    v[35][0] = 1;
    v[36][0] = 1;
    v[37][0] = 1;
    v[38][0] = 1;
    v[39][0] = 1;
    v[40][0] = 1;
    v[41][0] = 1;
    v[42][0] = 1;
    v[43][0] = 1;
    v[44][0] = 1;
    v[45][0] = 1;
    v[46][0] = 1;
    v[47][0] = 1;
    v[48][0] = 1;
    v[49][0] = 1;
    v[50][0] = 1;
    v[51][0] = 1;
    v[52][0] = 1;
    v[53][0] = 1;
    v[54][0] = 1;
    v[55][0] = 1;
    v[56][0] = 1;
    v[57][0] = 1;
    v[58][0] = 1;
    v[59][0] = 1;
    v[60][0] = 1;
    v[61][0] = 1;
    v[62][0] = 1;
    v[63][0] = 1;
    v[64][0] = 1;
    v[65][0] = 1;
    v[66][0] = 1;
    v[67][0] = 1;
    v[68][0] = 1;
    v[69][0] = 1;
    v[70][0] = 1;
    v[71][0] = 1;
    v[72][0] = 1;
    v[73][0] = 1;
    v[74][0] = 1;
    v[75][0] = 1;
    v[76][0] = 1;
    v[77][0] = 1;
    v[78][0] = 1;
    v[79][0] = 1;
    v[80][0] = 1;
    v[81][0] = 1;
    v[82][0] = 1;
    v[83][0] = 1;
    v[84][0] = 1;
    v[85][0] = 1;
    v[86][0] = 1;
    v[87][0] = 1;
    v[88][0] = 1;
    v[89][0] = 1;
    v[90][0] = 1;
    v[91][0] = 1;
    v[92][0] = 1;
    v[93][0] = 1;
    v[94][0] = 1;
    v[95][0] = 1;
    v[96][0] = 1;
    v[97][0] = 1;
    v[98][0] = 1;
    v[99][0] = 1;
    v[100][0] = 1;
    v[101][0] = 1;
    v[102][0] = 1;
    v[103][0] = 1;
    v[104][0] = 1;
    v[105][0] = 1;
    v[106][0] = 1;
    v[107][0] = 1;
    v[108][0] = 1;
    v[109][0] = 1;
    v[110][0] = 1;
    v[111][0] = 1;
    v[112][0] = 1;
    v[113][0] = 1;
    v[114][0] = 1;
    v[115][0] = 1;
    v[116][0] = 1;
    v[117][0] = 1;
    v[118][0] = 1;
    v[119][0] = 1;
    v[120][0] = 1;
    v[121][0] = 1;
    v[122][0] = 1;
    v[123][0] = 1;
    v[124][0] = 1;
    v[125][0] = 1;
    v[126][0] = 1;
    v[127][0] = 1;
    v[128][0] = 1;
    v[129][0] = 1;
    v[130][0] = 1;
    v[131][0] = 1;
    v[132][0] = 1;
    v[133][0] = 1;
    v[134][0] = 1;
    v[135][0] = 1;
    v[136][0] = 1;
    v[137][0] = 1;
    v[138][0] = 1;
    v[139][0] = 1;
    v[140][0] = 1;
    v[141][0] = 1;
    v[142][0] = 1;
    v[143][0] = 1;
    v[144][0] = 1;
    v[145][0] = 1;
    v[146][0] = 1;
    v[147][0] = 1;
    v[148][0] = 1;
    v[149][0] = 1;
    v[150][0] = 1;
    v[151][0] = 1;
    v[152][0] = 1;
    v[153][0] = 1;
    v[154][0] = 1;
    v[155][0] = 1;
    v[156][0] = 1;
    v[157][0] = 1;
    v[158][0] = 1;
    v[159][0] = 1;
    v[160][0] = 1;
    v[161][0] = 1;
    v[162][0] = 1;
    v[163][0] = 1;
    v[164][0] = 1;
    v[165][0] = 1;
    v[166][0] = 1;
    v[167][0] = 1;
    v[168][0] = 1;
    v[169][0] = 1;
    v[170][0] = 1;
    v[171][0] = 1;
    v[172][0] = 1;
    v[173][0] = 1;
    v[174][0] = 1;
    v[175][0] = 1;
    v[176][0] = 1;
    v[177][0] = 1;
    v[178][0] = 1;
    v[179][0] = 1;
    v[180][0] = 1;
    v[181][0] = 1;
    v[182][0] = 1;
    v[183][0] = 1;
    v[184][0] = 1;
    v[185][0] = 1;
    v[186][0] = 1;
    v[187][0] = 1;
    v[188][0] = 1;
    v[189][0] = 1;
    v[190][0] = 1;
    v[191][0] = 1;
    v[192][0] = 1;
    v[193][0] = 1;
    v[194][0] = 1;
    v[195][0] = 1;
    v[196][0] = 1;
    v[197][0] = 1;
    v[198][0] = 1;
    v[199][0] = 1;
    v[200][0] = 1;
    v[201][0] = 1;
    v[202][0] = 1;
    v[203][0] = 1;
    v[204][0] = 1;
    v[205][0] = 1;
    v[206][0] = 1;
    v[207][0] = 1;
    v[208][0] = 1;
    v[209][0] = 1;
    v[210][0] = 1;
    v[211][0] = 1;
    v[212][0] = 1;
    v[213][0] = 1;
    v[214][0] = 1;
    v[215][0] = 1;
    v[216][0] = 1;
    v[217][0] = 1;
    v[218][0] = 1;
    v[219][0] = 1;
    v[220][0] = 1;
    v[221][0] = 1;
    v[222][0] = 1;
    v[223][0] = 1;
    v[224][0] = 1;
    v[225][0] = 1;
    v[226][0] = 1;
    v[227][0] = 1;
    v[228][0] = 1;
    v[229][0] = 1;
    v[230][0] = 1;
    v[231][0] = 1;
    v[232][0] = 1;
    v[233][0] = 1;
    v[234][0] = 1;
    v[235][0] = 1;
    v[236][0] = 1;
    v[237][0] = 1;
    v[238][0] = 1;
    v[239][0] = 1;
    v[240][0] = 1;
    v[241][0] = 1;
    v[242][0] = 1;
    v[243][0] = 1;
    v[244][0] = 1;
    v[245][0] = 1;
    v[246][0] = 1;
    v[247][0] = 1;
    v[248][0] = 1;
    v[249][0] = 1;
    v[250][0] = 1;
    v[251][0] = 1;
    v[252][0] = 1;
    v[253][0] = 1;
    v[254][0] = 1;
    v[255][0] = 1;
    v[256][0] = 1;
    v[257][0] = 1;
    v[258][0] = 1;
    v[259][0] = 1;
    v[260][0] = 1;
    v[261][0] = 1;
    v[262][0] = 1;
    v[263][0] = 1;
    v[264][0] = 1;
    v[265][0] = 1;
    v[266][0] = 1;
    v[267][0] = 1;
    v[268][0] = 1;
    v[269][0] = 1;
    v[270][0] = 1;
    v[271][0] = 1;
    v[272][0] = 1;
    v[273][0] = 1;
    v[274][0] = 1;
    v[275][0] = 1;
    v[276][0] = 1;
    v[277][0] = 1;
    v[278][0] = 1;
    v[279][0] = 1;
    v[280][0] = 1;
    v[281][0] = 1;
    v[282][0] = 1;
    v[283][0] = 1;
    v[284][0] = 1;
    v[285][0] = 1;
    v[286][0] = 1;
    v[287][0] = 1;
    v[288][0] = 1;
    v[289][0] = 1;
    v[290][0] = 1;
    v[291][0] = 1;
    v[292][0] = 1;
    v[293][0] = 1;
    v[294][0] = 1;
    v[295][0] = 1;
    v[296][0] = 1;
    v[297][0] = 1;
    v[298][0] = 1;
    v[299][0] = 1;
    v[300][0] = 1;
    v[301][0] = 1;
    v[302][0] = 1;
    v[303][0] = 1;
    v[304][0] = 1;
    v[305][0] = 1;
    v[306][0] = 1;
    v[307][0] = 1;
    v[308][0] = 1;
    v[309][0] = 1;
    v[310][0] = 1;
    v[311][0] = 1;
    v[312][0] = 1;
    v[313][0] = 1;
    v[314][0] = 1;
    v[315][0] = 1;
    v[316][0] = 1;
    v[317][0] = 1;
    v[318][0] = 1;
    v[319][0] = 1;
    v[320][0] = 1;
    v[321][0] = 1;
    v[322][0] = 1;
    v[323][0] = 1;
    v[324][0] = 1;
    v[325][0] = 1;
    v[326][0] = 1;
    v[327][0] = 1;
    v[328][0] = 1;
    v[329][0] = 1;
    v[330][0] = 1;
    v[331][0] = 1;
    v[332][0] = 1;
    v[333][0] = 1;
    v[334][0] = 1;
    v[335][0] = 1;
    v[336][0] = 1;
    v[337][0] = 1;
    v[338][0] = 1;
    v[339][0] = 1;
    v[340][0] = 1;
    v[341][0] = 1;
    v[342][0] = 1;
    v[343][0] = 1;
    v[344][0] = 1;
    v[345][0] = 1;
    v[346][0] = 1;
    v[347][0] = 1;
    v[348][0] = 1;
    v[349][0] = 1;
    v[350][0] = 1;
    v[351][0] = 1;
    v[352][0] = 1;
    v[353][0] = 1;
    v[354][0] = 1;
    v[355][0] = 1;
    v[356][0] = 1;
    v[357][0] = 1;
    v[358][0] = 1;
    v[359][0] = 1;
    v[360][0] = 1;
    v[361][0] = 1;
    v[362][0] = 1;
    v[363][0] = 1;
    v[364][0] = 1;
    v[365][0] = 1;
    v[366][0] = 1;
    v[367][0] = 1;
    v[368][0] = 1;
    v[369][0] = 1;
    v[370][0] = 1;
    v[371][0] = 1;
    v[372][0] = 1;
    v[373][0] = 1;
    v[374][0] = 1;
    v[375][0] = 1;
    v[376][0] = 1;
    v[377][0] = 1;
    v[378][0] = 1;
    v[379][0] = 1;
    v[380][0] = 1;
    v[381][0] = 1;
    v[382][0] = 1;
    v[383][0] = 1;
    v[384][0] = 1;
    v[385][0] = 1;
    v[386][0] = 1;
    v[387][0] = 1;
    v[388][0] = 1;
    v[389][0] = 1;
    v[390][0] = 1;
    v[391][0] = 1;
    v[392][0] = 1;
    v[393][0] = 1;
    v[394][0] = 1;
    v[395][0] = 1;
    v[396][0] = 1;
    v[397][0] = 1;
    v[398][0] = 1;
    v[399][0] = 1;
    v[400][0] = 1;
    v[401][0] = 1;
    v[402][0] = 1;
    v[403][0] = 1;
    v[404][0] = 1;
    v[405][0] = 1;
    v[406][0] = 1;
    v[407][0] = 1;
    v[408][0] = 1;
    v[409][0] = 1;
    v[410][0] = 1;
    v[411][0] = 1;
    v[412][0] = 1;
    v[413][0] = 1;
    v[414][0] = 1;
    v[415][0] = 1;
    v[416][0] = 1;
    v[417][0] = 1;
    v[418][0] = 1;
    v[419][0] = 1;
    v[420][0] = 1;
    v[421][0] = 1;
    v[422][0] = 1;
    v[423][0] = 1;
    v[424][0] = 1;
    v[425][0] = 1;
    v[426][0] = 1;
    v[427][0] = 1;
    v[428][0] = 1;
    v[429][0] = 1;
    v[430][0] = 1;
    v[431][0] = 1;
    v[432][0] = 1;
    v[433][0] = 1;
    v[434][0] = 1;
    v[435][0] = 1;
    v[436][0] = 1;
    v[437][0] = 1;
    v[438][0] = 1;
    v[439][0] = 1;
    v[440][0] = 1;
    v[441][0] = 1;
    v[442][0] = 1;
    v[443][0] = 1;
    v[444][0] = 1;
    v[445][0] = 1;
    v[446][0] = 1;
    v[447][0] = 1;
    v[448][0] = 1;
    v[449][0] = 1;
    v[450][0] = 1;
    v[451][0] = 1;
    v[452][0] = 1;
    v[453][0] = 1;
    v[454][0] = 1;
    v[455][0] = 1;
    v[456][0] = 1;
    v[457][0] = 1;
    v[458][0] = 1;
    v[459][0] = 1;
    v[460][0] = 1;
    v[461][0] = 1;
    v[462][0] = 1;
    v[463][0] = 1;
    v[464][0] = 1;
    v[465][0] = 1;
    v[466][0] = 1;
    v[467][0] = 1;
    v[468][0] = 1;
    v[469][0] = 1;
    v[470][0] = 1;
    v[471][0] = 1;
    v[472][0] = 1;
    v[473][0] = 1;
    v[474][0] = 1;
    v[475][0] = 1;
    v[476][0] = 1;
    v[477][0] = 1;
    v[478][0] = 1;
    v[479][0] = 1;
    v[480][0] = 1;
    v[481][0] = 1;
    v[482][0] = 1;
    v[483][0] = 1;
    v[484][0] = 1;
    v[485][0] = 1;
    v[486][0] = 1;
    v[487][0] = 1;
    v[488][0] = 1;
    v[489][0] = 1;
    v[490][0] = 1;
    v[491][0] = 1;
    v[492][0] = 1;
    v[493][0] = 1;
    v[494][0] = 1;
    v[495][0] = 1;
    v[496][0] = 1;
    v[497][0] = 1;
    v[498][0] = 1;
    v[499][0] = 1;
    v[500][0] = 1;
    v[501][0] = 1;
    v[502][0] = 1;
    v[503][0] = 1;
    v[504][0] = 1;
    v[505][0] = 1;
    v[506][0] = 1;
    v[507][0] = 1;
    v[508][0] = 1;
    v[509][0] = 1;
    v[510][0] = 1;
    v[511][0] = 1;
    v[512][0] = 1;
    v[513][0] = 1;
    v[514][0] = 1;
    v[515][0] = 1;
    v[516][0] = 1;
    v[517][0] = 1;
    v[518][0] = 1;
    v[519][0] = 1;
    v[520][0] = 1;
    v[521][0] = 1;
    v[522][0] = 1;
    v[523][0] = 1;
    v[524][0] = 1;
    v[525][0] = 1;
    v[526][0] = 1;
    v[527][0] = 1;
    v[528][0] = 1;
    v[529][0] = 1;
    v[530][0] = 1;
    v[531][0] = 1;
    v[532][0] = 1;
    v[533][0] = 1;
    v[534][0] = 1;
    v[535][0] = 1;
    v[536][0] = 1;
    v[537][0] = 1;
    v[538][0] = 1;
    v[539][0] = 1;
    v[540][0] = 1;
    v[541][0] = 1;
    v[542][0] = 1;
    v[543][0] = 1;
    v[544][0] = 1;
    v[545][0] = 1;
    v[546][0] = 1;
    v[547][0] = 1;
    v[548][0] = 1;
    v[549][0] = 1;
    v[550][0] = 1;
    v[551][0] = 1;
    v[552][0] = 1;
    v[553][0] = 1;
    v[554][0] = 1;
    v[555][0] = 1;
    v[556][0] = 1;
    v[557][0] = 1;
    v[558][0] = 1;
    v[559][0] = 1;
    v[560][0] = 1;
    v[561][0] = 1;
    v[562][0] = 1;
    v[563][0] = 1;
    v[564][0] = 1;
    v[565][0] = 1;
    v[566][0] = 1;
    v[567][0] = 1;
    v[568][0] = 1;
    v[569][0] = 1;
    v[570][0] = 1;
    v[571][0] = 1;
    v[572][0] = 1;
    v[573][0] = 1;
    v[574][0] = 1;
    v[575][0] = 1;
    v[576][0] = 1;
    v[577][0] = 1;
    v[578][0] = 1;
    v[579][0] = 1;
    v[580][0] = 1;
    v[581][0] = 1;
    v[582][0] = 1;
    v[583][0] = 1;
    v[584][0] = 1;
    v[585][0] = 1;
    v[586][0] = 1;
    v[587][0] = 1;
    v[588][0] = 1;
    v[589][0] = 1;
    v[590][0] = 1;
    v[591][0] = 1;
    v[592][0] = 1;
    v[593][0] = 1;
    v[594][0] = 1;
    v[595][0] = 1;
    v[596][0] = 1;
    v[597][0] = 1;
    v[598][0] = 1;
    v[599][0] = 1;
    v[600][0] = 1;
    v[601][0] = 1;
    v[602][0] = 1;
    v[603][0] = 1;
    v[604][0] = 1;
    v[605][0] = 1;
    v[606][0] = 1;
    v[607][0] = 1;
    v[608][0] = 1;
    v[609][0] = 1;
    v[610][0] = 1;
    v[611][0] = 1;
    v[612][0] = 1;
    v[613][0] = 1;
    v[614][0] = 1;
    v[615][0] = 1;
    v[616][0] = 1;
    v[617][0] = 1;
    v[618][0] = 1;
    v[619][0] = 1;
    v[620][0] = 1;
    v[621][0] = 1;
    v[622][0] = 1;
    v[623][0] = 1;
    v[624][0] = 1;
    v[625][0] = 1;
    v[626][0] = 1;
    v[627][0] = 1;
    v[628][0] = 1;
    v[629][0] = 1;
    v[630][0] = 1;
    v[631][0] = 1;
    v[632][0] = 1;
    v[633][0] = 1;
    v[634][0] = 1;
    v[635][0] = 1;
    v[636][0] = 1;
    v[637][0] = 1;
    v[638][0] = 1;
    v[639][0] = 1;
    v[640][0] = 1;
    v[641][0] = 1;
    v[642][0] = 1;
    v[643][0] = 1;
    v[644][0] = 1;
    v[645][0] = 1;
    v[646][0] = 1;
    v[647][0] = 1;
    v[648][0] = 1;
    v[649][0] = 1;
    v[650][0] = 1;
    v[651][0] = 1;
    v[652][0] = 1;
    v[653][0] = 1;
    v[654][0] = 1;
    v[655][0] = 1;
    v[656][0] = 1;
    v[657][0] = 1;
    v[658][0] = 1;
    v[659][0] = 1;
    v[660][0] = 1;
    v[661][0] = 1;
    v[662][0] = 1;
    v[663][0] = 1;
    v[664][0] = 1;
    v[665][0] = 1;
    v[666][0] = 1;
    v[667][0] = 1;
    v[668][0] = 1;
    v[669][0] = 1;
    v[670][0] = 1;
    v[671][0] = 1;
    v[672][0] = 1;
    v[673][0] = 1;
    v[674][0] = 1;
    v[675][0] = 1;
    v[676][0] = 1;
    v[677][0] = 1;
    v[678][0] = 1;
    v[679][0] = 1;
    v[680][0] = 1;
    v[681][0] = 1;
    v[682][0] = 1;
    v[683][0] = 1;
    v[684][0] = 1;
    v[685][0] = 1;
    v[686][0] = 1;
    v[687][0] = 1;
    v[688][0] = 1;
    v[689][0] = 1;
    v[690][0] = 1;
    v[691][0] = 1;
    v[692][0] = 1;
    v[693][0] = 1;
    v[694][0] = 1;
    v[695][0] = 1;
    v[696][0] = 1;
    v[697][0] = 1;
    v[698][0] = 1;
    v[699][0] = 1;
    v[700][0] = 1;
    v[701][0] = 1;
    v[702][0] = 1;
    v[703][0] = 1;
    v[704][0] = 1;
    v[705][0] = 1;
    v[706][0] = 1;
    v[707][0] = 1;
    v[708][0] = 1;
    v[709][0] = 1;
    v[710][0] = 1;
    v[711][0] = 1;
    v[712][0] = 1;
    v[713][0] = 1;
    v[714][0] = 1;
    v[715][0] = 1;
    v[716][0] = 1;
    v[717][0] = 1;
    v[718][0] = 1;
    v[719][0] = 1;
    v[720][0] = 1;
    v[721][0] = 1;
    v[722][0] = 1;
    v[723][0] = 1;
    v[724][0] = 1;
    v[725][0] = 1;
    v[726][0] = 1;
    v[727][0] = 1;
    v[728][0] = 1;
    v[729][0] = 1;
    v[730][0] = 1;
    v[731][0] = 1;
    v[732][0] = 1;
    v[733][0] = 1;
    v[734][0] = 1;
    v[735][0] = 1;
    v[736][0] = 1;
    v[737][0] = 1;
    v[738][0] = 1;
    v[739][0] = 1;
    v[740][0] = 1;
    v[741][0] = 1;
    v[742][0] = 1;
    v[743][0] = 1;
    v[744][0] = 1;
    v[745][0] = 1;
    v[746][0] = 1;
    v[747][0] = 1;
    v[748][0] = 1;
    v[749][0] = 1;
    v[750][0] = 1;
    v[751][0] = 1;
    v[752][0] = 1;
    v[753][0] = 1;
    v[754][0] = 1;
    v[755][0] = 1;
    v[756][0] = 1;
    v[757][0] = 1;
    v[758][0] = 1;
    v[759][0] = 1;
    v[760][0] = 1;
    v[761][0] = 1;
    v[762][0] = 1;
    v[763][0] = 1;
    v[764][0] = 1;
    v[765][0] = 1;
    v[766][0] = 1;
    v[767][0] = 1;
    v[768][0] = 1;
    v[769][0] = 1;
    v[770][0] = 1;
    v[771][0] = 1;
    v[772][0] = 1;
    v[773][0] = 1;
    v[774][0] = 1;
    v[775][0] = 1;
    v[776][0] = 1;
    v[777][0] = 1;
    v[778][0] = 1;
    v[779][0] = 1;
    v[780][0] = 1;
    v[781][0] = 1;
    v[782][0] = 1;
    v[783][0] = 1;
    v[784][0] = 1;
    v[785][0] = 1;
    v[786][0] = 1;
    v[787][0] = 1;
    v[788][0] = 1;
    v[789][0] = 1;
    v[790][0] = 1;
    v[791][0] = 1;
    v[792][0] = 1;
    v[793][0] = 1;
    v[794][0] = 1;
    v[795][0] = 1;
    v[796][0] = 1;
    v[797][0] = 1;
    v[798][0] = 1;
    v[799][0] = 1;
    v[800][0] = 1;
    v[801][0] = 1;
    v[802][0] = 1;
    v[803][0] = 1;
    v[804][0] = 1;
    v[805][0] = 1;
    v[806][0] = 1;
    v[807][0] = 1;
    v[808][0] = 1;
    v[809][0] = 1;
    v[810][0] = 1;
    v[811][0] = 1;
    v[812][0] = 1;
    v[813][0] = 1;
    v[814][0] = 1;
    v[815][0] = 1;
    v[816][0] = 1;
    v[817][0] = 1;
    v[818][0] = 1;
    v[819][0] = 1;
    v[820][0] = 1;
    v[821][0] = 1;
    v[822][0] = 1;
    v[823][0] = 1;
    v[824][0] = 1;
    v[825][0] = 1;
    v[826][0] = 1;
    v[827][0] = 1;
    v[828][0] = 1;
    v[829][0] = 1;
    v[830][0] = 1;
    v[831][0] = 1;
    v[832][0] = 1;
    v[833][0] = 1;
    v[834][0] = 1;
    v[835][0] = 1;
    v[836][0] = 1;
    v[837][0] = 1;
    v[838][0] = 1;
    v[839][0] = 1;
    v[840][0] = 1;
    v[841][0] = 1;
    v[842][0] = 1;
    v[843][0] = 1;
    v[844][0] = 1;
    v[845][0] = 1;
    v[846][0] = 1;
    v[847][0] = 1;
    v[848][0] = 1;
    v[849][0] = 1;
    v[850][0] = 1;
    v[851][0] = 1;
    v[852][0] = 1;
    v[853][0] = 1;
    v[854][0] = 1;
    v[855][0] = 1;
    v[856][0] = 1;
    v[857][0] = 1;
    v[858][0] = 1;
    v[859][0] = 1;
    v[860][0] = 1;
    v[861][0] = 1;
    v[862][0] = 1;
    v[863][0] = 1;
    v[864][0] = 1;
    v[865][0] = 1;
    v[866][0] = 1;
    v[867][0] = 1;
    v[868][0] = 1;
    v[869][0] = 1;
    v[870][0] = 1;
    v[871][0] = 1;
    v[872][0] = 1;
    v[873][0] = 1;
    v[874][0] = 1;
    v[875][0] = 1;
    v[876][0] = 1;
    v[877][0] = 1;
    v[878][0] = 1;
    v[879][0] = 1;
    v[880][0] = 1;
    v[881][0] = 1;
    v[882][0] = 1;
    v[883][0] = 1;
    v[884][0] = 1;
    v[885][0] = 1;
    v[886][0] = 1;
    v[887][0] = 1;
    v[888][0] = 1;
    v[889][0] = 1;
    v[890][0] = 1;
    v[891][0] = 1;
    v[892][0] = 1;
    v[893][0] = 1;
    v[894][0] = 1;
    v[895][0] = 1;
    v[896][0] = 1;
    v[897][0] = 1;
    v[898][0] = 1;
    v[899][0] = 1;
    v[900][0] = 1;
    v[901][0] = 1;
    v[902][0] = 1;
    v[903][0] = 1;
    v[904][0] = 1;
    v[905][0] = 1;
    v[906][0] = 1;
    v[907][0] = 1;
    v[908][0] = 1;
    v[909][0] = 1;
    v[910][0] = 1;
    v[911][0] = 1;
    v[912][0] = 1;
    v[913][0] = 1;
    v[914][0] = 1;
    v[915][0] = 1;
    v[916][0] = 1;
    v[917][0] = 1;
    v[918][0] = 1;
    v[919][0] = 1;
    v[920][0] = 1;
    v[921][0] = 1;
    v[922][0] = 1;
    v[923][0] = 1;
    v[924][0] = 1;
    v[925][0] = 1;
    v[926][0] = 1;
    v[927][0] = 1;
    v[928][0] = 1;
    v[929][0] = 1;
    v[930][0] = 1;
    v[931][0] = 1;
    v[932][0] = 1;
    v[933][0] = 1;
    v[934][0] = 1;
    v[935][0] = 1;
    v[936][0] = 1;
    v[937][0] = 1;
    v[938][0] = 1;
    v[939][0] = 1;
    v[940][0] = 1;
    v[941][0] = 1;
    v[942][0] = 1;
    v[943][0] = 1;
    v[944][0] = 1;
    v[945][0] = 1;
    v[946][0] = 1;
    v[947][0] = 1;
    v[948][0] = 1;
    v[949][0] = 1;
    v[950][0] = 1;
    v[951][0] = 1;
    v[952][0] = 1;
    v[953][0] = 1;
    v[954][0] = 1;
    v[955][0] = 1;
    v[956][0] = 1;
    v[957][0] = 1;
    v[958][0] = 1;
    v[959][0] = 1;
    v[960][0] = 1;
    v[961][0] = 1;
    v[962][0] = 1;
    v[963][0] = 1;
    v[964][0] = 1;
    v[965][0] = 1;
    v[966][0] = 1;
    v[967][0] = 1;
    v[968][0] = 1;
    v[969][0] = 1;
    v[970][0] = 1;
    v[971][0] = 1;
    v[972][0] = 1;
    v[973][0] = 1;
    v[974][0] = 1;
    v[975][0] = 1;
    v[976][0] = 1;
    v[977][0] = 1;
    v[978][0] = 1;
    v[979][0] = 1;
    v[980][0] = 1;
    v[981][0] = 1;
    v[982][0] = 1;
    v[983][0] = 1;
    v[984][0] = 1;
    v[985][0] = 1;
    v[986][0] = 1;
    v[987][0] = 1;
    v[988][0] = 1;
    v[989][0] = 1;
    v[990][0] = 1;
    v[991][0] = 1;
    v[992][0] = 1;
    v[993][0] = 1;
    v[994][0] = 1;
    v[995][0] = 1;
    v[996][0] = 1;
    v[997][0] = 1;
    v[998][0] = 1;
    v[999][0] = 1;
    v[1000][0] = 1;
    v[1001][0] = 1;
    v[1002][0] = 1;
    v[1003][0] = 1;
    v[1004][0] = 1;
    v[1005][0] = 1;
    v[1006][0] = 1;
    v[1007][0] = 1;
    v[1008][0] = 1;
    v[1009][0] = 1;
    v[1010][0] = 1;
    v[1011][0] = 1;
    v[1012][0] = 1;
    v[1013][0] = 1;
    v[1014][0] = 1;
    v[1015][0] = 1;
    v[1016][0] = 1;
    v[1017][0] = 1;
    v[1018][0] = 1;
    v[1019][0] = 1;
    v[1020][0] = 1;
    v[1021][0] = 1;
    v[1022][0] = 1;
    v[1023][0] = 1;
    v[1024][0] = 1;
    v[1025][0] = 1;
    v[1026][0] = 1;
    v[1027][0] = 1;
    v[1028][0] = 1;
    v[1029][0] = 1;
    v[1030][0] = 1;
    v[1031][0] = 1;
    v[1032][0] = 1;
    v[1033][0] = 1;
    v[1034][0] = 1;
    v[1035][0] = 1;
    v[1036][0] = 1;
    v[1037][0] = 1;
    v[1038][0] = 1;
    v[1039][0] = 1;
    v[1040][0] = 1;
    v[1041][0] = 1;
    v[1042][0] = 1;
    v[1043][0] = 1;
    v[1044][0] = 1;
    v[1045][0] = 1;
    v[1046][0] = 1;
    v[1047][0] = 1;
    v[1048][0] = 1;
    v[1049][0] = 1;
    v[1050][0] = 1;
    v[1051][0] = 1;
    v[1052][0] = 1;
    v[1053][0] = 1;
    v[1054][0] = 1;
    v[1055][0] = 1;
    v[1056][0] = 1;
    v[1057][0] = 1;
    v[1058][0] = 1;
    v[1059][0] = 1;
    v[1060][0] = 1;
    v[1061][0] = 1;
    v[1062][0] = 1;
    v[1063][0] = 1;
    v[1064][0] = 1;
    v[1065][0] = 1;
    v[1066][0] = 1;
    v[1067][0] = 1;
    v[1068][0] = 1;
    v[1069][0] = 1;
    v[1070][0] = 1;
    v[1071][0] = 1;
    v[1072][0] = 1;
    v[1073][0] = 1;
    v[1074][0] = 1;
    v[1075][0] = 1;
    v[1076][0] = 1;
    v[1077][0] = 1;
    v[1078][0] = 1;
    v[1079][0] = 1;
    v[1080][0] = 1;
    v[1081][0] = 1;
    v[1082][0] = 1;
    v[1083][0] = 1;
    v[1084][0] = 1;
    v[1085][0] = 1;
    v[1086][0] = 1;
    v[1087][0] = 1;
    v[1088][0] = 1;
    v[1089][0] = 1;
    v[1090][0] = 1;
    v[1091][0] = 1;
    v[1092][0] = 1;
    v[1093][0] = 1;
    v[1094][0] = 1;
    v[1095][0] = 1;
    v[1096][0] = 1;
    v[1097][0] = 1;
    v[1098][0] = 1;
    v[1099][0] = 1;
    v[1100][0] = 1;
    v[1101][0] = 1;
    v[1102][0] = 1;
    v[1103][0] = 1;
    v[1104][0] = 1;
    v[1105][0] = 1;
    v[1106][0] = 1;
    v[1107][0] = 1;
    v[1108][0] = 1;
    v[1109][0] = 1;
    v[1110][0] = 1;

    v[2][1] = 1;
    v[3][1] = 3;
    v[4][1] = 1;
    v[5][1] = 3;
    v[6][1] = 1;
    v[7][1] = 3;
    v[8][1] = 3;
    v[9][1] = 1;
    v[10][1] = 3;
    v[11][1] = 1;
    v[12][1] = 3;
    v[13][1] = 1;
    v[14][1] = 3;
    v[15][1] = 1;
    v[16][1] = 1;
    v[17][1] = 3;
    v[18][1] = 1;
    v[19][1] = 3;
    v[20][1] = 1;
    v[21][1] = 3;
    v[22][1] = 1;
    v[23][1] = 3;
    v[24][1] = 3;
    v[25][1] = 1;
    v[26][1] = 1;
    v[27][1] = 1;
    v[28][1] = 3;
    v[29][1] = 1;
    v[30][1] = 3;
    v[31][1] = 1;
    v[32][1] = 3;
    v[33][1] = 3;
    v[34][1] = 1;
    v[35][1] = 3;
    v[36][1] = 1;
    v[37][1] = 1;
    v[38][1] = 1;
    v[39][1] = 3;
    v[40][1] = 1;
    v[41][1] = 3;
    v[42][1] = 1;
    v[43][1] = 1;
    v[44][1] = 1;
    v[45][1] = 3;
    v[46][1] = 3;
    v[47][1] = 1;
    v[48][1] = 3;
    v[49][1] = 3;
    v[50][1] = 1;
    v[51][1] = 1;
    v[52][1] = 3;
    v[53][1] = 3;
    v[54][1] = 1;
    v[55][1] = 3;
    v[56][1] = 3;
    v[57][1] = 3;
    v[58][1] = 1;
    v[59][1] = 3;
    v[60][1] = 1;
    v[61][1] = 3;
    v[62][1] = 1;
    v[63][1] = 1;
    v[64][1] = 3;
    v[65][1] = 3;
    v[66][1] = 1;
    v[67][1] = 1;
    v[68][1] = 1;
    v[69][1] = 1;
    v[70][1] = 3;
    v[71][1] = 1;
    v[72][1] = 1;
    v[73][1] = 3;
    v[74][1] = 1;
    v[75][1] = 1;
    v[76][1] = 1;
    v[77][1] = 3;
    v[78][1] = 3;
    v[79][1] = 1;
    v[80][1] = 3;
    v[81][1] = 3;
    v[82][1] = 1;
    v[83][1] = 3;
    v[84][1] = 3;
    v[85][1] = 3;
    v[86][1] = 1;
    v[87][1] = 3;
    v[88][1] = 3;
    v[89][1] = 3;
    v[90][1] = 1;
    v[91][1] = 3;
    v[92][1] = 3;
    v[93][1] = 1;
    v[94][1] = 3;
    v[95][1] = 3;
    v[96][1] = 3;
    v[97][1] = 1;
    v[98][1] = 3;
    v[99][1] = 1;
    v[100][1] = 3;
    v[101][1] = 1;
    v[102][1] = 1;
    v[103][1] = 3;
    v[104][1] = 3;
    v[105][1] = 1;
    v[106][1] = 3;
    v[107][1] = 3;
    v[108][1] = 1;
    v[109][1] = 1;
    v[110][1] = 1;
    v[111][1] = 3;
    v[112][1] = 3;
    v[113][1] = 1;
    v[114][1] = 3;
    v[115][1] = 3;
    v[116][1] = 1;
    v[117][1] = 3;
    v[118][1] = 1;
    v[119][1] = 1;
    v[120][1] = 3;
    v[121][1] = 3;
    v[122][1] = 3;
    v[123][1] = 1;
    v[124][1] = 1;
    v[125][1] = 1;
    v[126][1] = 3;
    v[127][1] = 1;
    v[128][1] = 1;
    v[129][1] = 3;
    v[130][1] = 1;
    v[131][1] = 1;
    v[132][1] = 3;
    v[133][1] = 3;
    v[134][1] = 1;
    v[135][1] = 3;
    v[136][1] = 1;
    v[137][1] = 3;
    v[138][1] = 3;
    v[139][1] = 3;
    v[140][1] = 3;
    v[141][1] = 1;
    v[142][1] = 1;
    v[143][1] = 1;
    v[144][1] = 3;
    v[145][1] = 3;
    v[146][1] = 1;
    v[147][1] = 1;
    v[148][1] = 3;
    v[149][1] = 1;
    v[150][1] = 1;
    v[151][1] = 1;
    v[152][1] = 1;
    v[153][1] = 1;
    v[154][1] = 1;
    v[155][1] = 3;
    v[156][1] = 1;
    v[157][1] = 3;
    v[158][1] = 1;
    v[159][1] = 1;
    v[160][1] = 1;
    v[161][1] = 3;
    v[162][1] = 1;
    v[163][1] = 3;
    v[164][1] = 1;
    v[165][1] = 3;
    v[166][1] = 3;
    v[167][1] = 3;
    v[168][1] = 1;
    v[169][1] = 1;
    v[170][1] = 3;
    v[171][1] = 3;
    v[172][1] = 1;
    v[173][1] = 3;
    v[174][1] = 1;
    v[175][1] = 3;
    v[176][1] = 1;
    v[177][1] = 1;
    v[178][1] = 3;
    v[179][1] = 1;
    v[180][1] = 3;
    v[181][1] = 1;
    v[182][1] = 3;
    v[183][1] = 1;
    v[184][1] = 3;
    v[185][1] = 1;
    v[186][1] = 1;
    v[187][1] = 1;
    v[188][1] = 3;
    v[189][1] = 3;
    v[190][1] = 1;
    v[191][1] = 3;
    v[192][1] = 3;
    v[193][1] = 1;
    v[194][1] = 3;
    v[195][1] = 1;
    v[196][1] = 1;
    v[197][1] = 1;
    v[198][1] = 3;
    v[199][1] = 1;
    v[200][1] = 3;
    v[201][1] = 1;
    v[202][1] = 1;
    v[203][1] = 3;
    v[204][1] = 1;
    v[205][1] = 1;
    v[206][1] = 3;
    v[207][1] = 3;
    v[208][1] = 1;
    v[209][1] = 1;
    v[210][1] = 3;
    v[211][1] = 3;
    v[212][1] = 3;
    v[213][1] = 1;
    v[214][1] = 3;
    v[215][1] = 3;
    v[216][1] = 3;
    v[217][1] = 1;
    v[218][1] = 3;
    v[219][1] = 1;
    v[220][1] = 3;
    v[221][1] = 1;
    v[222][1] = 1;
    v[223][1] = 1;
    v[224][1] = 3;
    v[225][1] = 1;
    v[226][1] = 1;
    v[227][1] = 1;
    v[228][1] = 3;
    v[229][1] = 1;
    v[230][1] = 1;
    v[231][1] = 1;
    v[232][1] = 1;
    v[233][1] = 1;
    v[234][1] = 3;
    v[235][1] = 3;
    v[236][1] = 3;
    v[237][1] = 1;
    v[238][1] = 1;
    v[239][1] = 1;
    v[240][1] = 1;
    v[241][1] = 3;
    v[242][1] = 3;
    v[243][1] = 3;
    v[244][1] = 1;
    v[245][1] = 3;
    v[246][1] = 3;
    v[247][1] = 1;
    v[248][1] = 1;
    v[249][1] = 1;
    v[250][1] = 1;
    v[251][1] = 3;
    v[252][1] = 1;
    v[253][1] = 1;
    v[254][1] = 3;
    v[255][1] = 1;
    v[256][1] = 3;
    v[257][1] = 3;
    v[258][1] = 1;
    v[259][1] = 1;
    v[260][1] = 3;
    v[261][1] = 3;
    v[262][1] = 1;
    v[263][1] = 1;
    v[264][1] = 1;
    v[265][1] = 1;
    v[266][1] = 3;
    v[267][1] = 1;
    v[268][1] = 3;
    v[269][1] = 3;
    v[270][1] = 1;
    v[271][1] = 3;
    v[272][1] = 3;
    v[273][1] = 1;
    v[274][1] = 1;
    v[275][1] = 1;
    v[276][1] = 3;
    v[277][1] = 3;
    v[278][1] = 3;
    v[279][1] = 1;
    v[280][1] = 3;
    v[281][1] = 3;
    v[282][1] = 1;
    v[283][1] = 3;
    v[284][1] = 3;
    v[285][1] = 1;
    v[286][1] = 3;
    v[287][1] = 1;
    v[288][1] = 3;
    v[289][1] = 3;
    v[290][1] = 3;
    v[291][1] = 1;
    v[292][1] = 3;
    v[293][1] = 1;
    v[294][1] = 1;
    v[295][1] = 3;
    v[296][1] = 1;
    v[297][1] = 3;
    v[298][1] = 1;
    v[299][1] = 1;
    v[300][1] = 1;
    v[301][1] = 3;
    v[302][1] = 3;
    v[303][1] = 3;
    v[304][1] = 1;
    v[305][1] = 1;
    v[306][1] = 3;
    v[307][1] = 1;
    v[308][1] = 3;
    v[309][1] = 1;
    v[310][1] = 1;
    v[311][1] = 1;
    v[312][1] = 1;
    v[313][1] = 1;
    v[314][1] = 1;
    v[315][1] = 3;
    v[316][1] = 1;
    v[317][1] = 1;
    v[318][1] = 3;
    v[319][1] = 1;
    v[320][1] = 3;
    v[321][1] = 3;
    v[322][1] = 1;
    v[323][1] = 1;
    v[324][1] = 1;
    v[325][1] = 1;
    v[326][1] = 3;
    v[327][1] = 1;
    v[328][1] = 3;
    v[329][1] = 1;
    v[330][1] = 3;
    v[331][1] = 1;
    v[332][1] = 1;
    v[333][1] = 1;
    v[334][1] = 1;
    v[335][1] = 3;
    v[336][1] = 3;
    v[337][1] = 1;
    v[338][1] = 1;
    v[339][1] = 1;
    v[340][1] = 1;
    v[341][1] = 1;
    v[342][1] = 3;
    v[343][1] = 3;
    v[344][1] = 3;
    v[345][1] = 1;
    v[346][1] = 1;
    v[347][1] = 3;
    v[348][1] = 3;
    v[349][1] = 3;
    v[350][1] = 3;
    v[351][1] = 3;
    v[352][1] = 1;
    v[353][1] = 3;
    v[354][1] = 3;
    v[355][1] = 1;
    v[356][1] = 3;
    v[357][1] = 3;
    v[358][1] = 3;
    v[359][1] = 3;
    v[360][1] = 1;
    v[361][1] = 1;
    v[362][1] = 1;
    v[363][1] = 1;
    v[364][1] = 1;
    v[365][1] = 1;
    v[366][1] = 3;
    v[367][1] = 1;
    v[368][1] = 1;
    v[369][1] = 3;
    v[370][1] = 1;
    v[371][1] = 1;
    v[372][1] = 1;
    v[373][1] = 3;
    v[374][1] = 1;
    v[375][1] = 1;
    v[376][1] = 1;
    v[377][1] = 3;
    v[378][1] = 3;
    v[379][1] = 3;
    v[380][1] = 1;
    v[381][1] = 3;
    v[382][1] = 1;
    v[383][1] = 1;
    v[384][1] = 3;
    v[385][1] = 3;
    v[386][1] = 3;
    v[387][1] = 1;
    v[388][1] = 3;
    v[389][1] = 3;
    v[390][1] = 1;
    v[391][1] = 3;
    v[392][1] = 1;
    v[393][1] = 3;
    v[394][1] = 3;
    v[395][1] = 1;
    v[396][1] = 3;
    v[397][1] = 3;
    v[398][1] = 3;
    v[399][1] = 1;
    v[400][1] = 1;
    v[401][1] = 3;
    v[402][1] = 3;
    v[403][1] = 1;
    v[404][1] = 3;
    v[405][1] = 1;
    v[406][1] = 3;
    v[407][1] = 1;
    v[408][1] = 1;
    v[409][1] = 1;
    v[410][1] = 3;
    v[411][1] = 3;
    v[412][1] = 3;
    v[413][1] = 3;
    v[414][1] = 1;
    v[415][1] = 3;
    v[416][1] = 1;
    v[417][1] = 1;
    v[418][1] = 3;
    v[419][1] = 1;
    v[420][1] = 3;
    v[421][1] = 1;
    v[422][1] = 1;
    v[423][1] = 1;
    v[424][1] = 3;
    v[425][1] = 1;
    v[426][1] = 3;
    v[427][1] = 1;
    v[428][1] = 3;
    v[429][1] = 1;
    v[430][1] = 3;
    v[431][1] = 3;
    v[432][1] = 3;
    v[433][1] = 3;
    v[434][1] = 3;
    v[435][1] = 3;
    v[436][1] = 3;
    v[437][1] = 3;
    v[438][1] = 1;
    v[439][1] = 3;
    v[440][1] = 3;
    v[441][1] = 3;
    v[442][1] = 3;
    v[443][1] = 3;
    v[444][1] = 1;
    v[445][1] = 3;
    v[446][1] = 1;
    v[447][1] = 3;
    v[448][1] = 3;
    v[449][1] = 3;
    v[450][1] = 1;
    v[451][1] = 3;
    v[452][1] = 1;
    v[453][1] = 3;
    v[454][1] = 1;
    v[455][1] = 3;
    v[456][1] = 3;
    v[457][1] = 1;
    v[458][1] = 3;
    v[459][1] = 3;
    v[460][1] = 3;
    v[461][1] = 3;
    v[462][1] = 3;
    v[463][1] = 3;
    v[464][1] = 3;
    v[465][1] = 3;
    v[466][1] = 3;
    v[467][1] = 1;
    v[468][1] = 1;
    v[469][1] = 1;
    v[470][1] = 1;
    v[471][1] = 1;
    v[472][1] = 1;
    v[473][1] = 3;
    v[474][1] = 3;
    v[475][1] = 1;
    v[476][1] = 1;
    v[477][1] = 3;
    v[478][1] = 3;
    v[479][1] = 1;
    v[480][1] = 1;
    v[481][1] = 1;
    v[482][1] = 3;
    v[483][1] = 3;
    v[484][1] = 1;
    v[485][1] = 1;
    v[486][1] = 3;
    v[487][1] = 3;
    v[488][1] = 3;
    v[489][1] = 3;
    v[490][1] = 1;
    v[491][1] = 1;
    v[492][1] = 3;
    v[493][1] = 1;
    v[494][1] = 3;
    v[495][1] = 3;
    v[496][1] = 1;
    v[497][1] = 3;
    v[498][1] = 3;
    v[499][1] = 1;
    v[500][1] = 1;
    v[501][1] = 1;
    v[502][1] = 3;
    v[503][1] = 3;
    v[504][1] = 3;
    v[505][1] = 1;
    v[506][1] = 1;
    v[507][1] = 3;
    v[508][1] = 3;
    v[509][1] = 3;
    v[510][1] = 3;
    v[511][1] = 3;
    v[512][1] = 1;
    v[513][1] = 1;
    v[514][1] = 1;
    v[515][1] = 3;
    v[516][1] = 1;
    v[517][1] = 3;
    v[518][1] = 3;
    v[519][1] = 1;
    v[520][1] = 3;
    v[521][1] = 3;
    v[522][1] = 3;
    v[523][1] = 3;
    v[524][1] = 1;
    v[525][1] = 1;
    v[526][1] = 3;
    v[527][1] = 1;
    v[528][1] = 1;
    v[529][1] = 3;
    v[530][1] = 1;
    v[531][1] = 3;
    v[532][1] = 1;
    v[533][1] = 3;
    v[534][1] = 1;
    v[535][1] = 3;
    v[536][1] = 3;
    v[537][1] = 1;
    v[538][1] = 1;
    v[539][1] = 3;
    v[540][1] = 3;
    v[541][1] = 1;
    v[542][1] = 3;
    v[543][1] = 3;
    v[544][1] = 1;
    v[545][1] = 3;
    v[546][1] = 3;
    v[547][1] = 1;
    v[548][1] = 1;
    v[549][1] = 3;
    v[550][1] = 1;
    v[551][1] = 3;
    v[552][1] = 3;
    v[553][1] = 1;
    v[554][1] = 1;
    v[555][1] = 3;
    v[556][1] = 1;
    v[557][1] = 3;
    v[558][1] = 1;
    v[559][1] = 3;
    v[560][1] = 1;
    v[561][1] = 1;
    v[562][1] = 3;
    v[563][1] = 3;
    v[564][1] = 1;
    v[565][1] = 1;
    v[566][1] = 1;
    v[567][1] = 3;
    v[568][1] = 3;
    v[569][1] = 1;
    v[570][1] = 3;
    v[571][1] = 1;
    v[572][1] = 1;
    v[573][1] = 3;
    v[574][1] = 3;
    v[575][1] = 1;
    v[576][1] = 1;
    v[577][1] = 3;
    v[578][1] = 1;
    v[579][1] = 3;
    v[580][1] = 1;
    v[581][1] = 1;
    v[582][1] = 1;
    v[583][1] = 1;
    v[584][1] = 1;
    v[585][1] = 3;
    v[586][1] = 1;
    v[587][1] = 1;
    v[588][1] = 1;
    v[589][1] = 1;
    v[590][1] = 3;
    v[591][1] = 1;
    v[592][1] = 3;
    v[593][1] = 1;
    v[594][1] = 1;
    v[595][1] = 3;
    v[596][1] = 3;
    v[597][1] = 1;
    v[598][1] = 1;
    v[599][1] = 3;
    v[600][1] = 1;
    v[601][1] = 3;
    v[602][1] = 1;
    v[603][1] = 3;
    v[604][1] = 3;
    v[605][1] = 3;
    v[606][1] = 1;
    v[607][1] = 3;
    v[608][1] = 3;
    v[609][1] = 3;
    v[610][1] = 1;
    v[611][1] = 1;
    v[612][1] = 3;
    v[613][1] = 3;
    v[614][1] = 3;
    v[615][1] = 1;
    v[616][1] = 1;
    v[617][1] = 1;
    v[618][1] = 1;
    v[619][1] = 3;
    v[620][1] = 1;
    v[621][1] = 3;
    v[622][1] = 1;
    v[623][1] = 3;
    v[624][1] = 1;
    v[625][1] = 1;
    v[626][1] = 3;
    v[627][1] = 3;
    v[628][1] = 1;
    v[629][1] = 1;
    v[630][1] = 1;
    v[631][1] = 3;
    v[632][1] = 3;
    v[633][1] = 1;
    v[634][1] = 3;
    v[635][1] = 1;
    v[636][1] = 3;
    v[637][1] = 1;
    v[638][1] = 1;
    v[639][1] = 1;
    v[640][1] = 1;
    v[641][1] = 1;
    v[642][1] = 1;
    v[643][1] = 3;
    v[644][1] = 1;
    v[645][1] = 3;
    v[646][1] = 3;
    v[647][1] = 1;
    v[648][1] = 3;
    v[649][1] = 3;
    v[650][1] = 3;
    v[651][1] = 1;
    v[652][1] = 3;
    v[653][1] = 1;
    v[654][1] = 1;
    v[655][1] = 3;
    v[656][1] = 3;
    v[657][1] = 1;
    v[658][1] = 1;
    v[659][1] = 3;
    v[660][1] = 3;
    v[661][1] = 1;
    v[662][1] = 1;
    v[663][1] = 1;
    v[664][1] = 3;
    v[665][1] = 1;
    v[666][1] = 3;
    v[667][1] = 3;
    v[668][1] = 1;
    v[669][1] = 1;
    v[670][1] = 3;
    v[671][1] = 1;
    v[672][1] = 1;
    v[673][1] = 3;
    v[674][1] = 1;
    v[675][1] = 3;
    v[676][1] = 1;
    v[677][1] = 1;
    v[678][1] = 1;
    v[679][1] = 3;
    v[680][1] = 3;
    v[681][1] = 3;
    v[682][1] = 3;
    v[683][1] = 1;
    v[684][1] = 1;
    v[685][1] = 3;
    v[686][1] = 3;
    v[687][1] = 1;
    v[688][1] = 1;
    v[689][1] = 1;
    v[690][1] = 1;
    v[691][1] = 3;
    v[692][1] = 1;
    v[693][1] = 1;
    v[694][1] = 3;
    v[695][1] = 3;
    v[696][1] = 3;
    v[697][1] = 1;
    v[698][1] = 1;
    v[699][1] = 3;
    v[700][1] = 3;
    v[701][1] = 1;
    v[702][1] = 3;
    v[703][1] = 3;
    v[704][1] = 1;
    v[705][1] = 1;
    v[706][1] = 3;
    v[707][1] = 3;
    v[708][1] = 3;
    v[709][1] = 3;
    v[710][1] = 3;
    v[711][1] = 3;
    v[712][1] = 3;
    v[713][1] = 1;
    v[714][1] = 3;
    v[715][1] = 3;
    v[716][1] = 1;
    v[717][1] = 3;
    v[718][1] = 1;
    v[719][1] = 3;
    v[720][1] = 1;
    v[721][1] = 1;
    v[722][1] = 3;
    v[723][1] = 3;
    v[724][1] = 1;
    v[725][1] = 1;
    v[726][1] = 1;
    v[727][1] = 3;
    v[728][1] = 1;
    v[729][1] = 3;
    v[730][1] = 3;
    v[731][1] = 1;
    v[732][1] = 3;
    v[733][1] = 3;
    v[734][1] = 1;
    v[735][1] = 3;
    v[736][1] = 1;
    v[737][1] = 1;
    v[738][1] = 3;
    v[739][1] = 3;
    v[740][1] = 3;
    v[741][1] = 1;
    v[742][1] = 1;
    v[743][1] = 1;
    v[744][1] = 3;
    v[745][1] = 1;
    v[746][1] = 1;
    v[747][1] = 1;
    v[748][1] = 3;
    v[749][1] = 3;
    v[750][1] = 3;
    v[751][1] = 1;
    v[752][1] = 3;
    v[753][1] = 3;
    v[754][1] = 1;
    v[755][1] = 3;
    v[756][1] = 1;
    v[757][1] = 1;
    v[758][1] = 3;
    v[759][1] = 3;
    v[760][1] = 3;
    v[761][1] = 1;
    v[762][1] = 3;
    v[763][1] = 3;
    v[764][1] = 1;
    v[765][1] = 1;
    v[766][1] = 1;
    v[767][1] = 3;
    v[768][1] = 1;
    v[769][1] = 3;
    v[770][1] = 3;
    v[771][1] = 3;
    v[772][1] = 3;
    v[773][1] = 3;
    v[774][1] = 3;
    v[775][1] = 3;
    v[776][1] = 3;
    v[777][1] = 1;
    v[778][1] = 3;
    v[779][1] = 3;
    v[780][1] = 1;
    v[781][1] = 3;
    v[782][1] = 1;
    v[783][1] = 1;
    v[784][1] = 3;
    v[785][1] = 3;
    v[786][1] = 3;
    v[787][1] = 1;
    v[788][1] = 3;
    v[789][1] = 3;
    v[790][1] = 3;
    v[791][1] = 3;
    v[792][1] = 3;
    v[793][1] = 1;
    v[794][1] = 3;
    v[795][1] = 3;
    v[796][1] = 3;
    v[797][1] = 1;
    v[798][1] = 1;
    v[799][1] = 1;
    v[800][1] = 3;
    v[801][1] = 3;
    v[802][1] = 1;
    v[803][1] = 3;
    v[804][1] = 3;
    v[805][1] = 1;
    v[806][1] = 3;
    v[807][1] = 1;
    v[808][1] = 3;
    v[809][1] = 1;
    v[810][1] = 3;
    v[811][1] = 1;
    v[812][1] = 3;
    v[813][1] = 3;
    v[814][1] = 3;
    v[815][1] = 3;
    v[816][1] = 3;
    v[817][1] = 3;
    v[818][1] = 1;
    v[819][1] = 1;
    v[820][1] = 3;
    v[821][1] = 1;
    v[822][1] = 3;
    v[823][1] = 1;
    v[824][1] = 1;
    v[825][1] = 1;
    v[826][1] = 1;
    v[827][1] = 1;
    v[828][1] = 3;
    v[829][1] = 1;
    v[830][1] = 1;
    v[831][1] = 1;
    v[832][1] = 3;
    v[833][1] = 1;
    v[834][1] = 3;
    v[835][1] = 1;
    v[836][1] = 1;
    v[837][1] = 3;
    v[838][1] = 3;
    v[839][1] = 3;
    v[840][1] = 1;
    v[841][1] = 3;
    v[842][1] = 1;
    v[843][1] = 3;
    v[844][1] = 1;
    v[845][1] = 1;
    v[846][1] = 3;
    v[847][1] = 1;
    v[848][1] = 3;
    v[849][1] = 3;
    v[850][1] = 1;
    v[851][1] = 3;
    v[852][1] = 1;
    v[853][1] = 3;
    v[854][1] = 3;
    v[855][1] = 1;
    v[856][1] = 3;
    v[857][1] = 3;
    v[858][1] = 1;
    v[859][1] = 3;
    v[860][1] = 3;
    v[861][1] = 3;
    v[862][1] = 3;
    v[863][1] = 3;
    v[864][1] = 3;
    v[865][1] = 1;
    v[866][1] = 3;
    v[867][1] = 1;
    v[868][1] = 1;
    v[869][1] = 3;
    v[870][1] = 3;
    v[871][1] = 3;
    v[872][1] = 1;
    v[873][1] = 1;
    v[874][1] = 3;
    v[875][1] = 3;
    v[876][1] = 3;
    v[877][1] = 3;
    v[878][1] = 3;
    v[879][1] = 3;
    v[880][1] = 3;
    v[881][1] = 1;
    v[882][1] = 3;
    v[883][1] = 3;
    v[884][1] = 3;
    v[885][1] = 3;
    v[886][1] = 1;
    v[887][1] = 3;
    v[888][1] = 1;
    v[889][1] = 3;
    v[890][1] = 3;
    v[891][1] = 3;
    v[892][1] = 1;
    v[893][1] = 3;
    v[894][1] = 1;
    v[895][1] = 3;
    v[896][1] = 1;
    v[897][1] = 1;
    v[898][1] = 1;
    v[899][1] = 3;
    v[900][1] = 3;
    v[901][1] = 1;
    v[902][1] = 3;
    v[903][1] = 1;
    v[904][1] = 1;
    v[905][1] = 3;
    v[906][1] = 3;
    v[907][1] = 1;
    v[908][1] = 3;
    v[909][1] = 1;
    v[910][1] = 1;
    v[911][1] = 1;
    v[912][1] = 1;
    v[913][1] = 3;
    v[914][1] = 1;
    v[915][1] = 3;
    v[916][1] = 1;
    v[917][1] = 1;
    v[918][1] = 3;
    v[919][1] = 1;
    v[920][1] = 3;
    v[921][1] = 1;
    v[922][1] = 3;
    v[923][1] = 3;
    v[924][1] = 3;
    v[925][1] = 3;
    v[926][1] = 3;
    v[927][1] = 3;
    v[928][1] = 1;
    v[929][1] = 3;
    v[930][1] = 3;
    v[931][1] = 3;
    v[932][1] = 3;
    v[933][1] = 1;
    v[934][1] = 3;
    v[935][1] = 3;
    v[936][1] = 1;
    v[937][1] = 3;
    v[938][1] = 3;
    v[939][1] = 3;
    v[940][1] = 3;
    v[941][1] = 3;
    v[942][1] = 1;
    v[943][1] = 1;
    v[944][1] = 1;
    v[945][1] = 1;
    v[946][1] = 3;
    v[947][1] = 3;
    v[948][1] = 3;
    v[949][1] = 1;
    v[950][1] = 3;
    v[951][1] = 3;
    v[952][1] = 1;
    v[953][1] = 1;
    v[954][1] = 3;
    v[955][1] = 3;
    v[956][1] = 1;
    v[957][1] = 1;
    v[958][1] = 3;
    v[959][1] = 3;
    v[960][1] = 1;
    v[961][1] = 3;
    v[962][1] = 1;
    v[963][1] = 1;
    v[964][1] = 3;
    v[965][1] = 1;
    v[966][1] = 3;
    v[967][1] = 3;
    v[968][1] = 3;
    v[969][1] = 3;
    v[970][1] = 3;
    v[971][1] = 1;
    v[972][1] = 3;
    v[973][1] = 1;
    v[974][1] = 1;
    v[975][1] = 3;
    v[976][1] = 3;
    v[977][1] = 3;
    v[978][1] = 3;
    v[979][1] = 1;
    v[980][1] = 3;
    v[981][1] = 1;
    v[982][1] = 1;
    v[983][1] = 3;
    v[984][1] = 3;
    v[985][1] = 3;
    v[986][1] = 3;
    v[987][1] = 3;
    v[988][1] = 3;
    v[989][1] = 1;
    v[990][1] = 1;
    v[991][1] = 3;
    v[992][1] = 1;
    v[993][1] = 3;
    v[994][1] = 1;
    v[995][1] = 1;
    v[996][1] = 3;
    v[997][1] = 1;
    v[998][1] = 1;
    v[999][1] = 1;
    v[1000][1] = 1;
    v[1001][1] = 3;
    v[1002][1] = 3;
    v[1003][1] = 1;
    v[1004][1] = 1;
    v[1005][1] = 3;
    v[1006][1] = 1;
    v[1007][1] = 1;
    v[1008][1] = 1;
    v[1009][1] = 3;
    v[1010][1] = 1;
    v[1011][1] = 3;
    v[1012][1] = 1;
    v[1013][1] = 1;
    v[1014][1] = 3;
    v[1015][1] = 3;
    v[1016][1] = 1;
    v[1017][1] = 3;
    v[1018][1] = 1;
    v[1019][1] = 1;
    v[1020][1] = 3;
    v[1021][1] = 3;
    v[1022][1] = 3;
    v[1023][1] = 3;
    v[1024][1] = 3;
    v[1025][1] = 1;
    v[1026][1] = 3;
    v[1027][1] = 1;
    v[1028][1] = 1;
    v[1029][1] = 1;
    v[1030][1] = 3;
    v[1031][1] = 1;
    v[1032][1] = 1;
    v[1033][1] = 1;
    v[1034][1] = 3;
    v[1035][1] = 1;
    v[1036][1] = 1;
    v[1037][1] = 3;
    v[1038][1] = 1;
    v[1039][1] = 3;
    v[1040][1] = 3;
    v[1041][1] = 3;
    v[1042][1] = 3;
    v[1043][1] = 3;
    v[1044][1] = 1;
    v[1045][1] = 1;
    v[1046][1] = 1;
    v[1047][1] = 3;
    v[1048][1] = 3;
    v[1049][1] = 3;
    v[1050][1] = 3;
    v[1051][1] = 1;
    v[1052][1] = 3;
    v[1053][1] = 3;
    v[1054][1] = 3;
    v[1055][1] = 3;
    v[1056][1] = 1;
    v[1057][1] = 1;
    v[1058][1] = 3;
    v[1059][1] = 3;
    v[1060][1] = 3;
    v[1061][1] = 1;
    v[1062][1] = 3;
    v[1063][1] = 1;
    v[1064][1] = 1;
    v[1065][1] = 3;
    v[1066][1] = 3;
    v[1067][1] = 1;
    v[1068][1] = 3;
    v[1069][1] = 3;
    v[1070][1] = 1;
    v[1071][1] = 1;
    v[1072][1] = 1;
    v[1073][1] = 1;
    v[1074][1] = 1;
    v[1075][1] = 3;
    v[1076][1] = 1;
    v[1077][1] = 1;
    v[1078][1] = 3;
    v[1079][1] = 3;
    v[1080][1] = 1;
    v[1081][1] = 1;
    v[1082][1] = 1;
    v[1083][1] = 3;
    v[1084][1] = 1;
    v[1085][1] = 1;
    v[1086][1] = 3;
    v[1087][1] = 3;
    v[1088][1] = 1;
    v[1089][1] = 3;
    v[1090][1] = 3;
    v[1091][1] = 3;
    v[1092][1] = 3;
    v[1093][1] = 3;
    v[1094][1] = 3;
    v[1095][1] = 3;
    v[1096][1] = 3;
    v[1097][1] = 1;
    v[1098][1] = 1;
    v[1099][1] = 3;
    v[1100][1] = 3;
    v[1101][1] = 1;
    v[1102][1] = 1;
    v[1103][1] = 3;
    v[1104][1] = 1;
    v[1105][1] = 3;
    v[1106][1] = 3;
    v[1107][1] = 3;
    v[1108][1] = 3;
    v[1109][1] = 3;
    v[1110][1] = 1;

    v[3][2] = 7;
    v[4][2] = 5;
    v[5][2] = 1;
    v[6][2] = 3;
    v[7][2] = 3;
    v[8][2] = 7;
    v[9][2] = 5;
    v[10][2] = 5;
    v[11][2] = 7;
    v[12][2] = 7;
    v[13][2] = 1;
    v[14][2] = 3;
    v[15][2] = 3;
    v[16][2] = 7;
    v[17][2] = 5;
    v[18][2] = 1;
    v[19][2] = 1;
    v[20][2] = 5;
    v[21][2] = 3;
    v[22][2] = 7;
    v[23][2] = 1;
    v[24][2] = 7;
    v[25][2] = 5;
    v[26][2] = 1;
    v[27][2] = 3;
    v[28][2] = 7;
    v[29][2] = 7;
    v[30][2] = 1;
    v[31][2] = 1;
    v[32][2] = 1;
    v[33][2] = 5;
    v[34][2] = 7;
    v[35][2] = 7;
    v[36][2] = 5;
    v[37][2] = 1;
    v[38][2] = 3;
    v[39][2] = 3;
    v[40][2] = 7;
    v[41][2] = 5;
    v[42][2] = 5;
    v[43][2] = 5;
    v[44][2] = 3;
    v[45][2] = 3;
    v[46][2] = 3;
    v[47][2] = 1;
    v[48][2] = 1;
    v[49][2] = 5;
    v[50][2] = 1;
    v[51][2] = 1;
    v[52][2] = 5;
    v[53][2] = 3;
    v[54][2] = 3;
    v[55][2] = 3;
    v[56][2] = 3;
    v[57][2] = 1;
    v[58][2] = 3;
    v[59][2] = 7;
    v[60][2] = 5;
    v[61][2] = 7;
    v[62][2] = 3;
    v[63][2] = 7;
    v[64][2] = 1;
    v[65][2] = 3;
    v[66][2] = 3;
    v[67][2] = 5;
    v[68][2] = 1;
    v[69][2] = 3;
    v[70][2] = 5;
    v[71][2] = 5;
    v[72][2] = 7;
    v[73][2] = 7;
    v[74][2] = 7;
    v[75][2] = 1;
    v[76][2] = 1;
    v[77][2] = 3;
    v[78][2] = 3;
    v[79][2] = 1;
    v[80][2] = 1;
    v[81][2] = 5;
    v[82][2] = 1;
    v[83][2] = 5;
    v[84][2] = 7;
    v[85][2] = 5;
    v[86][2] = 1;
    v[87][2] = 7;
    v[88][2] = 5;
    v[89][2] = 3;
    v[90][2] = 3;
    v[91][2] = 1;
    v[92][2] = 5;
    v[93][2] = 7;
    v[94][2] = 1;
    v[95][2] = 7;
    v[96][2] = 5;
    v[97][2] = 1;
    v[98][2] = 7;
    v[99][2] = 3;
    v[100][2] = 1;
    v[101][2] = 7;
    v[102][2] = 1;
    v[103][2] = 7;
    v[104][2] = 3;
    v[105][2] = 3;
    v[106][2] = 5;
    v[107][2] = 7;
    v[108][2] = 3;
    v[109][2] = 3;
    v[110][2] = 5;
    v[111][2] = 1;
    v[112][2] = 3;
    v[113][2] = 3;
    v[114][2] = 1;
    v[115][2] = 3;
    v[116][2] = 5;
    v[117][2] = 1;
    v[118][2] = 3;
    v[119][2] = 3;
    v[120][2] = 3;
    v[121][2] = 7;
    v[122][2] = 1;
    v[123][2] = 1;
    v[124][2] = 7;
    v[125][2] = 3;
    v[126][2] = 1;
    v[127][2] = 3;
    v[128][2] = 7;
    v[129][2] = 5;
    v[130][2] = 5;
    v[131][2] = 7;
    v[132][2] = 5;
    v[133][2] = 5;
    v[134][2] = 3;
    v[135][2] = 1;
    v[136][2] = 3;
    v[137][2] = 3;
    v[138][2] = 3;
    v[139][2] = 1;
    v[140][2] = 3;
    v[141][2] = 3;
    v[142][2] = 7;
    v[143][2] = 3;
    v[144][2] = 3;
    v[145][2] = 1;
    v[146][2] = 7;
    v[147][2] = 5;
    v[148][2] = 1;
    v[149][2] = 7;
    v[150][2] = 7;
    v[151][2] = 5;
    v[152][2] = 7;
    v[153][2] = 5;
    v[154][2] = 1;
    v[155][2] = 3;
    v[156][2] = 1;
    v[157][2] = 7;
    v[158][2] = 3;
    v[159][2] = 7;
    v[160][2] = 3;
    v[161][2] = 5;
    v[162][2] = 7;
    v[163][2] = 3;
    v[164][2] = 1;
    v[165][2] = 3;
    v[166][2] = 3;
    v[167][2] = 3;
    v[168][2] = 1;
    v[169][2] = 5;
    v[170][2] = 7;
    v[171][2] = 3;
    v[172][2] = 3;
    v[173][2] = 7;
    v[174][2] = 7;
    v[175][2] = 7;
    v[176][2] = 5;
    v[177][2] = 3;
    v[178][2] = 1;
    v[179][2] = 7;
    v[180][2] = 1;
    v[181][2] = 3;
    v[182][2] = 7;
    v[183][2] = 5;
    v[184][2] = 3;
    v[185][2] = 3;
    v[186][2] = 3;
    v[187][2] = 7;
    v[188][2] = 1;
    v[189][2] = 1;
    v[190][2] = 3;
    v[191][2] = 1;
    v[192][2] = 5;
    v[193][2] = 7;
    v[194][2] = 1;
    v[195][2] = 3;
    v[196][2] = 5;
    v[197][2] = 3;
    v[198][2] = 5;
    v[199][2] = 3;
    v[200][2] = 3;
    v[201][2] = 7;
    v[202][2] = 5;
    v[203][2] = 5;
    v[204][2] = 3;
    v[205][2] = 3;
    v[206][2] = 1;
    v[207][2] = 3;
    v[208][2] = 7;
    v[209][2] = 7;
    v[210][2] = 7;
    v[211][2] = 1;
    v[212][2] = 5;
    v[213][2] = 7;
    v[214][2] = 1;
    v[215][2] = 3;
    v[216][2] = 1;
    v[217][2] = 1;
    v[218][2] = 7;
    v[219][2] = 1;
    v[220][2] = 3;
    v[221][2] = 1;
    v[222][2] = 7;
    v[223][2] = 1;
    v[224][2] = 5;
    v[225][2] = 3;
    v[226][2] = 5;
    v[227][2] = 3;
    v[228][2] = 1;
    v[229][2] = 1;
    v[230][2] = 5;
    v[231][2] = 5;
    v[232][2] = 3;
    v[233][2] = 3;
    v[234][2] = 5;
    v[235][2] = 7;
    v[236][2] = 1;
    v[237][2] = 5;
    v[238][2] = 3;
    v[239][2] = 7;
    v[240][2] = 7;
    v[241][2] = 3;
    v[242][2] = 5;
    v[243][2] = 3;
    v[244][2] = 3;
    v[245][2] = 1;
    v[246][2] = 7;
    v[247][2] = 3;
    v[248][2] = 1;
    v[249][2] = 3;
    v[250][2] = 5;
    v[251][2] = 7;
    v[252][2] = 1;
    v[253][2] = 3;
    v[254][2] = 7;
    v[255][2] = 1;
    v[256][2] = 5;
    v[257][2] = 1;
    v[258][2] = 3;
    v[259][2] = 1;
    v[260][2] = 5;
    v[261][2] = 3;
    v[262][2] = 1;
    v[263][2] = 7;
    v[264][2] = 1;
    v[265][2] = 5;
    v[266][2] = 5;
    v[267][2] = 5;
    v[268][2] = 3;
    v[269][2] = 7;
    v[270][2] = 1;
    v[271][2] = 1;
    v[272][2] = 7;
    v[273][2] = 3;
    v[274][2] = 1;
    v[275][2] = 1;
    v[276][2] = 7;
    v[277][2] = 5;
    v[278][2] = 7;
    v[279][2] = 5;
    v[280][2] = 7;
    v[281][2] = 7;
    v[282][2] = 3;
    v[283][2] = 7;
    v[284][2] = 1;
    v[285][2] = 3;
    v[286][2] = 7;
    v[287][2] = 7;
    v[288][2] = 3;
    v[289][2] = 5;
    v[290][2] = 1;
    v[291][2] = 1;
    v[292][2] = 7;
    v[293][2] = 1;
    v[294][2] = 5;
    v[295][2] = 5;
    v[296][2] = 5;
    v[297][2] = 1;
    v[298][2] = 5;
    v[299][2] = 1;
    v[300][2] = 7;
    v[301][2] = 5;
    v[302][2] = 5;
    v[303][2] = 7;
    v[304][2] = 1;
    v[305][2] = 1;
    v[306][2] = 7;
    v[307][2] = 1;
    v[308][2] = 7;
    v[309][2] = 7;
    v[310][2] = 1;
    v[311][2] = 1;
    v[312][2] = 3;
    v[313][2] = 3;
    v[314][2] = 3;
    v[315][2] = 7;
    v[316][2] = 7;
    v[317][2] = 5;
    v[318][2] = 3;
    v[319][2] = 7;
    v[320][2] = 3;
    v[321][2] = 1;
    v[322][2] = 3;
    v[323][2] = 7;
    v[324][2] = 5;
    v[325][2] = 3;
    v[326][2] = 3;
    v[327][2] = 5;
    v[328][2] = 7;
    v[329][2] = 1;
    v[330][2] = 1;
    v[331][2] = 5;
    v[332][2] = 5;
    v[333][2] = 7;
    v[334][2] = 7;
    v[335][2] = 1;
    v[336][2] = 1;
    v[337][2] = 1;
    v[338][2] = 1;
    v[339][2] = 5;
    v[340][2] = 5;
    v[341][2] = 5;
    v[342][2] = 7;
    v[343][2] = 5;
    v[344][2] = 7;
    v[345][2] = 1;
    v[346][2] = 1;
    v[347][2] = 3;
    v[348][2] = 5;
    v[349][2] = 1;
    v[350][2] = 3;
    v[351][2] = 3;
    v[352][2] = 7;
    v[353][2] = 3;
    v[354][2] = 7;
    v[355][2] = 5;
    v[356][2] = 3;
    v[357][2] = 5;
    v[358][2] = 3;
    v[359][2] = 1;
    v[360][2] = 7;
    v[361][2] = 1;
    v[362][2] = 7;
    v[363][2] = 7;
    v[364][2] = 1;
    v[365][2] = 1;
    v[366][2] = 7;
    v[367][2] = 7;
    v[368][2] = 7;
    v[369][2] = 5;
    v[370][2] = 5;
    v[371][2] = 1;
    v[372][2] = 1;
    v[373][2] = 7;
    v[374][2] = 5;
    v[375][2] = 5;
    v[376][2] = 7;
    v[377][2] = 5;
    v[378][2] = 1;
    v[379][2] = 1;
    v[380][2] = 5;
    v[381][2] = 5;
    v[382][2] = 5;
    v[383][2] = 5;
    v[384][2] = 5;
    v[385][2] = 5;
    v[386][2] = 1;
    v[387][2] = 3;
    v[388][2] = 1;
    v[389][2] = 5;
    v[390][2] = 7;
    v[391][2] = 3;
    v[392][2] = 3;
    v[393][2] = 5;
    v[394][2] = 7;
    v[395][2] = 3;
    v[396][2] = 7;
    v[397][2] = 1;
    v[398][2] = 7;
    v[399][2] = 7;
    v[400][2] = 1;
    v[401][2] = 3;
    v[402][2] = 5;
    v[403][2] = 1;
    v[404][2] = 5;
    v[405][2] = 5;
    v[406][2] = 3;
    v[407][2] = 7;
    v[408][2] = 3;
    v[409][2] = 7;
    v[410][2] = 7;
    v[411][2] = 5;
    v[412][2] = 7;
    v[413][2] = 5;
    v[414][2] = 7;
    v[415][2] = 1;
    v[416][2] = 1;
    v[417][2] = 5;
    v[418][2] = 3;
    v[419][2] = 5;
    v[420][2] = 1;
    v[421][2] = 5;
    v[422][2] = 3;
    v[423][2] = 7;
    v[424][2] = 1;
    v[425][2] = 5;
    v[426][2] = 7;
    v[427][2] = 7;
    v[428][2] = 3;
    v[429][2] = 5;
    v[430][2] = 1;
    v[431][2] = 3;
    v[432][2] = 5;
    v[433][2] = 1;
    v[434][2] = 5;
    v[435][2] = 3;
    v[436][2] = 3;
    v[437][2] = 3;
    v[438][2] = 7;
    v[439][2] = 3;
    v[440][2] = 5;
    v[441][2] = 1;
    v[442][2] = 3;
    v[443][2] = 7;
    v[444][2] = 7;
    v[445][2] = 3;
    v[446][2] = 7;
    v[447][2] = 5;
    v[448][2] = 3;
    v[449][2] = 3;
    v[450][2] = 1;
    v[451][2] = 7;
    v[452][2] = 5;
    v[453][2] = 1;
    v[454][2] = 1;
    v[455][2] = 3;
    v[456][2] = 7;
    v[457][2] = 1;
    v[458][2] = 7;
    v[459][2] = 1;
    v[460][2] = 7;
    v[461][2] = 3;
    v[462][2] = 7;
    v[463][2] = 3;
    v[464][2] = 5;
    v[465][2] = 7;
    v[466][2] = 3;
    v[467][2] = 5;
    v[468][2] = 3;
    v[469][2] = 1;
    v[470][2] = 1;
    v[471][2] = 1;
    v[472][2] = 5;
    v[473][2] = 7;
    v[474][2] = 7;
    v[475][2] = 3;
    v[476][2] = 3;
    v[477][2] = 1;
    v[478][2] = 1;
    v[479][2] = 1;
    v[480][2] = 5;
    v[481][2] = 5;
    v[482][2] = 7;
    v[483][2] = 3;
    v[484][2] = 1;
    v[485][2] = 1;
    v[486][2] = 3;
    v[487][2] = 3;
    v[488][2] = 7;
    v[489][2] = 3;
    v[490][2] = 3;
    v[491][2] = 5;
    v[492][2] = 1;
    v[493][2] = 3;
    v[494][2] = 7;
    v[495][2] = 3;
    v[496][2] = 3;
    v[497][2] = 7;
    v[498][2] = 3;
    v[499][2] = 5;
    v[500][2] = 7;
    v[501][2] = 5;
    v[502][2] = 7;
    v[503][2] = 7;
    v[504][2] = 3;
    v[505][2] = 3;
    v[506][2] = 5;
    v[507][2] = 1;
    v[508][2] = 3;
    v[509][2] = 5;
    v[510][2] = 3;
    v[511][2] = 1;
    v[512][2] = 3;
    v[513][2] = 5;
    v[514][2] = 1;
    v[515][2] = 1;
    v[516][2] = 3;
    v[517][2] = 7;
    v[518][2] = 7;
    v[519][2] = 1;
    v[520][2] = 5;
    v[521][2] = 1;
    v[522][2] = 3;
    v[523][2] = 7;
    v[524][2] = 3;
    v[525][2] = 7;
    v[526][2] = 3;
    v[527][2] = 5;
    v[528][2] = 1;
    v[529][2] = 7;
    v[530][2] = 1;
    v[531][2] = 1;
    v[532][2] = 3;
    v[533][2] = 5;
    v[534][2] = 3;
    v[535][2] = 7;
    v[536][2] = 1;
    v[537][2] = 5;
    v[538][2] = 5;
    v[539][2] = 1;
    v[540][2] = 1;
    v[541][2] = 3;
    v[542][2] = 1;
    v[543][2] = 3;
    v[544][2] = 3;
    v[545][2] = 7;
    v[546][2] = 1;
    v[547][2] = 7;
    v[548][2] = 3;
    v[549][2] = 1;
    v[550][2] = 7;
    v[551][2] = 3;
    v[552][2] = 1;
    v[553][2] = 7;
    v[554][2] = 3;
    v[555][2] = 5;
    v[556][2] = 3;
    v[557][2] = 5;
    v[558][2] = 7;
    v[559][2] = 3;
    v[560][2] = 3;
    v[561][2] = 3;
    v[562][2] = 5;
    v[563][2] = 1;
    v[564][2] = 7;
    v[565][2] = 7;
    v[566][2] = 1;
    v[567][2] = 3;
    v[568][2] = 1;
    v[569][2] = 3;
    v[570][2] = 7;
    v[571][2] = 7;
    v[572][2] = 1;
    v[573][2] = 3;
    v[574][2] = 7;
    v[575][2] = 3;
    v[576][2] = 1;
    v[577][2] = 5;
    v[578][2] = 3;
    v[579][2] = 1;
    v[580][2] = 1;
    v[581][2] = 1;
    v[582][2] = 5;
    v[583][2] = 3;
    v[584][2] = 3;
    v[585][2] = 7;
    v[586][2] = 1;
    v[587][2] = 5;
    v[588][2] = 3;
    v[589][2] = 5;
    v[590][2] = 1;
    v[591][2] = 3;
    v[592][2] = 1;
    v[593][2] = 3;
    v[594][2] = 1;
    v[595][2] = 5;
    v[596][2] = 7;
    v[597][2] = 7;
    v[598][2] = 1;
    v[599][2] = 1;
    v[600][2] = 5;
    v[601][2] = 3;
    v[602][2] = 1;
    v[603][2] = 5;
    v[604][2] = 1;
    v[605][2] = 1;
    v[606][2] = 7;
    v[607][2] = 7;
    v[608][2] = 3;
    v[609][2] = 5;
    v[610][2] = 5;
    v[611][2] = 1;
    v[612][2] = 7;
    v[613][2] = 1;
    v[614][2] = 5;
    v[615][2] = 1;
    v[616][2] = 1;
    v[617][2] = 3;
    v[618][2] = 1;
    v[619][2] = 5;
    v[620][2] = 7;
    v[621][2] = 5;
    v[622][2] = 7;
    v[623][2] = 7;
    v[624][2] = 1;
    v[625][2] = 5;
    v[626][2] = 1;
    v[627][2] = 1;
    v[628][2] = 3;
    v[629][2] = 5;
    v[630][2] = 1;
    v[631][2] = 5;
    v[632][2] = 5;
    v[633][2] = 3;
    v[634][2] = 1;
    v[635][2] = 3;
    v[636][2] = 1;
    v[637][2] = 5;
    v[638][2] = 5;
    v[639][2] = 3;
    v[640][2] = 3;
    v[641][2] = 3;
    v[642][2] = 3;
    v[643][2] = 1;
    v[644][2] = 1;
    v[645][2] = 3;
    v[646][2] = 1;
    v[647][2] = 3;
    v[648][2] = 5;
    v[649][2] = 5;
    v[650][2] = 7;
    v[651][2] = 5;
    v[652][2] = 5;
    v[653][2] = 7;
    v[654][2] = 5;
    v[655][2] = 7;
    v[656][2] = 1;
    v[657][2] = 3;
    v[658][2] = 7;
    v[659][2] = 7;
    v[660][2] = 3;
    v[661][2] = 5;
    v[662][2] = 5;
    v[663][2] = 7;
    v[664][2] = 5;
    v[665][2] = 5;
    v[666][2] = 3;
    v[667][2] = 3;
    v[668][2] = 3;
    v[669][2] = 1;
    v[670][2] = 7;
    v[671][2] = 1;
    v[672][2] = 5;
    v[673][2] = 5;
    v[674][2] = 5;
    v[675][2] = 3;
    v[676][2] = 3;
    v[677][2] = 5;
    v[678][2] = 1;
    v[679][2] = 3;
    v[680][2] = 1;
    v[681][2] = 3;
    v[682][2] = 3;
    v[683][2] = 3;
    v[684][2] = 7;
    v[685][2] = 1;
    v[686][2] = 7;
    v[687][2] = 7;
    v[688][2] = 3;
    v[689][2] = 7;
    v[690][2] = 1;
    v[691][2] = 1;
    v[692][2] = 5;
    v[693][2] = 7;
    v[694][2] = 1;
    v[695][2] = 7;
    v[696][2] = 1;
    v[697][2] = 7;
    v[698][2] = 7;
    v[699][2] = 1;
    v[700][2] = 3;
    v[701][2] = 7;
    v[702][2] = 5;
    v[703][2] = 1;
    v[704][2] = 3;
    v[705][2] = 5;
    v[706][2] = 5;
    v[707][2] = 5;
    v[708][2] = 1;
    v[709][2] = 1;
    v[710][2] = 7;
    v[711][2] = 1;
    v[712][2] = 7;
    v[713][2] = 1;
    v[714][2] = 7;
    v[715][2] = 7;
    v[716][2] = 3;
    v[717][2] = 1;
    v[718][2] = 1;
    v[719][2] = 5;
    v[720][2] = 1;
    v[721][2] = 5;
    v[722][2] = 1;
    v[723][2] = 5;
    v[724][2] = 3;
    v[725][2] = 5;
    v[726][2] = 5;
    v[727][2] = 5;
    v[728][2] = 5;
    v[729][2] = 5;
    v[730][2] = 3;
    v[731][2] = 3;
    v[732][2] = 7;
    v[733][2] = 3;
    v[734][2] = 3;
    v[735][2] = 5;
    v[736][2] = 5;
    v[737][2] = 3;
    v[738][2] = 7;
    v[739][2] = 1;
    v[740][2] = 5;
    v[741][2] = 7;
    v[742][2] = 5;
    v[743][2] = 1;
    v[744][2] = 5;
    v[745][2] = 5;
    v[746][2] = 3;
    v[747][2] = 5;
    v[748][2] = 5;
    v[749][2] = 7;
    v[750][2] = 5;
    v[751][2] = 3;
    v[752][2] = 5;
    v[753][2] = 5;
    v[754][2] = 5;
    v[755][2] = 1;
    v[756][2] = 5;
    v[757][2] = 5;
    v[758][2] = 5;
    v[759][2] = 5;
    v[760][2] = 1;
    v[761][2] = 3;
    v[762][2] = 5;
    v[763][2] = 3;
    v[764][2] = 1;
    v[765][2] = 7;
    v[766][2] = 5;
    v[767][2] = 5;
    v[768][2] = 7;
    v[769][2] = 1;
    v[770][2] = 5;
    v[771][2] = 3;
    v[772][2] = 3;
    v[773][2] = 1;
    v[774][2] = 5;
    v[775][2] = 3;
    v[776][2] = 7;
    v[777][2] = 1;
    v[778][2] = 7;
    v[779][2] = 5;
    v[780][2] = 1;
    v[781][2] = 1;
    v[782][2] = 3;
    v[783][2] = 1;
    v[784][2] = 1;
    v[785][2] = 7;
    v[786][2] = 1;
    v[787][2] = 5;
    v[788][2] = 5;
    v[789][2] = 3;
    v[790][2] = 7;
    v[791][2] = 3;
    v[792][2] = 7;
    v[793][2] = 5;
    v[794][2] = 3;
    v[795][2] = 1;
    v[796][2] = 1;
    v[797][2] = 3;
    v[798][2] = 1;
    v[799][2] = 3;
    v[800][2] = 5;
    v[801][2] = 5;
    v[802][2] = 7;
    v[803][2] = 5;
    v[804][2] = 3;
    v[805][2] = 7;
    v[806][2] = 7;
    v[807][2] = 7;
    v[808][2] = 3;
    v[809][2] = 7;
    v[810][2] = 3;
    v[811][2] = 7;
    v[812][2] = 1;
    v[813][2] = 3;
    v[814][2] = 1;
    v[815][2] = 7;
    v[816][2] = 7;
    v[817][2] = 1;
    v[818][2] = 7;
    v[819][2] = 3;
    v[820][2] = 7;
    v[821][2] = 3;
    v[822][2] = 7;
    v[823][2] = 3;
    v[824][2] = 7;
    v[825][2] = 3;
    v[826][2] = 5;
    v[827][2] = 1;
    v[828][2] = 1;
    v[829][2] = 7;
    v[830][2] = 3;
    v[831][2] = 1;
    v[832][2] = 5;
    v[833][2] = 5;
    v[834][2] = 7;
    v[835][2] = 1;
    v[836][2] = 5;
    v[837][2] = 5;
    v[838][2] = 5;
    v[839][2] = 7;
    v[840][2] = 1;
    v[841][2] = 5;
    v[842][2] = 5;
    v[843][2] = 1;
    v[844][2] = 5;
    v[845][2] = 5;
    v[846][2] = 3;
    v[847][2] = 1;
    v[848][2] = 3;
    v[849][2] = 1;
    v[850][2] = 7;
    v[851][2] = 3;
    v[852][2] = 1;
    v[853][2] = 3;
    v[854][2] = 5;
    v[855][2] = 7;
    v[856][2] = 7;
    v[857][2] = 7;
    v[858][2] = 1;
    v[859][2] = 1;
    v[860][2] = 7;
    v[861][2] = 3;
    v[862][2] = 1;
    v[863][2] = 5;
    v[864][2] = 5;
    v[865][2] = 5;
    v[866][2] = 1;
    v[867][2] = 1;
    v[868][2] = 1;
    v[869][2] = 1;
    v[870][2] = 1;
    v[871][2] = 5;
    v[872][2] = 3;
    v[873][2] = 5;
    v[874][2] = 1;
    v[875][2] = 3;
    v[876][2] = 5;
    v[877][2] = 3;
    v[878][2] = 1;
    v[879][2] = 1;
    v[880][2] = 1;
    v[881][2] = 1;
    v[882][2] = 3;
    v[883][2] = 7;
    v[884][2] = 3;
    v[885][2] = 7;
    v[886][2] = 5;
    v[887][2] = 7;
    v[888][2] = 1;
    v[889][2] = 5;
    v[890][2] = 5;
    v[891][2] = 7;
    v[892][2] = 5;
    v[893][2] = 3;
    v[894][2] = 3;
    v[895][2] = 7;
    v[896][2] = 5;
    v[897][2] = 3;
    v[898][2] = 1;
    v[899][2] = 1;
    v[900][2] = 3;
    v[901][2] = 1;
    v[902][2] = 3;
    v[903][2] = 1;
    v[904][2] = 1;
    v[905][2] = 3;
    v[906][2] = 7;
    v[907][2] = 1;
    v[908][2] = 7;
    v[909][2] = 1;
    v[910][2] = 1;
    v[911][2] = 5;
    v[912][2] = 1;
    v[913][2] = 7;
    v[914][2] = 5;
    v[915][2] = 3;
    v[916][2] = 7;
    v[917][2] = 3;
    v[918][2] = 5;
    v[919][2] = 3;
    v[920][2] = 1;
    v[921][2] = 1;
    v[922][2] = 5;
    v[923][2] = 5;
    v[924][2] = 1;
    v[925][2] = 7;
    v[926][2] = 7;
    v[927][2] = 3;
    v[928][2] = 7;
    v[929][2] = 3;
    v[930][2] = 7;
    v[931][2] = 1;
    v[932][2] = 5;
    v[933][2] = 1;
    v[934][2] = 5;
    v[935][2] = 3;
    v[936][2] = 7;
    v[937][2] = 3;
    v[938][2] = 5;
    v[939][2] = 7;
    v[940][2] = 7;
    v[941][2] = 7;
    v[942][2] = 3;
    v[943][2] = 3;
    v[944][2] = 1;
    v[945][2] = 1;
    v[946][2] = 5;
    v[947][2] = 5;
    v[948][2] = 3;
    v[949][2] = 7;
    v[950][2] = 1;
    v[951][2] = 1;
    v[952][2] = 1;
    v[953][2] = 3;
    v[954][2] = 5;
    v[955][2] = 3;
    v[956][2] = 1;
    v[957][2] = 1;
    v[958][2] = 3;
    v[959][2] = 3;
    v[960][2] = 7;
    v[961][2] = 5;
    v[962][2] = 1;
    v[963][2] = 1;
    v[964][2] = 3;
    v[965][2] = 7;
    v[966][2] = 1;
    v[967][2] = 5;
    v[968][2] = 7;
    v[969][2] = 3;
    v[970][2] = 7;
    v[971][2] = 5;
    v[972][2] = 5;
    v[973][2] = 7;
    v[974][2] = 3;
    v[975][2] = 5;
    v[976][2] = 3;
    v[977][2] = 1;
    v[978][2] = 5;
    v[979][2] = 3;
    v[980][2] = 1;
    v[981][2] = 1;
    v[982][2] = 7;
    v[983][2] = 5;
    v[984][2] = 1;
    v[985][2] = 7;
    v[986][2] = 3;
    v[987][2] = 7;
    v[988][2] = 5;
    v[989][2] = 1;
    v[990][2] = 7;
    v[991][2] = 1;
    v[992][2] = 7;
    v[993][2] = 7;
    v[994][2] = 1;
    v[995][2] = 1;
    v[996][2] = 7;
    v[997][2] = 1;
    v[998][2] = 5;
    v[999][2] = 5;
    v[1000][2] = 1;
    v[1001][2] = 1;
    v[1002][2] = 7;
    v[1003][2] = 5;
    v[1004][2] = 7;
    v[1005][2] = 1;
    v[1006][2] = 5;
    v[1007][2] = 3;
    v[1008][2] = 5;
    v[1009][2] = 3;
    v[1010][2] = 3;
    v[1011][2] = 7;
    v[1012][2] = 1;
    v[1013][2] = 5;
    v[1014][2] = 1;
    v[1015][2] = 1;
    v[1016][2] = 5;
    v[1017][2] = 5;
    v[1018][2] = 3;
    v[1019][2] = 3;
    v[1020][2] = 7;
    v[1021][2] = 5;
    v[1022][2] = 5;
    v[1023][2] = 1;
    v[1024][2] = 1;
    v[1025][2] = 1;
    v[1026][2] = 3;
    v[1027][2] = 1;
    v[1028][2] = 5;
    v[1029][2] = 7;
    v[1030][2] = 7;
    v[1031][2] = 1;
    v[1032][2] = 7;
    v[1033][2] = 5;
    v[1034][2] = 7;
    v[1035][2] = 3;
    v[1036][2] = 7;
    v[1037][2] = 3;
    v[1038][2] = 1;
    v[1039][2] = 3;
    v[1040][2] = 7;
    v[1041][2] = 3;
    v[1042][2] = 1;
    v[1043][2] = 5;
    v[1044][2] = 5;
    v[1045][2] = 3;
    v[1046][2] = 5;
    v[1047][2] = 1;
    v[1048][2] = 3;
    v[1049][2] = 5;
    v[1050][2] = 5;
    v[1051][2] = 5;
    v[1052][2] = 1;
    v[1053][2] = 1;
    v[1054][2] = 7;
    v[1055][2] = 7;
    v[1056][2] = 1;
    v[1057][2] = 5;
    v[1058][2] = 5;
    v[1059][2] = 1;
    v[1060][2] = 3;
    v[1061][2] = 5;
    v[1062][2] = 1;
    v[1063][2] = 5;
    v[1064][2] = 3;
    v[1065][2] = 5;
    v[1066][2] = 3;
    v[1067][2] = 3;
    v[1068][2] = 7;
    v[1069][2] = 5;
    v[1070][2] = 7;
    v[1071][2] = 3;
    v[1072][2] = 7;
    v[1073][2] = 3;
    v[1074][2] = 1;
    v[1075][2] = 3;
    v[1076][2] = 7;
    v[1077][2] = 7;
    v[1078][2] = 3;
    v[1079][2] = 3;
    v[1080][2] = 1;
    v[1081][2] = 1;
    v[1082][2] = 3;
    v[1083][2] = 3;
    v[1084][2] = 3;
    v[1085][2] = 3;
    v[1086][2] = 3;
    v[1087][2] = 5;
    v[1088][2] = 5;
    v[1089][2] = 3;
    v[1090][2] = 3;
    v[1091][2] = 3;
    v[1092][2] = 1;
    v[1093][2] = 3;
    v[1094][2] = 5;
    v[1095][2] = 7;
    v[1096][2] = 7;
    v[1097][2] = 1;
    v[1098][2] = 5;
    v[1099][2] = 7;
    v[1100][2] = 3;
    v[1101][2] = 7;
    v[1102][2] = 1;
    v[1103][2] = 1;
    v[1104][2] = 3;
    v[1105][2] = 5;
    v[1106][2] = 7;
    v[1107][2] = 5;
    v[1108][2] = 3;
    v[1109][2] = 3;
    v[1110][2] = 3;

    v[5][3] = 1;
    v[6][3] = 7;
    v[7][3] = 9;
    v[8][3] = 13;
    v[9][3] = 11;
    v[10][3] = 1;
    v[11][3] = 3;
    v[12][3] = 7;
    v[13][3] = 9;
    v[14][3] = 5;
    v[15][3] = 13;
    v[16][3] = 13;
    v[17][3] = 11;
    v[18][3] = 3;
    v[19][3] = 15;
    v[20][3] = 5;
    v[21][3] = 3;
    v[22][3] = 15;
    v[23][3] = 7;
    v[24][3] = 9;
    v[25][3] = 13;
    v[26][3] = 9;
    v[27][3] = 1;
    v[28][3] = 11;
    v[29][3] = 7;
    v[30][3] = 5;
    v[31][3] = 15;
    v[32][3] = 1;
    v[33][3] = 15;
    v[34][3] = 11;
    v[35][3] = 5;
    v[36][3] = 11;
    v[37][3] = 1;
    v[38][3] = 7;
    v[39][3] = 9;
    v[40][3] = 7;
    v[41][3] = 7;
    v[42][3] = 1;
    v[43][3] = 15;
    v[44][3] = 15;
    v[45][3] = 15;
    v[46][3] = 13;
    v[47][3] = 3;
    v[48][3] = 3;
    v[49][3] = 15;
    v[50][3] = 5;
    v[51][3] = 9;
    v[52][3] = 7;
    v[53][3] = 13;
    v[54][3] = 3;
    v[55][3] = 7;
    v[56][3] = 5;
    v[57][3] = 11;
    v[58][3] = 9;
    v[59][3] = 1;
    v[60][3] = 9;
    v[61][3] = 1;
    v[62][3] = 5;
    v[63][3] = 7;
    v[64][3] = 13;
    v[65][3] = 9;
    v[66][3] = 9;
    v[67][3] = 1;
    v[68][3] = 7;
    v[69][3] = 3;
    v[70][3] = 5;
    v[71][3] = 1;
    v[72][3] = 11;
    v[73][3] = 11;
    v[74][3] = 13;
    v[75][3] = 7;
    v[76][3] = 7;
    v[77][3] = 9;
    v[78][3] = 9;
    v[79][3] = 1;
    v[80][3] = 1;
    v[81][3] = 3;
    v[82][3] = 9;
    v[83][3] = 15;
    v[84][3] = 1;
    v[85][3] = 5;
    v[86][3] = 13;
    v[87][3] = 1;
    v[88][3] = 9;
    v[89][3] = 9;
    v[90][3] = 9;
    v[91][3] = 9;
    v[92][3] = 9;
    v[93][3] = 13;
    v[94][3] = 11;
    v[95][3] = 3;
    v[96][3] = 5;
    v[97][3] = 11;
    v[98][3] = 11;
    v[99][3] = 13;
    v[100][3] = 5;
    v[101][3] = 3;
    v[102][3] = 15;
    v[103][3] = 1;
    v[104][3] = 11;
    v[105][3] = 11;
    v[106][3] = 7;
    v[107][3] = 13;
    v[108][3] = 15;
    v[109][3] = 11;
    v[110][3] = 13;
    v[111][3] = 9;
    v[112][3] = 11;
    v[113][3] = 15;
    v[114][3] = 15;
    v[115][3] = 13;
    v[116][3] = 3;
    v[117][3] = 15;
    v[118][3] = 7;
    v[119][3] = 9;
    v[120][3] = 11;
    v[121][3] = 13;
    v[122][3] = 11;
    v[123][3] = 9;
    v[124][3] = 9;
    v[125][3] = 5;
    v[126][3] = 13;
    v[127][3] = 9;
    v[128][3] = 1;
    v[129][3] = 13;
    v[130][3] = 7;
    v[131][3] = 7;
    v[132][3] = 7;
    v[133][3] = 7;
    v[134][3] = 7;
    v[135][3] = 5;
    v[136][3] = 9;
    v[137][3] = 7;
    v[138][3] = 13;
    v[139][3] = 11;
    v[140][3] = 9;
    v[141][3] = 11;
    v[142][3] = 15;
    v[143][3] = 3;
    v[144][3] = 13;
    v[145][3] = 11;
    v[146][3] = 1;
    v[147][3] = 11;
    v[148][3] = 3;
    v[149][3] = 3;
    v[150][3] = 9;
    v[151][3] = 11;
    v[152][3] = 1;
    v[153][3] = 7;
    v[154][3] = 1;
    v[155][3] = 15;
    v[156][3] = 15;
    v[157][3] = 3;
    v[158][3] = 1;
    v[159][3] = 9;
    v[160][3] = 1;
    v[161][3] = 7;
    v[162][3] = 13;
    v[163][3] = 11;
    v[164][3] = 3;
    v[165][3] = 13;
    v[166][3] = 11;
    v[167][3] = 7;
    v[168][3] = 3;
    v[169][3] = 3;
    v[170][3] = 5;
    v[171][3] = 13;
    v[172][3] = 11;
    v[173][3] = 5;
    v[174][3] = 11;
    v[175][3] = 1;
    v[176][3] = 3;
    v[177][3] = 9;
    v[178][3] = 7;
    v[179][3] = 15;
    v[180][3] = 7;
    v[181][3] = 5;
    v[182][3] = 13;
    v[183][3] = 7;
    v[184][3] = 9;
    v[185][3] = 13;
    v[186][3] = 15;
    v[187][3] = 13;
    v[188][3] = 9;
    v[189][3] = 7;
    v[190][3] = 15;
    v[191][3] = 7;
    v[192][3] = 9;
    v[193][3] = 5;
    v[194][3] = 11;
    v[195][3] = 11;
    v[196][3] = 13;
    v[197][3] = 13;
    v[198][3] = 9;
    v[199][3] = 3;
    v[200][3] = 5;
    v[201][3] = 13;
    v[202][3] = 9;
    v[203][3] = 11;
    v[204][3] = 15;
    v[205][3] = 11;
    v[206][3] = 7;
    v[207][3] = 1;
    v[208][3] = 7;
    v[209][3] = 13;
    v[210][3] = 3;
    v[211][3] = 13;
    v[212][3] = 3;
    v[213][3] = 13;
    v[214][3] = 9;
    v[215][3] = 15;
    v[216][3] = 7;
    v[217][3] = 13;
    v[218][3] = 13;
    v[219][3] = 3;
    v[220][3] = 13;
    v[221][3] = 15;
    v[222][3] = 15;
    v[223][3] = 11;
    v[224][3] = 9;
    v[225][3] = 13;
    v[226][3] = 9;
    v[227][3] = 15;
    v[228][3] = 1;
    v[229][3] = 1;
    v[230][3] = 15;
    v[231][3] = 11;
    v[232][3] = 11;
    v[233][3] = 7;
    v[234][3] = 1;
    v[235][3] = 11;
    v[236][3] = 13;
    v[237][3] = 9;
    v[238][3] = 13;
    v[239][3] = 3;
    v[240][3] = 5;
    v[241][3] = 11;
    v[242][3] = 13;
    v[243][3] = 9;
    v[244][3] = 9;
    v[245][3] = 13;
    v[246][3] = 1;
    v[247][3] = 11;
    v[248][3] = 15;
    v[249][3] = 13;
    v[250][3] = 3;
    v[251][3] = 13;
    v[252][3] = 7;
    v[253][3] = 15;
    v[254][3] = 1;
    v[255][3] = 15;
    v[256][3] = 3;
    v[257][3] = 3;
    v[258][3] = 11;
    v[259][3] = 7;
    v[260][3] = 13;
    v[261][3] = 7;
    v[262][3] = 7;
    v[263][3] = 9;
    v[264][3] = 7;
    v[265][3] = 5;
    v[266][3] = 15;
    v[267][3] = 9;
    v[268][3] = 5;
    v[269][3] = 5;
    v[270][3] = 7;
    v[271][3] = 15;
    v[272][3] = 13;
    v[273][3] = 15;
    v[274][3] = 5;
    v[275][3] = 15;
    v[276][3] = 5;
    v[277][3] = 3;
    v[278][3] = 1;
    v[279][3] = 11;
    v[280][3] = 7;
    v[281][3] = 1;
    v[282][3] = 5;
    v[283][3] = 7;
    v[284][3] = 9;
    v[285][3] = 3;
    v[286][3] = 11;
    v[287][3] = 1;
    v[288][3] = 15;
    v[289][3] = 1;
    v[290][3] = 3;
    v[291][3] = 15;
    v[292][3] = 11;
    v[293][3] = 13;
    v[294][3] = 5;
    v[295][3] = 13;
    v[296][3] = 1;
    v[297][3] = 7;
    v[298][3] = 1;
    v[299][3] = 15;
    v[300][3] = 7;
    v[301][3] = 5;
    v[302][3] = 1;
    v[303][3] = 1;
    v[304][3] = 15;
    v[305][3] = 13;
    v[306][3] = 11;
    v[307][3] = 11;
    v[308][3] = 13;
    v[309][3] = 5;
    v[310][3] = 11;
    v[311][3] = 7;
    v[312][3] = 9;
    v[313][3] = 7;
    v[314][3] = 1;
    v[315][3] = 5;
    v[316][3] = 3;
    v[317][3] = 9;
    v[318][3] = 5;
    v[319][3] = 5;
    v[320][3] = 11;
    v[321][3] = 5;
    v[322][3] = 1;
    v[323][3] = 7;
    v[324][3] = 1;
    v[325][3] = 11;
    v[326][3] = 7;
    v[327][3] = 9;
    v[328][3] = 13;
    v[329][3] = 15;
    v[330][3] = 13;
    v[331][3] = 3;
    v[332][3] = 1;
    v[333][3] = 11;
    v[334][3] = 13;
    v[335][3] = 15;
    v[336][3] = 1;
    v[337][3] = 1;
    v[338][3] = 11;
    v[339][3] = 9;
    v[340][3] = 13;
    v[341][3] = 3;
    v[342][3] = 13;
    v[343][3] = 11;
    v[344][3] = 15;
    v[345][3] = 13;
    v[346][3] = 9;
    v[347][3] = 9;
    v[348][3] = 9;
    v[349][3] = 5;
    v[350][3] = 5;
    v[351][3] = 5;
    v[352][3] = 5;
    v[353][3] = 1;
    v[354][3] = 15;
    v[355][3] = 5;
    v[356][3] = 9;
    v[357][3] = 11;
    v[358][3] = 7;
    v[359][3] = 15;
    v[360][3] = 5;
    v[361][3] = 3;
    v[362][3] = 13;
    v[363][3] = 5;
    v[364][3] = 3;
    v[365][3] = 11;
    v[366][3] = 5;
    v[367][3] = 1;
    v[368][3] = 11;
    v[369][3] = 13;
    v[370][3] = 9;
    v[371][3] = 11;
    v[372][3] = 3;
    v[373][3] = 7;
    v[374][3] = 13;
    v[375][3] = 15;
    v[376][3] = 1;
    v[377][3] = 7;
    v[378][3] = 11;
    v[379][3] = 1;
    v[380][3] = 13;
    v[381][3] = 1;
    v[382][3] = 15;
    v[383][3] = 1;
    v[384][3] = 9;
    v[385][3] = 7;
    v[386][3] = 3;
    v[387][3] = 9;
    v[388][3] = 11;
    v[389][3] = 1;
    v[390][3] = 9;
    v[391][3] = 13;
    v[392][3] = 13;
    v[393][3] = 3;
    v[394][3] = 11;
    v[395][3] = 7;
    v[396][3] = 9;
    v[397][3] = 1;
    v[398][3] = 7;
    v[399][3] = 15;
    v[400][3] = 9;
    v[401][3] = 1;
    v[402][3] = 5;
    v[403][3] = 13;
    v[404][3] = 5;
    v[405][3] = 11;
    v[406][3] = 3;
    v[407][3] = 9;
    v[408][3] = 15;
    v[409][3] = 11;
    v[410][3] = 13;
    v[411][3] = 5;
    v[412][3] = 1;
    v[413][3] = 7;
    v[414][3] = 7;
    v[415][3] = 5;
    v[416][3] = 13;
    v[417][3] = 7;
    v[418][3] = 7;
    v[419][3] = 9;
    v[420][3] = 5;
    v[421][3] = 11;
    v[422][3] = 11;
    v[423][3] = 1;
    v[424][3] = 1;
    v[425][3] = 15;
    v[426][3] = 3;
    v[427][3] = 13;
    v[428][3] = 9;
    v[429][3] = 13;
    v[430][3] = 9;
    v[431][3] = 9;
    v[432][3] = 11;
    v[433][3] = 5;
    v[434][3] = 5;
    v[435][3] = 13;
    v[436][3] = 15;
    v[437][3] = 3;
    v[438][3] = 9;
    v[439][3] = 15;
    v[440][3] = 3;
    v[441][3] = 11;
    v[442][3] = 11;
    v[443][3] = 15;
    v[444][3] = 15;
    v[445][3] = 3;
    v[446][3] = 11;
    v[447][3] = 15;
    v[448][3] = 15;
    v[449][3] = 3;
    v[450][3] = 1;
    v[451][3] = 3;
    v[452][3] = 1;
    v[453][3] = 3;
    v[454][3] = 3;
    v[455][3] = 1;
    v[456][3] = 3;
    v[457][3] = 13;
    v[458][3] = 1;
    v[459][3] = 11;
    v[460][3] = 5;
    v[461][3] = 15;
    v[462][3] = 7;
    v[463][3] = 15;
    v[464][3] = 9;
    v[465][3] = 1;
    v[466][3] = 7;
    v[467][3] = 1;
    v[468][3] = 9;
    v[469][3] = 11;
    v[470][3] = 15;
    v[471][3] = 1;
    v[472][3] = 13;
    v[473][3] = 9;
    v[474][3] = 13;
    v[475][3] = 11;
    v[476][3] = 7;
    v[477][3] = 3;
    v[478][3] = 7;
    v[479][3] = 3;
    v[480][3] = 13;
    v[481][3] = 7;
    v[482][3] = 9;
    v[483][3] = 7;
    v[484][3] = 7;
    v[485][3] = 3;
    v[486][3] = 3;
    v[487][3] = 9;
    v[488][3] = 9;
    v[489][3] = 7;
    v[490][3] = 5;
    v[491][3] = 11;
    v[492][3] = 13;
    v[493][3] = 13;
    v[494][3] = 7;
    v[495][3] = 7;
    v[496][3] = 15;
    v[497][3] = 9;
    v[498][3] = 5;
    v[499][3] = 5;
    v[500][3] = 3;
    v[501][3] = 3;
    v[502][3] = 13;
    v[503][3] = 3;
    v[504][3] = 9;
    v[505][3] = 3;
    v[506][3] = 1;
    v[507][3] = 11;
    v[508][3] = 1;
    v[509][3] = 3;
    v[510][3] = 11;
    v[511][3] = 15;
    v[512][3] = 11;
    v[513][3] = 11;
    v[514][3] = 11;
    v[515][3] = 9;
    v[516][3] = 13;
    v[517][3] = 7;
    v[518][3] = 9;
    v[519][3] = 15;
    v[520][3] = 9;
    v[521][3] = 11;
    v[522][3] = 1;
    v[523][3] = 3;
    v[524][3] = 3;
    v[525][3] = 9;
    v[526][3] = 7;
    v[527][3] = 15;
    v[528][3] = 13;
    v[529][3] = 13;
    v[530][3] = 7;
    v[531][3] = 15;
    v[532][3] = 9;
    v[533][3] = 13;
    v[534][3] = 9;
    v[535][3] = 15;
    v[536][3] = 13;
    v[537][3] = 15;
    v[538][3] = 9;
    v[539][3] = 13;
    v[540][3] = 1;
    v[541][3] = 11;
    v[542][3] = 7;
    v[543][3] = 11;
    v[544][3] = 3;
    v[545][3] = 13;
    v[546][3] = 5;
    v[547][3] = 1;
    v[548][3] = 7;
    v[549][3] = 15;
    v[550][3] = 3;
    v[551][3] = 13;
    v[552][3] = 7;
    v[553][3] = 13;
    v[554][3] = 13;
    v[555][3] = 11;
    v[556][3] = 3;
    v[557][3] = 5;
    v[558][3] = 3;
    v[559][3] = 13;
    v[560][3] = 11;
    v[561][3] = 9;
    v[562][3] = 9;
    v[563][3] = 3;
    v[564][3] = 11;
    v[565][3] = 11;
    v[566][3] = 7;
    v[567][3] = 9;
    v[568][3] = 13;
    v[569][3] = 11;
    v[570][3] = 7;
    v[571][3] = 15;
    v[572][3] = 13;
    v[573][3] = 7;
    v[574][3] = 5;
    v[575][3] = 3;
    v[576][3] = 1;
    v[577][3] = 5;
    v[578][3] = 15;
    v[579][3] = 15;
    v[580][3] = 3;
    v[581][3] = 11;
    v[582][3] = 1;
    v[583][3] = 7;
    v[584][3] = 3;
    v[585][3] = 15;
    v[586][3] = 11;
    v[587][3] = 5;
    v[588][3] = 5;
    v[589][3] = 3;
    v[590][3] = 5;
    v[591][3] = 5;
    v[592][3] = 1;
    v[593][3] = 15;
    v[594][3] = 5;
    v[595][3] = 1;
    v[596][3] = 5;
    v[597][3] = 3;
    v[598][3] = 7;
    v[599][3] = 5;
    v[600][3] = 11;
    v[601][3] = 3;
    v[602][3] = 13;
    v[603][3] = 9;
    v[604][3] = 13;
    v[605][3] = 15;
    v[606][3] = 5;
    v[607][3] = 3;
    v[608][3] = 5;
    v[609][3] = 9;
    v[610][3] = 5;
    v[611][3] = 3;
    v[612][3] = 11;
    v[613][3] = 1;
    v[614][3] = 13;
    v[615][3] = 9;
    v[616][3] = 15;
    v[617][3] = 3;
    v[618][3] = 5;
    v[619][3] = 11;
    v[620][3] = 9;
    v[621][3] = 1;
    v[622][3] = 3;
    v[623][3] = 15;
    v[624][3] = 9;
    v[625][3] = 9;
    v[626][3] = 9;
    v[627][3] = 11;
    v[628][3] = 7;
    v[629][3] = 5;
    v[630][3] = 13;
    v[631][3] = 1;
    v[632][3] = 15;
    v[633][3] = 3;
    v[634][3] = 13;
    v[635][3] = 9;
    v[636][3] = 13;
    v[637][3] = 5;
    v[638][3] = 1;
    v[639][3] = 5;
    v[640][3] = 1;
    v[641][3] = 13;
    v[642][3] = 13;
    v[643][3] = 7;
    v[644][3] = 7;
    v[645][3] = 1;
    v[646][3] = 9;
    v[647][3] = 5;
    v[648][3] = 11;
    v[649][3] = 9;
    v[650][3] = 11;
    v[651][3] = 13;
    v[652][3] = 3;
    v[653][3] = 15;
    v[654][3] = 15;
    v[655][3] = 13;
    v[656][3] = 15;
    v[657][3] = 7;
    v[658][3] = 5;
    v[659][3] = 7;
    v[660][3] = 9;
    v[661][3] = 7;
    v[662][3] = 9;
    v[663][3] = 9;
    v[664][3] = 9;
    v[665][3] = 11;
    v[666][3] = 9;
    v[667][3] = 3;
    v[668][3] = 11;
    v[669][3] = 15;
    v[670][3] = 13;
    v[671][3] = 13;
    v[672][3] = 5;
    v[673][3] = 9;
    v[674][3] = 15;
    v[675][3] = 1;
    v[676][3] = 1;
    v[677][3] = 9;
    v[678][3] = 5;
    v[679][3] = 13;
    v[680][3] = 3;
    v[681][3] = 13;
    v[682][3] = 15;
    v[683][3] = 3;
    v[684][3] = 1;
    v[685][3] = 3;
    v[686][3] = 11;
    v[687][3] = 13;
    v[688][3] = 1;
    v[689][3] = 15;
    v[690][3] = 9;
    v[691][3] = 9;
    v[692][3] = 3;
    v[693][3] = 1;
    v[694][3] = 9;
    v[695][3] = 1;
    v[696][3] = 9;
    v[697][3] = 1;
    v[698][3] = 13;
    v[699][3] = 11;
    v[700][3] = 15;
    v[701][3] = 7;
    v[702][3] = 11;
    v[703][3] = 15;
    v[704][3] = 13;
    v[705][3] = 15;
    v[706][3] = 1;
    v[707][3] = 9;
    v[708][3] = 9;
    v[709][3] = 7;
    v[710][3] = 3;
    v[711][3] = 5;
    v[712][3] = 11;
    v[713][3] = 7;
    v[714][3] = 3;
    v[715][3] = 9;
    v[716][3] = 5;
    v[717][3] = 15;
    v[718][3] = 7;
    v[719][3] = 5;
    v[720][3] = 3;
    v[721][3] = 13;
    v[722][3] = 7;
    v[723][3] = 1;
    v[724][3] = 1;
    v[725][3] = 9;
    v[726][3] = 15;
    v[727][3] = 15;
    v[728][3] = 15;
    v[729][3] = 11;
    v[730][3] = 3;
    v[731][3] = 5;
    v[732][3] = 15;
    v[733][3] = 13;
    v[734][3] = 7;
    v[735][3] = 15;
    v[736][3] = 15;
    v[737][3] = 11;
    v[738][3] = 11;
    v[739][3] = 9;
    v[740][3] = 5;
    v[741][3] = 15;
    v[742][3] = 9;
    v[743][3] = 7;
    v[744][3] = 3;
    v[745][3] = 13;
    v[746][3] = 1;
    v[747][3] = 1;
    v[748][3] = 5;
    v[749][3] = 1;
    v[750][3] = 3;
    v[751][3] = 1;
    v[752][3] = 7;
    v[753][3] = 1;
    v[754][3] = 1;
    v[755][3] = 5;
    v[756][3] = 1;
    v[757][3] = 11;
    v[758][3] = 11;
    v[759][3] = 9;
    v[760][3] = 9;
    v[761][3] = 5;
    v[762][3] = 13;
    v[763][3] = 7;
    v[764][3] = 7;
    v[765][3] = 7;
    v[766][3] = 1;
    v[767][3] = 1;
    v[768][3] = 9;
    v[769][3] = 9;
    v[770][3] = 11;
    v[771][3] = 11;
    v[772][3] = 15;
    v[773][3] = 7;
    v[774][3] = 5;
    v[775][3] = 5;
    v[776][3] = 3;
    v[777][3] = 11;
    v[778][3] = 1;
    v[779][3] = 3;
    v[780][3] = 7;
    v[781][3] = 13;
    v[782][3] = 7;
    v[783][3] = 7;
    v[784][3] = 7;
    v[785][3] = 3;
    v[786][3] = 15;
    v[787][3] = 15;
    v[788][3] = 11;
    v[789][3] = 9;
    v[790][3] = 3;
    v[791][3] = 9;
    v[792][3] = 3;
    v[793][3] = 15;
    v[794][3] = 13;
    v[795][3] = 5;
    v[796][3] = 3;
    v[797][3] = 3;
    v[798][3] = 3;
    v[799][3] = 5;
    v[800][3] = 9;
    v[801][3] = 15;
    v[802][3] = 9;
    v[803][3] = 9;
    v[804][3] = 1;
    v[805][3] = 5;
    v[806][3] = 9;
    v[807][3] = 9;
    v[808][3] = 15;
    v[809][3] = 5;
    v[810][3] = 15;
    v[811][3] = 7;
    v[812][3] = 9;
    v[813][3] = 1;
    v[814][3] = 9;
    v[815][3] = 9;
    v[816][3] = 5;
    v[817][3] = 11;
    v[818][3] = 5;
    v[819][3] = 15;
    v[820][3] = 15;
    v[821][3] = 11;
    v[822][3] = 7;
    v[823][3] = 7;
    v[824][3] = 7;
    v[825][3] = 1;
    v[826][3] = 1;
    v[827][3] = 11;
    v[828][3] = 11;
    v[829][3] = 13;
    v[830][3] = 15;
    v[831][3] = 3;
    v[832][3] = 13;
    v[833][3] = 5;
    v[834][3] = 1;
    v[835][3] = 7;
    v[836][3] = 1;
    v[837][3] = 11;
    v[838][3] = 3;
    v[839][3] = 13;
    v[840][3] = 15;
    v[841][3] = 3;
    v[842][3] = 5;
    v[843][3] = 3;
    v[844][3] = 5;
    v[845][3] = 7;
    v[846][3] = 3;
    v[847][3] = 9;
    v[848][3] = 9;
    v[849][3] = 5;
    v[850][3] = 1;
    v[851][3] = 7;
    v[852][3] = 11;
    v[853][3] = 9;
    v[854][3] = 3;
    v[855][3] = 5;
    v[856][3] = 11;
    v[857][3] = 13;
    v[858][3] = 13;
    v[859][3] = 13;
    v[860][3] = 9;
    v[861][3] = 15;
    v[862][3] = 5;
    v[863][3] = 7;
    v[864][3] = 1;
    v[865][3] = 15;
    v[866][3] = 11;
    v[867][3] = 9;
    v[868][3] = 15;
    v[869][3] = 15;
    v[870][3] = 13;
    v[871][3] = 13;
    v[872][3] = 13;
    v[873][3] = 1;
    v[874][3] = 11;
    v[875][3] = 9;
    v[876][3] = 15;
    v[877][3] = 9;
    v[878][3] = 5;
    v[879][3] = 15;
    v[880][3] = 5;
    v[881][3] = 7;
    v[882][3] = 3;
    v[883][3] = 11;
    v[884][3] = 3;
    v[885][3] = 15;
    v[886][3] = 7;
    v[887][3] = 13;
    v[888][3] = 11;
    v[889][3] = 7;
    v[890][3] = 3;
    v[891][3] = 7;
    v[892][3] = 13;
    v[893][3] = 5;
    v[894][3] = 13;
    v[895][3] = 15;
    v[896][3] = 5;
    v[897][3] = 13;
    v[898][3] = 9;
    v[899][3] = 1;
    v[900][3] = 15;
    v[901][3] = 11;
    v[902][3] = 5;
    v[903][3] = 5;
    v[904][3] = 1;
    v[905][3] = 11;
    v[906][3] = 3;
    v[907][3] = 3;
    v[908][3] = 7;
    v[909][3] = 1;
    v[910][3] = 9;
    v[911][3] = 7;
    v[912][3] = 15;
    v[913][3] = 9;
    v[914][3] = 9;
    v[915][3] = 3;
    v[916][3] = 11;
    v[917][3] = 15;
    v[918][3] = 7;
    v[919][3] = 1;
    v[920][3] = 3;
    v[921][3] = 1;
    v[922][3] = 1;
    v[923][3] = 1;
    v[924][3] = 9;
    v[925][3] = 1;
    v[926][3] = 5;
    v[927][3] = 15;
    v[928][3] = 15;
    v[929][3] = 7;
    v[930][3] = 5;
    v[931][3] = 5;
    v[932][3] = 7;
    v[933][3] = 9;
    v[934][3] = 7;
    v[935][3] = 15;
    v[936][3] = 13;
    v[937][3] = 13;
    v[938][3] = 11;
    v[939][3] = 1;
    v[940][3] = 9;
    v[941][3] = 11;
    v[942][3] = 1;
    v[943][3] = 13;
    v[944][3] = 1;
    v[945][3] = 7;
    v[946][3] = 15;
    v[947][3] = 15;
    v[948][3] = 5;
    v[949][3] = 5;
    v[950][3] = 1;
    v[951][3] = 11;
    v[952][3] = 3;
    v[953][3] = 9;
    v[954][3] = 11;
    v[955][3] = 9;
    v[956][3] = 9;
    v[957][3] = 9;
    v[958][3] = 1;
    v[959][3] = 9;
    v[960][3] = 3;
    v[961][3] = 5;
    v[962][3] = 15;
    v[963][3] = 1;
    v[964][3] = 1;
    v[965][3] = 9;
    v[966][3] = 7;
    v[967][3] = 3;
    v[968][3] = 3;
    v[969][3] = 1;
    v[970][3] = 9;
    v[971][3] = 9;
    v[972][3] = 11;
    v[973][3] = 9;
    v[974][3] = 9;
    v[975][3] = 13;
    v[976][3] = 13;
    v[977][3] = 3;
    v[978][3] = 13;
    v[979][3] = 11;
    v[980][3] = 13;
    v[981][3] = 5;
    v[982][3] = 1;
    v[983][3] = 5;
    v[984][3] = 5;
    v[985][3] = 9;
    v[986][3] = 9;
    v[987][3] = 3;
    v[988][3] = 13;
    v[989][3] = 13;
    v[990][3] = 9;
    v[991][3] = 15;
    v[992][3] = 9;
    v[993][3] = 11;
    v[994][3] = 7;
    v[995][3] = 11;
    v[996][3] = 9;
    v[997][3] = 13;
    v[998][3] = 9;
    v[999][3] = 1;
    v[1000][3] = 15;
    v[1001][3] = 9;
    v[1002][3] = 7;
    v[1003][3] = 7;
    v[1004][3] = 1;
    v[1005][3] = 7;
    v[1006][3] = 9;
    v[1007][3] = 9;
    v[1008][3] = 15;
    v[1009][3] = 1;
    v[1010][3] = 11;
    v[1011][3] = 1;
    v[1012][3] = 13;
    v[1013][3] = 13;
    v[1014][3] = 15;
    v[1015][3] = 9;
    v[1016][3] = 13;
    v[1017][3] = 7;
    v[1018][3] = 15;
    v[1019][3] = 3;
    v[1020][3] = 9;
    v[1021][3] = 3;
    v[1022][3] = 1;
    v[1023][3] = 13;
    v[1024][3] = 7;
    v[1025][3] = 5;
    v[1026][3] = 9;
    v[1027][3] = 3;
    v[1028][3] = 1;
    v[1029][3] = 7;
    v[1030][3] = 1;
    v[1031][3] = 1;
    v[1032][3] = 13;
    v[1033][3] = 3;
    v[1034][3] = 3;
    v[1035][3] = 11;
    v[1036][3] = 1;
    v[1037][3] = 7;
    v[1038][3] = 13;
    v[1039][3] = 15;
    v[1040][3] = 15;
    v[1041][3] = 5;
    v[1042][3] = 7;
    v[1043][3] = 13;
    v[1044][3] = 13;
    v[1045][3] = 15;
    v[1046][3] = 11;
    v[1047][3] = 13;
    v[1048][3] = 1;
    v[1049][3] = 13;
    v[1050][3] = 13;
    v[1051][3] = 3;
    v[1052][3] = 9;
    v[1053][3] = 15;
    v[1054][3] = 15;
    v[1055][3] = 11;
    v[1056][3] = 15;
    v[1057][3] = 9;
    v[1058][3] = 15;
    v[1059][3] = 1;
    v[1060][3] = 13;
    v[1061][3] = 15;
    v[1062][3] = 1;
    v[1063][3] = 1;
    v[1064][3] = 5;
    v[1065][3] = 11;
    v[1066][3] = 5;
    v[1067][3] = 1;
    v[1068][3] = 11;
    v[1069][3] = 11;
    v[1070][3] = 5;
    v[1071][3] = 3;
    v[1072][3] = 9;
    v[1073][3] = 1;
    v[1074][3] = 3;
    v[1075][3] = 5;
    v[1076][3] = 13;
    v[1077][3] = 9;
    v[1078][3] = 7;
    v[1079][3] = 7;
    v[1080][3] = 1;
    v[1081][3] = 9;
    v[1082][3] = 9;
    v[1083][3] = 15;
    v[1084][3] = 7;
    v[1085][3] = 5;
    v[1086][3] = 5;
    v[1087][3] = 15;
    v[1088][3] = 13;
    v[1089][3] = 9;
    v[1090][3] = 7;
    v[1091][3] = 13;
    v[1092][3] = 3;
    v[1093][3] = 13;
    v[1094][3] = 11;
    v[1095][3] = 13;
    v[1096][3] = 7;
    v[1097][3] = 9;
    v[1098][3] = 13;
    v[1099][3] = 13;
    v[1100][3] = 13;
    v[1101][3] = 15;
    v[1102][3] = 9;
    v[1103][3] = 5;
    v[1104][3] = 5;
    v[1105][3] = 3;
    v[1106][3] = 3;
    v[1107][3] = 3;
    v[1108][3] = 1;
    v[1109][3] = 3;
    v[1110][3] = 15;

    v[7][4] = 9;
    v[8][4] = 3;
    v[9][4] = 27;
    v[10][4] = 15;
    v[11][4] = 29;
    v[12][4] = 21;
    v[13][4] = 23;
    v[14][4] = 19;
    v[15][4] = 11;
    v[16][4] = 25;
    v[17][4] = 7;
    v[18][4] = 13;
    v[19][4] = 17;
    v[20][4] = 1;
    v[21][4] = 25;
    v[22][4] = 29;
    v[23][4] = 3;
    v[24][4] = 31;
    v[25][4] = 11;
    v[26][4] = 5;
    v[27][4] = 23;
    v[28][4] = 27;
    v[29][4] = 19;
    v[30][4] = 21;
    v[31][4] = 5;
    v[32][4] = 1;
    v[33][4] = 17;
    v[34][4] = 13;
    v[35][4] = 7;
    v[36][4] = 15;
    v[37][4] = 9;
    v[38][4] = 31;
    v[39][4] = 25;
    v[40][4] = 3;
    v[41][4] = 5;
    v[42][4] = 23;
    v[43][4] = 7;
    v[44][4] = 3;
    v[45][4] = 17;
    v[46][4] = 23;
    v[47][4] = 3;
    v[48][4] = 3;
    v[49][4] = 21;
    v[50][4] = 25;
    v[51][4] = 25;
    v[52][4] = 23;
    v[53][4] = 11;
    v[54][4] = 19;
    v[55][4] = 3;
    v[56][4] = 11;
    v[57][4] = 31;
    v[58][4] = 7;
    v[59][4] = 9;
    v[60][4] = 5;
    v[61][4] = 17;
    v[62][4] = 23;
    v[63][4] = 17;
    v[64][4] = 17;
    v[65][4] = 25;
    v[66][4] = 13;
    v[67][4] = 11;
    v[68][4] = 31;
    v[69][4] = 27;
    v[70][4] = 19;
    v[71][4] = 17;
    v[72][4] = 23;
    v[73][4] = 7;
    v[74][4] = 5;
    v[75][4] = 11;
    v[76][4] = 19;
    v[77][4] = 19;
    v[78][4] = 7;
    v[79][4] = 13;
    v[80][4] = 21;
    v[81][4] = 21;
    v[82][4] = 7;
    v[83][4] = 9;
    v[84][4] = 11;
    v[85][4] = 1;
    v[86][4] = 5;
    v[87][4] = 21;
    v[88][4] = 11;
    v[89][4] = 13;
    v[90][4] = 25;
    v[91][4] = 9;
    v[92][4] = 7;
    v[93][4] = 7;
    v[94][4] = 27;
    v[95][4] = 15;
    v[96][4] = 25;
    v[97][4] = 15;
    v[98][4] = 21;
    v[99][4] = 17;
    v[100][4] = 19;
    v[101][4] = 19;
    v[102][4] = 21;
    v[103][4] = 5;
    v[104][4] = 11;
    v[105][4] = 3;
    v[106][4] = 5;
    v[107][4] = 29;
    v[108][4] = 31;
    v[109][4] = 29;
    v[110][4] = 5;
    v[111][4] = 5;
    v[112][4] = 1;
    v[113][4] = 31;
    v[114][4] = 27;
    v[115][4] = 11;
    v[116][4] = 13;
    v[117][4] = 1;
    v[118][4] = 3;
    v[119][4] = 7;
    v[120][4] = 11;
    v[121][4] = 7;
    v[122][4] = 3;
    v[123][4] = 23;
    v[124][4] = 13;
    v[125][4] = 31;
    v[126][4] = 17;
    v[127][4] = 1;
    v[128][4] = 27;
    v[129][4] = 11;
    v[130][4] = 25;
    v[131][4] = 1;
    v[132][4] = 23;
    v[133][4] = 29;
    v[134][4] = 17;
    v[135][4] = 25;
    v[136][4] = 7;
    v[137][4] = 25;
    v[138][4] = 27;
    v[139][4] = 17;
    v[140][4] = 13;
    v[141][4] = 17;
    v[142][4] = 23;
    v[143][4] = 5;
    v[144][4] = 17;
    v[145][4] = 5;
    v[146][4] = 13;
    v[147][4] = 11;
    v[148][4] = 21;
    v[149][4] = 5;
    v[150][4] = 11;
    v[151][4] = 5;
    v[152][4] = 9;
    v[153][4] = 31;
    v[154][4] = 19;
    v[155][4] = 17;
    v[156][4] = 9;
    v[157][4] = 9;
    v[158][4] = 27;
    v[159][4] = 21;
    v[160][4] = 15;
    v[161][4] = 15;
    v[162][4] = 1;
    v[163][4] = 1;
    v[164][4] = 29;
    v[165][4] = 5;
    v[166][4] = 31;
    v[167][4] = 11;
    v[168][4] = 17;
    v[169][4] = 23;
    v[170][4] = 19;
    v[171][4] = 21;
    v[172][4] = 25;
    v[173][4] = 15;
    v[174][4] = 11;
    v[175][4] = 5;
    v[176][4] = 5;
    v[177][4] = 1;
    v[178][4] = 19;
    v[179][4] = 19;
    v[180][4] = 19;
    v[181][4] = 7;
    v[182][4] = 13;
    v[183][4] = 21;
    v[184][4] = 17;
    v[185][4] = 17;
    v[186][4] = 25;
    v[187][4] = 23;
    v[188][4] = 19;
    v[189][4] = 23;
    v[190][4] = 15;
    v[191][4] = 13;
    v[192][4] = 5;
    v[193][4] = 19;
    v[194][4] = 25;
    v[195][4] = 9;
    v[196][4] = 7;
    v[197][4] = 3;
    v[198][4] = 21;
    v[199][4] = 17;
    v[200][4] = 25;
    v[201][4] = 1;
    v[202][4] = 27;
    v[203][4] = 25;
    v[204][4] = 27;
    v[205][4] = 25;
    v[206][4] = 9;
    v[207][4] = 13;
    v[208][4] = 3;
    v[209][4] = 17;
    v[210][4] = 25;
    v[211][4] = 23;
    v[212][4] = 9;
    v[213][4] = 25;
    v[214][4] = 9;
    v[215][4] = 13;
    v[216][4] = 17;
    v[217][4] = 17;
    v[218][4] = 3;
    v[219][4] = 15;
    v[220][4] = 7;
    v[221][4] = 7;
    v[222][4] = 29;
    v[223][4] = 3;
    v[224][4] = 19;
    v[225][4] = 29;
    v[226][4] = 29;
    v[227][4] = 19;
    v[228][4] = 29;
    v[229][4] = 13;
    v[230][4] = 15;
    v[231][4] = 25;
    v[232][4] = 27;
    v[233][4] = 1;
    v[234][4] = 3;
    v[235][4] = 9;
    v[236][4] = 9;
    v[237][4] = 13;
    v[238][4] = 31;
    v[239][4] = 29;
    v[240][4] = 31;
    v[241][4] = 5;
    v[242][4] = 15;
    v[243][4] = 29;
    v[244][4] = 1;
    v[245][4] = 19;
    v[246][4] = 5;
    v[247][4] = 9;
    v[248][4] = 19;
    v[249][4] = 5;
    v[250][4] = 15;
    v[251][4] = 3;
    v[252][4] = 5;
    v[253][4] = 7;
    v[254][4] = 15;
    v[255][4] = 17;
    v[256][4] = 17;
    v[257][4] = 23;
    v[258][4] = 11;
    v[259][4] = 9;
    v[260][4] = 23;
    v[261][4] = 19;
    v[262][4] = 3;
    v[263][4] = 17;
    v[264][4] = 1;
    v[265][4] = 27;
    v[266][4] = 9;
    v[267][4] = 9;
    v[268][4] = 17;
    v[269][4] = 13;
    v[270][4] = 25;
    v[271][4] = 29;
    v[272][4] = 23;
    v[273][4] = 29;
    v[274][4] = 11;
    v[275][4] = 31;
    v[276][4] = 25;
    v[277][4] = 21;
    v[278][4] = 29;
    v[279][4] = 19;
    v[280][4] = 27;
    v[281][4] = 31;
    v[282][4] = 3;
    v[283][4] = 5;
    v[284][4] = 3;
    v[285][4] = 3;
    v[286][4] = 13;
    v[287][4] = 21;
    v[288][4] = 9;
    v[289][4] = 29;
    v[290][4] = 3;
    v[291][4] = 17;
    v[292][4] = 11;
    v[293][4] = 11;
    v[294][4] = 9;
    v[295][4] = 21;
    v[296][4] = 19;
    v[297][4] = 7;
    v[298][4] = 17;
    v[299][4] = 31;
    v[300][4] = 25;
    v[301][4] = 1;
    v[302][4] = 27;
    v[303][4] = 5;
    v[304][4] = 15;
    v[305][4] = 27;
    v[306][4] = 29;
    v[307][4] = 29;
    v[308][4] = 29;
    v[309][4] = 25;
    v[310][4] = 27;
    v[311][4] = 25;
    v[312][4] = 3;
    v[313][4] = 21;
    v[314][4] = 17;
    v[315][4] = 25;
    v[316][4] = 13;
    v[317][4] = 15;
    v[318][4] = 17;
    v[319][4] = 13;
    v[320][4] = 23;
    v[321][4] = 9;
    v[322][4] = 3;
    v[323][4] = 11;
    v[324][4] = 7;
    v[325][4] = 9;
    v[326][4] = 9;
    v[327][4] = 7;
    v[328][4] = 17;
    v[329][4] = 7;
    v[330][4] = 1;
    v[331][4] = 27;
    v[332][4] = 1;
    v[333][4] = 9;
    v[334][4] = 5;
    v[335][4] = 31;
    v[336][4] = 21;
    v[337][4] = 25;
    v[338][4] = 25;
    v[339][4] = 21;
    v[340][4] = 11;
    v[341][4] = 1;
    v[342][4] = 23;
    v[343][4] = 19;
    v[344][4] = 27;
    v[345][4] = 15;
    v[346][4] = 3;
    v[347][4] = 5;
    v[348][4] = 23;
    v[349][4] = 9;
    v[350][4] = 25;
    v[351][4] = 7;
    v[352][4] = 29;
    v[353][4] = 11;
    v[354][4] = 9;
    v[355][4] = 13;
    v[356][4] = 5;
    v[357][4] = 11;
    v[358][4] = 1;
    v[359][4] = 3;
    v[360][4] = 31;
    v[361][4] = 27;
    v[362][4] = 3;
    v[363][4] = 17;
    v[364][4] = 27;
    v[365][4] = 11;
    v[366][4] = 13;
    v[367][4] = 15;
    v[368][4] = 29;
    v[369][4] = 15;
    v[370][4] = 1;
    v[371][4] = 15;
    v[372][4] = 23;
    v[373][4] = 25;
    v[374][4] = 13;
    v[375][4] = 21;
    v[376][4] = 15;
    v[377][4] = 3;
    v[378][4] = 29;
    v[379][4] = 29;
    v[380][4] = 5;
    v[381][4] = 25;
    v[382][4] = 17;
    v[383][4] = 11;
    v[384][4] = 7;
    v[385][4] = 15;
    v[386][4] = 5;
    v[387][4] = 21;
    v[388][4] = 7;
    v[389][4] = 31;
    v[390][4] = 13;
    v[391][4] = 11;
    v[392][4] = 23;
    v[393][4] = 5;
    v[394][4] = 7;
    v[395][4] = 23;
    v[396][4] = 27;
    v[397][4] = 21;
    v[398][4] = 29;
    v[399][4] = 15;
    v[400][4] = 7;
    v[401][4] = 27;
    v[402][4] = 27;
    v[403][4] = 19;
    v[404][4] = 7;
    v[405][4] = 15;
    v[406][4] = 27;
    v[407][4] = 27;
    v[408][4] = 19;
    v[409][4] = 19;
    v[410][4] = 9;
    v[411][4] = 15;
    v[412][4] = 1;
    v[413][4] = 3;
    v[414][4] = 29;
    v[415][4] = 29;
    v[416][4] = 5;
    v[417][4] = 27;
    v[418][4] = 31;
    v[419][4] = 9;
    v[420][4] = 1;
    v[421][4] = 7;
    v[422][4] = 3;
    v[423][4] = 19;
    v[424][4] = 19;
    v[425][4] = 29;
    v[426][4] = 9;
    v[427][4] = 3;
    v[428][4] = 21;
    v[429][4] = 31;
    v[430][4] = 29;
    v[431][4] = 25;
    v[432][4] = 1;
    v[433][4] = 3;
    v[434][4] = 9;
    v[435][4] = 27;
    v[436][4] = 5;
    v[437][4] = 27;
    v[438][4] = 25;
    v[439][4] = 21;
    v[440][4] = 11;
    v[441][4] = 29;
    v[442][4] = 31;
    v[443][4] = 27;
    v[444][4] = 21;
    v[445][4] = 29;
    v[446][4] = 17;
    v[447][4] = 9;
    v[448][4] = 17;
    v[449][4] = 13;
    v[450][4] = 11;
    v[451][4] = 25;
    v[452][4] = 15;
    v[453][4] = 21;
    v[454][4] = 11;
    v[455][4] = 19;
    v[456][4] = 31;
    v[457][4] = 3;
    v[458][4] = 19;
    v[459][4] = 5;
    v[460][4] = 3;
    v[461][4] = 3;
    v[462][4] = 9;
    v[463][4] = 13;
    v[464][4] = 13;
    v[465][4] = 3;
    v[466][4] = 29;
    v[467][4] = 7;
    v[468][4] = 5;
    v[469][4] = 9;
    v[470][4] = 23;
    v[471][4] = 13;
    v[472][4] = 21;
    v[473][4] = 23;
    v[474][4] = 21;
    v[475][4] = 31;
    v[476][4] = 11;
    v[477][4] = 7;
    v[478][4] = 7;
    v[479][4] = 3;
    v[480][4] = 23;
    v[481][4] = 1;
    v[482][4] = 23;
    v[483][4] = 5;
    v[484][4] = 9;
    v[485][4] = 17;
    v[486][4] = 21;
    v[487][4] = 1;
    v[488][4] = 17;
    v[489][4] = 29;
    v[490][4] = 7;
    v[491][4] = 5;
    v[492][4] = 17;
    v[493][4] = 13;
    v[494][4] = 25;
    v[495][4] = 17;
    v[496][4] = 9;
    v[497][4] = 19;
    v[498][4] = 9;
    v[499][4] = 5;
    v[500][4] = 7;
    v[501][4] = 21;
    v[502][4] = 19;
    v[503][4] = 13;
    v[504][4] = 9;
    v[505][4] = 7;
    v[506][4] = 3;
    v[507][4] = 9;
    v[508][4] = 3;
    v[509][4] = 15;
    v[510][4] = 31;
    v[511][4] = 29;
    v[512][4] = 29;
    v[513][4] = 25;
    v[514][4] = 13;
    v[515][4] = 9;
    v[516][4] = 21;
    v[517][4] = 9;
    v[518][4] = 31;
    v[519][4] = 7;
    v[520][4] = 15;
    v[521][4] = 5;
    v[522][4] = 31;
    v[523][4] = 7;
    v[524][4] = 15;
    v[525][4] = 27;
    v[526][4] = 25;
    v[527][4] = 19;
    v[528][4] = 9;
    v[529][4] = 9;
    v[530][4] = 25;
    v[531][4] = 25;
    v[532][4] = 23;
    v[533][4] = 1;
    v[534][4] = 9;
    v[535][4] = 7;
    v[536][4] = 11;
    v[537][4] = 15;
    v[538][4] = 19;
    v[539][4] = 15;
    v[540][4] = 27;
    v[541][4] = 17;
    v[542][4] = 11;
    v[543][4] = 11;
    v[544][4] = 31;
    v[545][4] = 13;
    v[546][4] = 25;
    v[547][4] = 25;
    v[548][4] = 9;
    v[549][4] = 7;
    v[550][4] = 13;
    v[551][4] = 29;
    v[552][4] = 19;
    v[553][4] = 5;
    v[554][4] = 19;
    v[555][4] = 31;
    v[556][4] = 25;
    v[557][4] = 13;
    v[558][4] = 25;
    v[559][4] = 15;
    v[560][4] = 5;
    v[561][4] = 9;
    v[562][4] = 29;
    v[563][4] = 31;
    v[564][4] = 9;
    v[565][4] = 29;
    v[566][4] = 27;
    v[567][4] = 25;
    v[568][4] = 27;
    v[569][4] = 11;
    v[570][4] = 17;
    v[571][4] = 5;
    v[572][4] = 17;
    v[573][4] = 3;
    v[574][4] = 23;
    v[575][4] = 15;
    v[576][4] = 9;
    v[577][4] = 9;
    v[578][4] = 17;
    v[579][4] = 17;
    v[580][4] = 31;
    v[581][4] = 11;
    v[582][4] = 19;
    v[583][4] = 25;
    v[584][4] = 13;
    v[585][4] = 23;
    v[586][4] = 15;
    v[587][4] = 25;
    v[588][4] = 21;
    v[589][4] = 31;
    v[590][4] = 19;
    v[591][4] = 3;
    v[592][4] = 11;
    v[593][4] = 25;
    v[594][4] = 7;
    v[595][4] = 15;
    v[596][4] = 19;
    v[597][4] = 7;
    v[598][4] = 5;
    v[599][4] = 3;
    v[600][4] = 13;
    v[601][4] = 13;
    v[602][4] = 1;
    v[603][4] = 23;
    v[604][4] = 5;
    v[605][4] = 25;
    v[606][4] = 11;
    v[607][4] = 25;
    v[608][4] = 15;
    v[609][4] = 13;
    v[610][4] = 21;
    v[611][4] = 11;
    v[612][4] = 23;
    v[613][4] = 29;
    v[614][4] = 5;
    v[615][4] = 17;
    v[616][4] = 27;
    v[617][4] = 9;
    v[618][4] = 19;
    v[619][4] = 15;
    v[620][4] = 5;
    v[621][4] = 29;
    v[622][4] = 23;
    v[623][4] = 19;
    v[624][4] = 1;
    v[625][4] = 27;
    v[626][4] = 3;
    v[627][4] = 23;
    v[628][4] = 21;
    v[629][4] = 19;
    v[630][4] = 27;
    v[631][4] = 11;
    v[632][4] = 17;
    v[633][4] = 13;
    v[634][4] = 27;
    v[635][4] = 11;
    v[636][4] = 31;
    v[637][4] = 23;
    v[638][4] = 5;
    v[639][4] = 9;
    v[640][4] = 21;
    v[641][4] = 31;
    v[642][4] = 29;
    v[643][4] = 11;
    v[644][4] = 21;
    v[645][4] = 17;
    v[646][4] = 15;
    v[647][4] = 7;
    v[648][4] = 15;
    v[649][4] = 7;
    v[650][4] = 9;
    v[651][4] = 21;
    v[652][4] = 27;
    v[653][4] = 25;
    v[654][4] = 29;
    v[655][4] = 11;
    v[656][4] = 3;
    v[657][4] = 21;
    v[658][4] = 13;
    v[659][4] = 23;
    v[660][4] = 19;
    v[661][4] = 27;
    v[662][4] = 17;
    v[663][4] = 29;
    v[664][4] = 25;
    v[665][4] = 17;
    v[666][4] = 9;
    v[667][4] = 1;
    v[668][4] = 19;
    v[669][4] = 23;
    v[670][4] = 5;
    v[671][4] = 23;
    v[672][4] = 1;
    v[673][4] = 17;
    v[674][4] = 17;
    v[675][4] = 13;
    v[676][4] = 27;
    v[677][4] = 23;
    v[678][4] = 7;
    v[679][4] = 7;
    v[680][4] = 11;
    v[681][4] = 13;
    v[682][4] = 17;
    v[683][4] = 13;
    v[684][4] = 11;
    v[685][4] = 21;
    v[686][4] = 13;
    v[687][4] = 23;
    v[688][4] = 1;
    v[689][4] = 27;
    v[690][4] = 13;
    v[691][4] = 9;
    v[692][4] = 7;
    v[693][4] = 1;
    v[694][4] = 27;
    v[695][4] = 29;
    v[696][4] = 5;
    v[697][4] = 13;
    v[698][4] = 25;
    v[699][4] = 21;
    v[700][4] = 3;
    v[701][4] = 31;
    v[702][4] = 15;
    v[703][4] = 13;
    v[704][4] = 3;
    v[705][4] = 19;
    v[706][4] = 13;
    v[707][4] = 1;
    v[708][4] = 27;
    v[709][4] = 15;
    v[710][4] = 17;
    v[711][4] = 1;
    v[712][4] = 3;
    v[713][4] = 13;
    v[714][4] = 13;
    v[715][4] = 13;
    v[716][4] = 31;
    v[717][4] = 29;
    v[718][4] = 27;
    v[719][4] = 7;
    v[720][4] = 7;
    v[721][4] = 21;
    v[722][4] = 29;
    v[723][4] = 15;
    v[724][4] = 17;
    v[725][4] = 17;
    v[726][4] = 21;
    v[727][4] = 19;
    v[728][4] = 17;
    v[729][4] = 3;
    v[730][4] = 15;
    v[731][4] = 5;
    v[732][4] = 27;
    v[733][4] = 27;
    v[734][4] = 3;
    v[735][4] = 31;
    v[736][4] = 31;
    v[737][4] = 7;
    v[738][4] = 21;
    v[739][4] = 3;
    v[740][4] = 13;
    v[741][4] = 11;
    v[742][4] = 17;
    v[743][4] = 27;
    v[744][4] = 25;
    v[745][4] = 1;
    v[746][4] = 9;
    v[747][4] = 7;
    v[748][4] = 29;
    v[749][4] = 27;
    v[750][4] = 21;
    v[751][4] = 23;
    v[752][4] = 13;
    v[753][4] = 25;
    v[754][4] = 29;
    v[755][4] = 15;
    v[756][4] = 17;
    v[757][4] = 29;
    v[758][4] = 9;
    v[759][4] = 15;
    v[760][4] = 3;
    v[761][4] = 21;
    v[762][4] = 15;
    v[763][4] = 17;
    v[764][4] = 17;
    v[765][4] = 31;
    v[766][4] = 9;
    v[767][4] = 9;
    v[768][4] = 23;
    v[769][4] = 19;
    v[770][4] = 25;
    v[771][4] = 3;
    v[772][4] = 1;
    v[773][4] = 11;
    v[774][4] = 27;
    v[775][4] = 29;
    v[776][4] = 1;
    v[777][4] = 31;
    v[778][4] = 29;
    v[779][4] = 25;
    v[780][4] = 29;
    v[781][4] = 1;
    v[782][4] = 23;
    v[783][4] = 29;
    v[784][4] = 25;
    v[785][4] = 13;
    v[786][4] = 3;
    v[787][4] = 31;
    v[788][4] = 25;
    v[789][4] = 5;
    v[790][4] = 5;
    v[791][4] = 11;
    v[792][4] = 3;
    v[793][4] = 21;
    v[794][4] = 9;
    v[795][4] = 23;
    v[796][4] = 7;
    v[797][4] = 11;
    v[798][4] = 23;
    v[799][4] = 11;
    v[800][4] = 1;
    v[801][4] = 1;
    v[802][4] = 3;
    v[803][4] = 23;
    v[804][4] = 25;
    v[805][4] = 23;
    v[806][4] = 1;
    v[807][4] = 23;
    v[808][4] = 3;
    v[809][4] = 27;
    v[810][4] = 9;
    v[811][4] = 27;
    v[812][4] = 3;
    v[813][4] = 23;
    v[814][4] = 25;
    v[815][4] = 19;
    v[816][4] = 29;
    v[817][4] = 29;
    v[818][4] = 13;
    v[819][4] = 27;
    v[820][4] = 5;
    v[821][4] = 9;
    v[822][4] = 29;
    v[823][4] = 29;
    v[824][4] = 13;
    v[825][4] = 17;
    v[826][4] = 3;
    v[827][4] = 23;
    v[828][4] = 19;
    v[829][4] = 7;
    v[830][4] = 13;
    v[831][4] = 3;
    v[832][4] = 19;
    v[833][4] = 23;
    v[834][4] = 5;
    v[835][4] = 29;
    v[836][4] = 29;
    v[837][4] = 13;
    v[838][4] = 13;
    v[839][4] = 5;
    v[840][4] = 19;
    v[841][4] = 5;
    v[842][4] = 17;
    v[843][4] = 9;
    v[844][4] = 11;
    v[845][4] = 11;
    v[846][4] = 29;
    v[847][4] = 27;
    v[848][4] = 23;
    v[849][4] = 19;
    v[850][4] = 17;
    v[851][4] = 25;
    v[852][4] = 13;
    v[853][4] = 1;
    v[854][4] = 13;
    v[855][4] = 3;
    v[856][4] = 11;
    v[857][4] = 1;
    v[858][4] = 17;
    v[859][4] = 29;
    v[860][4] = 1;
    v[861][4] = 13;
    v[862][4] = 17;
    v[863][4] = 9;
    v[864][4] = 17;
    v[865][4] = 21;
    v[866][4] = 1;
    v[867][4] = 11;
    v[868][4] = 1;
    v[869][4] = 1;
    v[870][4] = 25;
    v[871][4] = 5;
    v[872][4] = 7;
    v[873][4] = 29;
    v[874][4] = 29;
    v[875][4] = 19;
    v[876][4] = 19;
    v[877][4] = 1;
    v[878][4] = 29;
    v[879][4] = 13;
    v[880][4] = 3;
    v[881][4] = 1;
    v[882][4] = 31;
    v[883][4] = 15;
    v[884][4] = 13;
    v[885][4] = 3;
    v[886][4] = 1;
    v[887][4] = 11;
    v[888][4] = 19;
    v[889][4] = 5;
    v[890][4] = 29;
    v[891][4] = 13;
    v[892][4] = 29;
    v[893][4] = 23;
    v[894][4] = 3;
    v[895][4] = 1;
    v[896][4] = 31;
    v[897][4] = 13;
    v[898][4] = 19;
    v[899][4] = 17;
    v[900][4] = 5;
    v[901][4] = 5;
    v[902][4] = 1;
    v[903][4] = 29;
    v[904][4] = 23;
    v[905][4] = 3;
    v[906][4] = 19;
    v[907][4] = 25;
    v[908][4] = 19;
    v[909][4] = 27;
    v[910][4] = 9;
    v[911][4] = 27;
    v[912][4] = 13;
    v[913][4] = 15;
    v[914][4] = 29;
    v[915][4] = 23;
    v[916][4] = 13;
    v[917][4] = 25;
    v[918][4] = 25;
    v[919][4] = 17;
    v[920][4] = 19;
    v[921][4] = 17;
    v[922][4] = 15;
    v[923][4] = 27;
    v[924][4] = 3;
    v[925][4] = 25;
    v[926][4] = 17;
    v[927][4] = 27;
    v[928][4] = 3;
    v[929][4] = 27;
    v[930][4] = 31;
    v[931][4] = 23;
    v[932][4] = 13;
    v[933][4] = 31;
    v[934][4] = 11;
    v[935][4] = 15;
    v[936][4] = 7;
    v[937][4] = 21;
    v[938][4] = 19;
    v[939][4] = 27;
    v[940][4] = 19;
    v[941][4] = 21;
    v[942][4] = 29;
    v[943][4] = 7;
    v[944][4] = 31;
    v[945][4] = 13;
    v[946][4] = 9;
    v[947][4] = 9;
    v[948][4] = 7;
    v[949][4] = 21;
    v[950][4] = 13;
    v[951][4] = 11;
    v[952][4] = 9;
    v[953][4] = 11;
    v[954][4] = 29;
    v[955][4] = 19;
    v[956][4] = 11;
    v[957][4] = 19;
    v[958][4] = 21;
    v[959][4] = 5;
    v[960][4] = 29;
    v[961][4] = 13;
    v[962][4] = 7;
    v[963][4] = 19;
    v[964][4] = 19;
    v[965][4] = 27;
    v[966][4] = 23;
    v[967][4] = 31;
    v[968][4] = 1;
    v[969][4] = 27;
    v[970][4] = 21;
    v[971][4] = 7;
    v[972][4] = 3;
    v[973][4] = 7;
    v[974][4] = 11;
    v[975][4] = 23;
    v[976][4] = 13;
    v[977][4] = 29;
    v[978][4] = 11;
    v[979][4] = 31;
    v[980][4] = 19;
    v[981][4] = 1;
    v[982][4] = 5;
    v[983][4] = 5;
    v[984][4] = 11;
    v[985][4] = 5;
    v[986][4] = 3;
    v[987][4] = 27;
    v[988][4] = 5;
    v[989][4] = 7;
    v[990][4] = 11;
    v[991][4] = 31;
    v[992][4] = 1;
    v[993][4] = 27;
    v[994][4] = 31;
    v[995][4] = 31;
    v[996][4] = 23;
    v[997][4] = 5;
    v[998][4] = 21;
    v[999][4] = 27;
    v[1000][4] = 9;
    v[1001][4] = 25;
    v[1002][4] = 3;
    v[1003][4] = 15;
    v[1004][4] = 19;
    v[1005][4] = 1;
    v[1006][4] = 19;
    v[1007][4] = 9;
    v[1008][4] = 5;
    v[1009][4] = 25;
    v[1010][4] = 21;
    v[1011][4] = 15;
    v[1012][4] = 25;
    v[1013][4] = 29;
    v[1014][4] = 15;
    v[1015][4] = 21;
    v[1016][4] = 11;
    v[1017][4] = 19;
    v[1018][4] = 15;
    v[1019][4] = 3;
    v[1020][4] = 7;
    v[1021][4] = 13;
    v[1022][4] = 11;
    v[1023][4] = 25;
    v[1024][4] = 17;
    v[1025][4] = 1;
    v[1026][4] = 5;
    v[1027][4] = 31;
    v[1028][4] = 13;
    v[1029][4] = 29;
    v[1030][4] = 23;
    v[1031][4] = 9;
    v[1032][4] = 5;
    v[1033][4] = 29;
    v[1034][4] = 7;
    v[1035][4] = 17;
    v[1036][4] = 27;
    v[1037][4] = 7;
    v[1038][4] = 17;
    v[1039][4] = 31;
    v[1040][4] = 9;
    v[1041][4] = 31;
    v[1042][4] = 9;
    v[1043][4] = 9;
    v[1044][4] = 7;
    v[1045][4] = 21;
    v[1046][4] = 3;
    v[1047][4] = 3;
    v[1048][4] = 3;
    v[1049][4] = 9;
    v[1050][4] = 11;
    v[1051][4] = 21;
    v[1052][4] = 11;
    v[1053][4] = 31;
    v[1054][4] = 9;
    v[1055][4] = 25;
    v[1056][4] = 5;
    v[1057][4] = 1;
    v[1058][4] = 31;
    v[1059][4] = 13;
    v[1060][4] = 29;
    v[1061][4] = 9;
    v[1062][4] = 29;
    v[1063][4] = 1;
    v[1064][4] = 11;
    v[1065][4] = 19;
    v[1066][4] = 7;
    v[1067][4] = 27;
    v[1068][4] = 13;
    v[1069][4] = 31;
    v[1070][4] = 7;
    v[1071][4] = 31;
    v[1072][4] = 7;
    v[1073][4] = 25;
    v[1074][4] = 23;
    v[1075][4] = 21;
    v[1076][4] = 29;
    v[1077][4] = 11;
    v[1078][4] = 11;
    v[1079][4] = 13;
    v[1080][4] = 11;
    v[1081][4] = 27;
    v[1082][4] = 1;
    v[1083][4] = 23;
    v[1084][4] = 31;
    v[1085][4] = 21;
    v[1086][4] = 23;
    v[1087][4] = 21;
    v[1088][4] = 19;
    v[1089][4] = 31;
    v[1090][4] = 5;
    v[1091][4] = 31;
    v[1092][4] = 25;
    v[1093][4] = 25;
    v[1094][4] = 19;
    v[1095][4] = 17;
    v[1096][4] = 11;
    v[1097][4] = 25;
    v[1098][4] = 7;
    v[1099][4] = 13;
    v[1100][4] = 1;
    v[1101][4] = 29;
    v[1102][4] = 17;
    v[1103][4] = 23;
    v[1104][4] = 15;
    v[1105][4] = 7;
    v[1106][4] = 29;
    v[1107][4] = 17;
    v[1108][4] = 13;
    v[1109][4] = 3;
    v[1110][4] = 17;

    v[13][5] = 37;
    v[14][5] = 33;
    v[15][5] = 7;
    v[16][5] = 5;
    v[17][5] = 11;
    v[18][5] = 39;
    v[19][5] = 63;
    v[20][5] = 59;
    v[21][5] = 17;
    v[22][5] = 15;
    v[23][5] = 23;
    v[24][5] = 29;
    v[25][5] = 3;
    v[26][5] = 21;
    v[27][5] = 13;
    v[28][5] = 31;
    v[29][5] = 25;
    v[30][5] = 9;
    v[31][5] = 49;
    v[32][5] = 33;
    v[33][5] = 19;
    v[34][5] = 29;
    v[35][5] = 11;
    v[36][5] = 19;
    v[37][5] = 27;
    v[38][5] = 15;
    v[39][5] = 25;
    v[40][5] = 63;
    v[41][5] = 55;
    v[42][5] = 17;
    v[43][5] = 63;
    v[44][5] = 49;
    v[45][5] = 19;
    v[46][5] = 41;
    v[47][5] = 59;
    v[48][5] = 3;
    v[49][5] = 57;
    v[50][5] = 33;
    v[51][5] = 49;
    v[52][5] = 53;
    v[53][5] = 57;
    v[54][5] = 57;
    v[55][5] = 39;
    v[56][5] = 21;
    v[57][5] = 7;
    v[58][5] = 53;
    v[59][5] = 9;
    v[60][5] = 55;
    v[61][5] = 15;
    v[62][5] = 59;
    v[63][5] = 19;
    v[64][5] = 49;
    v[65][5] = 31;
    v[66][5] = 3;
    v[67][5] = 39;
    v[68][5] = 5;
    v[69][5] = 5;
    v[70][5] = 41;
    v[71][5] = 9;
    v[72][5] = 19;
    v[73][5] = 9;
    v[74][5] = 57;
    v[75][5] = 25;
    v[76][5] = 1;
    v[77][5] = 15;
    v[78][5] = 51;
    v[79][5] = 11;
    v[80][5] = 19;
    v[81][5] = 61;
    v[82][5] = 53;
    v[83][5] = 29;
    v[84][5] = 19;
    v[85][5] = 11;
    v[86][5] = 9;
    v[87][5] = 21;
    v[88][5] = 19;
    v[89][5] = 43;
    v[90][5] = 13;
    v[91][5] = 13;
    v[92][5] = 41;
    v[93][5] = 25;
    v[94][5] = 31;
    v[95][5] = 9;
    v[96][5] = 11;
    v[97][5] = 19;
    v[98][5] = 5;
    v[99][5] = 53;
    v[100][5] = 37;
    v[101][5] = 7;
    v[102][5] = 51;
    v[103][5] = 45;
    v[104][5] = 7;
    v[105][5] = 7;
    v[106][5] = 61;
    v[107][5] = 23;
    v[108][5] = 45;
    v[109][5] = 7;
    v[110][5] = 59;
    v[111][5] = 41;
    v[112][5] = 1;
    v[113][5] = 29;
    v[114][5] = 61;
    v[115][5] = 37;
    v[116][5] = 27;
    v[117][5] = 47;
    v[118][5] = 15;
    v[119][5] = 31;
    v[120][5] = 35;
    v[121][5] = 31;
    v[122][5] = 17;
    v[123][5] = 51;
    v[124][5] = 13;
    v[125][5] = 25;
    v[126][5] = 45;
    v[127][5] = 5;
    v[128][5] = 5;
    v[129][5] = 33;
    v[130][5] = 39;
    v[131][5] = 5;
    v[132][5] = 47;
    v[133][5] = 29;
    v[134][5] = 35;
    v[135][5] = 47;
    v[136][5] = 63;
    v[137][5] = 45;
    v[138][5] = 37;
    v[139][5] = 47;
    v[140][5] = 59;
    v[141][5] = 21;
    v[142][5] = 59;
    v[143][5] = 33;
    v[144][5] = 51;
    v[145][5] = 9;
    v[146][5] = 27;
    v[147][5] = 13;
    v[148][5] = 25;
    v[149][5] = 43;
    v[150][5] = 3;
    v[151][5] = 17;
    v[152][5] = 21;
    v[153][5] = 59;
    v[154][5] = 61;
    v[155][5] = 27;
    v[156][5] = 47;
    v[157][5] = 57;
    v[158][5] = 11;
    v[159][5] = 17;
    v[160][5] = 39;
    v[161][5] = 1;
    v[162][5] = 63;
    v[163][5] = 21;
    v[164][5] = 59;
    v[165][5] = 17;
    v[166][5] = 13;
    v[167][5] = 31;
    v[168][5] = 3;
    v[169][5] = 31;
    v[170][5] = 7;
    v[171][5] = 9;
    v[172][5] = 27;
    v[173][5] = 37;
    v[174][5] = 23;
    v[175][5] = 31;
    v[176][5] = 9;
    v[177][5] = 45;
    v[178][5] = 43;
    v[179][5] = 31;
    v[180][5] = 63;
    v[181][5] = 21;
    v[182][5] = 39;
    v[183][5] = 51;
    v[184][5] = 27;
    v[185][5] = 7;
    v[186][5] = 53;
    v[187][5] = 11;
    v[188][5] = 1;
    v[189][5] = 59;
    v[190][5] = 39;
    v[191][5] = 23;
    v[192][5] = 49;
    v[193][5] = 23;
    v[194][5] = 7;
    v[195][5] = 55;
    v[196][5] = 59;
    v[197][5] = 3;
    v[198][5] = 19;
    v[199][5] = 35;
    v[200][5] = 13;
    v[201][5] = 9;
    v[202][5] = 13;
    v[203][5] = 15;
    v[204][5] = 23;
    v[205][5] = 9;
    v[206][5] = 7;
    v[207][5] = 43;
    v[208][5] = 55;
    v[209][5] = 3;
    v[210][5] = 19;
    v[211][5] = 9;
    v[212][5] = 27;
    v[213][5] = 33;
    v[214][5] = 27;
    v[215][5] = 49;
    v[216][5] = 23;
    v[217][5] = 47;
    v[218][5] = 19;
    v[219][5] = 7;
    v[220][5] = 11;
    v[221][5] = 55;
    v[222][5] = 27;
    v[223][5] = 35;
    v[224][5] = 5;
    v[225][5] = 5;
    v[226][5] = 55;
    v[227][5] = 35;
    v[228][5] = 37;
    v[229][5] = 9;
    v[230][5] = 33;
    v[231][5] = 29;
    v[232][5] = 47;
    v[233][5] = 25;
    v[234][5] = 11;
    v[235][5] = 47;
    v[236][5] = 53;
    v[237][5] = 61;
    v[238][5] = 59;
    v[239][5] = 3;
    v[240][5] = 53;
    v[241][5] = 47;
    v[242][5] = 5;
    v[243][5] = 19;
    v[244][5] = 59;
    v[245][5] = 5;
    v[246][5] = 47;
    v[247][5] = 23;
    v[248][5] = 45;
    v[249][5] = 53;
    v[250][5] = 3;
    v[251][5] = 49;
    v[252][5] = 61;
    v[253][5] = 47;
    v[254][5] = 39;
    v[255][5] = 29;
    v[256][5] = 17;
    v[257][5] = 57;
    v[258][5] = 5;
    v[259][5] = 17;
    v[260][5] = 31;
    v[261][5] = 23;
    v[262][5] = 41;
    v[263][5] = 39;
    v[264][5] = 5;
    v[265][5] = 27;
    v[266][5] = 7;
    v[267][5] = 29;
    v[268][5] = 29;
    v[269][5] = 33;
    v[270][5] = 31;
    v[271][5] = 41;
    v[272][5] = 31;
    v[273][5] = 29;
    v[274][5] = 17;
    v[275][5] = 29;
    v[276][5] = 29;
    v[277][5] = 9;
    v[278][5] = 9;
    v[279][5] = 31;
    v[280][5] = 27;
    v[281][5] = 53;
    v[282][5] = 35;
    v[283][5] = 5;
    v[284][5] = 61;
    v[285][5] = 1;
    v[286][5] = 49;
    v[287][5] = 13;
    v[288][5] = 57;
    v[289][5] = 29;
    v[290][5] = 5;
    v[291][5] = 21;
    v[292][5] = 43;
    v[293][5] = 25;
    v[294][5] = 57;
    v[295][5] = 49;
    v[296][5] = 37;
    v[297][5] = 27;
    v[298][5] = 11;
    v[299][5] = 61;
    v[300][5] = 37;
    v[301][5] = 49;
    v[302][5] = 5;
    v[303][5] = 63;
    v[304][5] = 63;
    v[305][5] = 3;
    v[306][5] = 45;
    v[307][5] = 37;
    v[308][5] = 63;
    v[309][5] = 21;
    v[310][5] = 21;
    v[311][5] = 19;
    v[312][5] = 27;
    v[313][5] = 59;
    v[314][5] = 21;
    v[315][5] = 45;
    v[316][5] = 23;
    v[317][5] = 13;
    v[318][5] = 15;
    v[319][5] = 3;
    v[320][5] = 43;
    v[321][5] = 63;
    v[322][5] = 39;
    v[323][5] = 19;
    v[324][5] = 63;
    v[325][5] = 31;
    v[326][5] = 41;
    v[327][5] = 41;
    v[328][5] = 15;
    v[329][5] = 43;
    v[330][5] = 63;
    v[331][5] = 53;
    v[332][5] = 1;
    v[333][5] = 63;
    v[334][5] = 31;
    v[335][5] = 7;
    v[336][5] = 17;
    v[337][5] = 11;
    v[338][5] = 61;
    v[339][5] = 31;
    v[340][5] = 51;
    v[341][5] = 37;
    v[342][5] = 29;
    v[343][5] = 59;
    v[344][5] = 25;
    v[345][5] = 63;
    v[346][5] = 59;
    v[347][5] = 47;
    v[348][5] = 15;
    v[349][5] = 27;
    v[350][5] = 19;
    v[351][5] = 29;
    v[352][5] = 45;
    v[353][5] = 35;
    v[354][5] = 55;
    v[355][5] = 39;
    v[356][5] = 19;
    v[357][5] = 43;
    v[358][5] = 21;
    v[359][5] = 19;
    v[360][5] = 13;
    v[361][5] = 17;
    v[362][5] = 51;
    v[363][5] = 37;
    v[364][5] = 5;
    v[365][5] = 33;
    v[366][5] = 35;
    v[367][5] = 49;
    v[368][5] = 25;
    v[369][5] = 45;
    v[370][5] = 1;
    v[371][5] = 63;
    v[372][5] = 47;
    v[373][5] = 9;
    v[374][5] = 63;
    v[375][5] = 15;
    v[376][5] = 25;
    v[377][5] = 25;
    v[378][5] = 15;
    v[379][5] = 41;
    v[380][5] = 13;
    v[381][5] = 3;
    v[382][5] = 19;
    v[383][5] = 51;
    v[384][5] = 49;
    v[385][5] = 37;
    v[386][5] = 25;
    v[387][5] = 49;
    v[388][5] = 13;
    v[389][5] = 53;
    v[390][5] = 47;
    v[391][5] = 23;
    v[392][5] = 35;
    v[393][5] = 29;
    v[394][5] = 33;
    v[395][5] = 21;
    v[396][5] = 35;
    v[397][5] = 23;
    v[398][5] = 3;
    v[399][5] = 43;
    v[400][5] = 31;
    v[401][5] = 63;
    v[402][5] = 9;
    v[403][5] = 1;
    v[404][5] = 61;
    v[405][5] = 43;
    v[406][5] = 3;
    v[407][5] = 11;
    v[408][5] = 55;
    v[409][5] = 11;
    v[410][5] = 35;
    v[411][5] = 1;
    v[412][5] = 63;
    v[413][5] = 35;
    v[414][5] = 49;
    v[415][5] = 19;
    v[416][5] = 45;
    v[417][5] = 9;
    v[418][5] = 57;
    v[419][5] = 51;
    v[420][5] = 1;
    v[421][5] = 47;
    v[422][5] = 41;
    v[423][5] = 9;
    v[424][5] = 11;
    v[425][5] = 37;
    v[426][5] = 19;
    v[427][5] = 55;
    v[428][5] = 23;
    v[429][5] = 55;
    v[430][5] = 55;
    v[431][5] = 13;
    v[432][5] = 7;
    v[433][5] = 47;
    v[434][5] = 37;
    v[435][5] = 11;
    v[436][5] = 43;
    v[437][5] = 17;
    v[438][5] = 3;
    v[439][5] = 25;
    v[440][5] = 19;
    v[441][5] = 55;
    v[442][5] = 59;
    v[443][5] = 37;
    v[444][5] = 33;
    v[445][5] = 43;
    v[446][5] = 1;
    v[447][5] = 5;
    v[448][5] = 21;
    v[449][5] = 5;
    v[450][5] = 63;
    v[451][5] = 49;
    v[452][5] = 61;
    v[453][5] = 21;
    v[454][5] = 51;
    v[455][5] = 15;
    v[456][5] = 19;
    v[457][5] = 43;
    v[458][5] = 47;
    v[459][5] = 17;
    v[460][5] = 9;
    v[461][5] = 53;
    v[462][5] = 45;
    v[463][5] = 11;
    v[464][5] = 51;
    v[465][5] = 25;
    v[466][5] = 11;
    v[467][5] = 25;
    v[468][5] = 47;
    v[469][5] = 47;
    v[470][5] = 1;
    v[471][5] = 43;
    v[472][5] = 29;
    v[473][5] = 17;
    v[474][5] = 31;
    v[475][5] = 15;
    v[476][5] = 59;
    v[477][5] = 27;
    v[478][5] = 63;
    v[479][5] = 11;
    v[480][5] = 41;
    v[481][5] = 51;
    v[482][5] = 29;
    v[483][5] = 7;
    v[484][5] = 27;
    v[485][5] = 63;
    v[486][5] = 31;
    v[487][5] = 43;
    v[488][5] = 3;
    v[489][5] = 29;
    v[490][5] = 39;
    v[491][5] = 3;
    v[492][5] = 59;
    v[493][5] = 59;
    v[494][5] = 1;
    v[495][5] = 53;
    v[496][5] = 63;
    v[497][5] = 23;
    v[498][5] = 63;
    v[499][5] = 47;
    v[500][5] = 51;
    v[501][5] = 23;
    v[502][5] = 61;
    v[503][5] = 39;
    v[504][5] = 47;
    v[505][5] = 21;
    v[506][5] = 39;
    v[507][5] = 15;
    v[508][5] = 3;
    v[509][5] = 9;
    v[510][5] = 57;
    v[511][5] = 61;
    v[512][5] = 39;
    v[513][5] = 37;
    v[514][5] = 21;
    v[515][5] = 51;
    v[516][5] = 1;
    v[517][5] = 23;
    v[518][5] = 43;
    v[519][5] = 27;
    v[520][5] = 25;
    v[521][5] = 11;
    v[522][5] = 13;
    v[523][5] = 21;
    v[524][5] = 43;
    v[525][5] = 7;
    v[526][5] = 11;
    v[527][5] = 33;
    v[528][5] = 55;
    v[529][5] = 1;
    v[530][5] = 37;
    v[531][5] = 35;
    v[532][5] = 27;
    v[533][5] = 61;
    v[534][5] = 39;
    v[535][5] = 5;
    v[536][5] = 19;
    v[537][5] = 61;
    v[538][5] = 61;
    v[539][5] = 57;
    v[540][5] = 59;
    v[541][5] = 21;
    v[542][5] = 59;
    v[543][5] = 61;
    v[544][5] = 57;
    v[545][5] = 25;
    v[546][5] = 55;
    v[547][5] = 27;
    v[548][5] = 31;
    v[549][5] = 41;
    v[550][5] = 33;
    v[551][5] = 63;
    v[552][5] = 19;
    v[553][5] = 57;
    v[554][5] = 35;
    v[555][5] = 13;
    v[556][5] = 63;
    v[557][5] = 35;
    v[558][5] = 17;
    v[559][5] = 11;
    v[560][5] = 11;
    v[561][5] = 49;
    v[562][5] = 41;
    v[563][5] = 55;
    v[564][5] = 5;
    v[565][5] = 45;
    v[566][5] = 17;
    v[567][5] = 35;
    v[568][5] = 5;
    v[569][5] = 31;
    v[570][5] = 31;
    v[571][5] = 37;
    v[572][5] = 17;
    v[573][5] = 45;
    v[574][5] = 51;
    v[575][5] = 1;
    v[576][5] = 39;
    v[577][5] = 49;
    v[578][5] = 55;
    v[579][5] = 19;
    v[580][5] = 41;
    v[581][5] = 13;
    v[582][5] = 5;
    v[583][5] = 51;
    v[584][5] = 5;
    v[585][5] = 49;
    v[586][5] = 1;
    v[587][5] = 21;
    v[588][5] = 13;
    v[589][5] = 17;
    v[590][5] = 59;
    v[591][5] = 51;
    v[592][5] = 11;
    v[593][5] = 3;
    v[594][5] = 61;
    v[595][5] = 1;
    v[596][5] = 33;
    v[597][5] = 37;
    v[598][5] = 33;
    v[599][5] = 61;
    v[600][5] = 25;
    v[601][5] = 27;
    v[602][5] = 59;
    v[603][5] = 7;
    v[604][5] = 49;
    v[605][5] = 13;
    v[606][5] = 63;
    v[607][5] = 3;
    v[608][5] = 33;
    v[609][5] = 3;
    v[610][5] = 15;
    v[611][5] = 9;
    v[612][5] = 13;
    v[613][5] = 35;
    v[614][5] = 39;
    v[615][5] = 11;
    v[616][5] = 59;
    v[617][5] = 59;
    v[618][5] = 1;
    v[619][5] = 57;
    v[620][5] = 11;
    v[621][5] = 5;
    v[622][5] = 57;
    v[623][5] = 13;
    v[624][5] = 31;
    v[625][5] = 13;
    v[626][5] = 11;
    v[627][5] = 55;
    v[628][5] = 45;
    v[629][5] = 9;
    v[630][5] = 55;
    v[631][5] = 55;
    v[632][5] = 19;
    v[633][5] = 25;
    v[634][5] = 41;
    v[635][5] = 23;
    v[636][5] = 45;
    v[637][5] = 29;
    v[638][5] = 63;
    v[639][5] = 59;
    v[640][5] = 27;
    v[641][5] = 39;
    v[642][5] = 21;
    v[643][5] = 37;
    v[644][5] = 7;
    v[645][5] = 61;
    v[646][5] = 49;
    v[647][5] = 35;
    v[648][5] = 39;
    v[649][5] = 9;
    v[650][5] = 29;
    v[651][5] = 7;
    v[652][5] = 25;
    v[653][5] = 23;
    v[654][5] = 57;
    v[655][5] = 5;
    v[656][5] = 19;
    v[657][5] = 15;
    v[658][5] = 33;
    v[659][5] = 49;
    v[660][5] = 37;
    v[661][5] = 25;
    v[662][5] = 17;
    v[663][5] = 45;
    v[664][5] = 29;
    v[665][5] = 15;
    v[666][5] = 25;
    v[667][5] = 3;
    v[668][5] = 3;
    v[669][5] = 49;
    v[670][5] = 11;
    v[671][5] = 39;
    v[672][5] = 15;
    v[673][5] = 19;
    v[674][5] = 57;
    v[675][5] = 39;
    v[676][5] = 15;
    v[677][5] = 11;
    v[678][5] = 3;
    v[679][5] = 57;
    v[680][5] = 31;
    v[681][5] = 55;
    v[682][5] = 61;
    v[683][5] = 19;
    v[684][5] = 5;
    v[685][5] = 41;
    v[686][5] = 35;
    v[687][5] = 59;
    v[688][5] = 61;
    v[689][5] = 39;
    v[690][5] = 41;
    v[691][5] = 53;
    v[692][5] = 53;
    v[693][5] = 63;
    v[694][5] = 31;
    v[695][5] = 9;
    v[696][5] = 59;
    v[697][5] = 13;
    v[698][5] = 35;
    v[699][5] = 55;
    v[700][5] = 41;
    v[701][5] = 49;
    v[702][5] = 5;
    v[703][5] = 41;
    v[704][5] = 25;
    v[705][5] = 27;
    v[706][5] = 43;
    v[707][5] = 5;
    v[708][5] = 5;
    v[709][5] = 43;
    v[710][5] = 5;
    v[711][5] = 5;
    v[712][5] = 17;
    v[713][5] = 5;
    v[714][5] = 15;
    v[715][5] = 27;
    v[716][5] = 29;
    v[717][5] = 17;
    v[718][5] = 9;
    v[719][5] = 3;
    v[720][5] = 55;
    v[721][5] = 31;
    v[722][5] = 1;
    v[723][5] = 45;
    v[724][5] = 45;
    v[725][5] = 13;
    v[726][5] = 57;
    v[727][5] = 17;
    v[728][5] = 3;
    v[729][5] = 61;
    v[730][5] = 15;
    v[731][5] = 49;
    v[732][5] = 15;
    v[733][5] = 47;
    v[734][5] = 9;
    v[735][5] = 37;
    v[736][5] = 45;
    v[737][5] = 9;
    v[738][5] = 51;
    v[739][5] = 61;
    v[740][5] = 21;
    v[741][5] = 33;
    v[742][5] = 11;
    v[743][5] = 21;
    v[744][5] = 63;
    v[745][5] = 63;
    v[746][5] = 47;
    v[747][5] = 57;
    v[748][5] = 61;
    v[749][5] = 49;
    v[750][5] = 9;
    v[751][5] = 59;
    v[752][5] = 19;
    v[753][5] = 29;
    v[754][5] = 21;
    v[755][5] = 23;
    v[756][5] = 55;
    v[757][5] = 23;
    v[758][5] = 43;
    v[759][5] = 41;
    v[760][5] = 57;
    v[761][5] = 9;
    v[762][5] = 39;
    v[763][5] = 27;
    v[764][5] = 41;
    v[765][5] = 35;
    v[766][5] = 61;
    v[767][5] = 29;
    v[768][5] = 57;
    v[769][5] = 63;
    v[770][5] = 21;
    v[771][5] = 31;
    v[772][5] = 59;
    v[773][5] = 35;
    v[774][5] = 49;
    v[775][5] = 3;
    v[776][5] = 49;
    v[777][5] = 47;
    v[778][5] = 49;
    v[779][5] = 33;
    v[780][5] = 21;
    v[781][5] = 19;
    v[782][5] = 21;
    v[783][5] = 35;
    v[784][5] = 11;
    v[785][5] = 17;
    v[786][5] = 37;
    v[787][5] = 23;
    v[788][5] = 59;
    v[789][5] = 13;
    v[790][5] = 37;
    v[791][5] = 35;
    v[792][5] = 55;
    v[793][5] = 57;
    v[794][5] = 1;
    v[795][5] = 29;
    v[796][5] = 45;
    v[797][5] = 11;
    v[798][5] = 1;
    v[799][5] = 15;
    v[800][5] = 9;
    v[801][5] = 33;
    v[802][5] = 19;
    v[803][5] = 53;
    v[804][5] = 43;
    v[805][5] = 39;
    v[806][5] = 23;
    v[807][5] = 7;
    v[808][5] = 13;
    v[809][5] = 13;
    v[810][5] = 1;
    v[811][5] = 19;
    v[812][5] = 41;
    v[813][5] = 55;
    v[814][5] = 1;
    v[815][5] = 13;
    v[816][5] = 15;
    v[817][5] = 59;
    v[818][5] = 55;
    v[819][5] = 15;
    v[820][5] = 3;
    v[821][5] = 57;
    v[822][5] = 37;
    v[823][5] = 31;
    v[824][5] = 17;
    v[825][5] = 1;
    v[826][5] = 3;
    v[827][5] = 21;
    v[828][5] = 29;
    v[829][5] = 25;
    v[830][5] = 55;
    v[831][5] = 9;
    v[832][5] = 37;
    v[833][5] = 33;
    v[834][5] = 53;
    v[835][5] = 41;
    v[836][5] = 51;
    v[837][5] = 19;
    v[838][5] = 57;
    v[839][5] = 13;
    v[840][5] = 63;
    v[841][5] = 43;
    v[842][5] = 19;
    v[843][5] = 7;
    v[844][5] = 13;
    v[845][5] = 37;
    v[846][5] = 33;
    v[847][5] = 19;
    v[848][5] = 15;
    v[849][5] = 63;
    v[850][5] = 51;
    v[851][5] = 11;
    v[852][5] = 49;
    v[853][5] = 23;
    v[854][5] = 57;
    v[855][5] = 47;
    v[856][5] = 51;
    v[857][5] = 15;
    v[858][5] = 53;
    v[859][5] = 41;
    v[860][5] = 1;
    v[861][5] = 15;
    v[862][5] = 37;
    v[863][5] = 61;
    v[864][5] = 11;
    v[865][5] = 35;
    v[866][5] = 29;
    v[867][5] = 33;
    v[868][5] = 23;
    v[869][5] = 55;
    v[870][5] = 11;
    v[871][5] = 59;
    v[872][5] = 19;
    v[873][5] = 61;
    v[874][5] = 61;
    v[875][5] = 45;
    v[876][5] = 13;
    v[877][5] = 49;
    v[878][5] = 13;
    v[879][5] = 63;
    v[880][5] = 5;
    v[881][5] = 61;
    v[882][5] = 5;
    v[883][5] = 31;
    v[884][5] = 17;
    v[885][5] = 61;
    v[886][5] = 63;
    v[887][5] = 13;
    v[888][5] = 27;
    v[889][5] = 57;
    v[890][5] = 1;
    v[891][5] = 21;
    v[892][5] = 5;
    v[893][5] = 11;
    v[894][5] = 39;
    v[895][5] = 57;
    v[896][5] = 51;
    v[897][5] = 53;
    v[898][5] = 39;
    v[899][5] = 25;
    v[900][5] = 41;
    v[901][5] = 39;
    v[902][5] = 37;
    v[903][5] = 23;
    v[904][5] = 31;
    v[905][5] = 25;
    v[906][5] = 33;
    v[907][5] = 17;
    v[908][5] = 57;
    v[909][5] = 29;
    v[910][5] = 27;
    v[911][5] = 23;
    v[912][5] = 47;
    v[913][5] = 41;
    v[914][5] = 29;
    v[915][5] = 19;
    v[916][5] = 47;
    v[917][5] = 41;
    v[918][5] = 25;
    v[919][5] = 5;
    v[920][5] = 51;
    v[921][5] = 43;
    v[922][5] = 39;
    v[923][5] = 29;
    v[924][5] = 7;
    v[925][5] = 31;
    v[926][5] = 45;
    v[927][5] = 51;
    v[928][5] = 49;
    v[929][5] = 55;
    v[930][5] = 17;
    v[931][5] = 43;
    v[932][5] = 49;
    v[933][5] = 45;
    v[934][5] = 9;
    v[935][5] = 29;
    v[936][5] = 3;
    v[937][5] = 5;
    v[938][5] = 47;
    v[939][5] = 9;
    v[940][5] = 15;
    v[941][5] = 19;
    v[942][5] = 51;
    v[943][5] = 45;
    v[944][5] = 57;
    v[945][5] = 63;
    v[946][5] = 9;
    v[947][5] = 21;
    v[948][5] = 59;
    v[949][5] = 3;
    v[950][5] = 9;
    v[951][5] = 13;
    v[952][5] = 45;
    v[953][5] = 23;
    v[954][5] = 15;
    v[955][5] = 31;
    v[956][5] = 21;
    v[957][5] = 15;
    v[958][5] = 51;
    v[959][5] = 35;
    v[960][5] = 9;
    v[961][5] = 11;
    v[962][5] = 61;
    v[963][5] = 23;
    v[964][5] = 53;
    v[965][5] = 29;
    v[966][5] = 51;
    v[967][5] = 45;
    v[968][5] = 31;
    v[969][5] = 29;
    v[970][5] = 5;
    v[971][5] = 35;
    v[972][5] = 29;
    v[973][5] = 53;
    v[974][5] = 35;
    v[975][5] = 17;
    v[976][5] = 59;
    v[977][5] = 55;
    v[978][5] = 27;
    v[979][5] = 51;
    v[980][5] = 59;
    v[981][5] = 27;
    v[982][5] = 47;
    v[983][5] = 15;
    v[984][5] = 29;
    v[985][5] = 37;
    v[986][5] = 7;
    v[987][5] = 49;
    v[988][5] = 55;
    v[989][5] = 5;
    v[990][5] = 19;
    v[991][5] = 45;
    v[992][5] = 29;
    v[993][5] = 19;
    v[994][5] = 57;
    v[995][5] = 33;
    v[996][5] = 53;
    v[997][5] = 45;
    v[998][5] = 21;
    v[999][5] = 9;
    v[1000][5] = 3;
    v[1001][5] = 35;
    v[1002][5] = 29;
    v[1003][5] = 43;
    v[1004][5] = 31;
    v[1005][5] = 39;
    v[1006][5] = 3;
    v[1007][5] = 45;
    v[1008][5] = 1;
    v[1009][5] = 41;
    v[1010][5] = 29;
    v[1011][5] = 5;
    v[1012][5] = 59;
    v[1013][5] = 41;
    v[1014][5] = 33;
    v[1015][5] = 35;
    v[1016][5] = 27;
    v[1017][5] = 19;
    v[1018][5] = 13;
    v[1019][5] = 25;
    v[1020][5] = 27;
    v[1021][5] = 43;
    v[1022][5] = 33;
    v[1023][5] = 35;
    v[1024][5] = 17;
    v[1025][5] = 17;
    v[1026][5] = 23;
    v[1027][5] = 7;
    v[1028][5] = 35;
    v[1029][5] = 15;
    v[1030][5] = 61;
    v[1031][5] = 61;
    v[1032][5] = 53;
    v[1033][5] = 5;
    v[1034][5] = 15;
    v[1035][5] = 23;
    v[1036][5] = 11;
    v[1037][5] = 13;
    v[1038][5] = 43;
    v[1039][5] = 55;
    v[1040][5] = 47;
    v[1041][5] = 25;
    v[1042][5] = 43;
    v[1043][5] = 15;
    v[1044][5] = 57;
    v[1045][5] = 45;
    v[1046][5] = 1;
    v[1047][5] = 49;
    v[1048][5] = 63;
    v[1049][5] = 57;
    v[1050][5] = 15;
    v[1051][5] = 31;
    v[1052][5] = 31;
    v[1053][5] = 7;
    v[1054][5] = 53;
    v[1055][5] = 27;
    v[1056][5] = 15;
    v[1057][5] = 47;
    v[1058][5] = 23;
    v[1059][5] = 7;
    v[1060][5] = 29;
    v[1061][5] = 53;
    v[1062][5] = 47;
    v[1063][5] = 9;
    v[1064][5] = 53;
    v[1065][5] = 3;
    v[1066][5] = 25;
    v[1067][5] = 55;
    v[1068][5] = 45;
    v[1069][5] = 63;
    v[1070][5] = 21;
    v[1071][5] = 17;
    v[1072][5] = 23;
    v[1073][5] = 31;
    v[1074][5] = 27;
    v[1075][5] = 27;
    v[1076][5] = 43;
    v[1077][5] = 63;
    v[1078][5] = 55;
    v[1079][5] = 63;
    v[1080][5] = 45;
    v[1081][5] = 51;
    v[1082][5] = 15;
    v[1083][5] = 27;
    v[1084][5] = 5;
    v[1085][5] = 37;
    v[1086][5] = 43;
    v[1087][5] = 11;
    v[1088][5] = 27;
    v[1089][5] = 5;
    v[1090][5] = 27;
    v[1091][5] = 59;
    v[1092][5] = 21;
    v[1093][5] = 7;
    v[1094][5] = 39;
    v[1095][5] = 27;
    v[1096][5] = 63;
    v[1097][5] = 35;
    v[1098][5] = 47;
    v[1099][5] = 55;
    v[1100][5] = 17;
    v[1101][5] = 17;
    v[1102][5] = 17;
    v[1103][5] = 3;
    v[1104][5] = 19;
    v[1105][5] = 21;
    v[1106][5] = 13;
    v[1107][5] = 49;
    v[1108][5] = 61;
    v[1109][5] = 39;
    v[1110][5] = 15;

    v[19][6] = 13;
    v[20][6] = 33;
    v[21][6] = 115;
    v[22][6] = 41;
    v[23][6] = 79;
    v[24][6] = 17;
    v[25][6] = 29;
    v[26][6] = 119;
    v[27][6] = 75;
    v[28][6] = 73;
    v[29][6] = 105;
    v[30][6] = 7;
    v[31][6] = 59;
    v[32][6] = 65;
    v[33][6] = 21;
    v[34][6] = 3;
    v[35][6] = 113;
    v[36][6] = 61;
    v[37][6] = 89;
    v[38][6] = 45;
    v[39][6] = 107;
    v[40][6] = 21;
    v[41][6] = 71;
    v[42][6] = 79;
    v[43][6] = 19;
    v[44][6] = 71;
    v[45][6] = 61;
    v[46][6] = 41;
    v[47][6] = 57;
    v[48][6] = 121;
    v[49][6] = 87;
    v[50][6] = 119;
    v[51][6] = 55;
    v[52][6] = 85;
    v[53][6] = 121;
    v[54][6] = 119;
    v[55][6] = 11;
    v[56][6] = 23;
    v[57][6] = 61;
    v[58][6] = 11;
    v[59][6] = 35;
    v[60][6] = 33;
    v[61][6] = 43;
    v[62][6] = 107;
    v[63][6] = 113;
    v[64][6] = 101;
    v[65][6] = 29;
    v[66][6] = 87;
    v[67][6] = 119;
    v[68][6] = 97;
    v[69][6] = 29;
    v[70][6] = 17;
    v[71][6] = 89;
    v[72][6] = 5;
    v[73][6] = 127;
    v[74][6] = 89;
    v[75][6] = 119;
    v[76][6] = 117;
    v[77][6] = 103;
    v[78][6] = 105;
    v[79][6] = 41;
    v[80][6] = 83;
    v[81][6] = 25;
    v[82][6] = 41;
    v[83][6] = 55;
    v[84][6] = 69;
    v[85][6] = 117;
    v[86][6] = 49;
    v[87][6] = 127;
    v[88][6] = 29;
    v[89][6] = 1;
    v[90][6] = 99;
    v[91][6] = 53;
    v[92][6] = 83;
    v[93][6] = 15;
    v[94][6] = 31;
    v[95][6] = 73;
    v[96][6] = 115;
    v[97][6] = 35;
    v[98][6] = 21;
    v[99][6] = 89;
    v[100][6] = 5;
    v[101][6] = 1;
    v[102][6] = 91;
    v[103][6] = 53;
    v[104][6] = 35;
    v[105][6] = 95;
    v[106][6] = 83;
    v[107][6] = 19;
    v[108][6] = 85;
    v[109][6] = 55;
    v[110][6] = 51;
    v[111][6] = 101;
    v[112][6] = 33;
    v[113][6] = 41;
    v[114][6] = 55;
    v[115][6] = 45;
    v[116][6] = 95;
    v[117][6] = 61;
    v[118][6] = 27;
    v[119][6] = 37;
    v[120][6] = 89;
    v[121][6] = 75;
    v[122][6] = 57;
    v[123][6] = 61;
    v[124][6] = 15;
    v[125][6] = 117;
    v[126][6] = 15;
    v[127][6] = 21;
    v[128][6] = 27;
    v[129][6] = 25;
    v[130][6] = 27;
    v[131][6] = 123;
    v[132][6] = 39;
    v[133][6] = 109;
    v[134][6] = 93;
    v[135][6] = 51;
    v[136][6] = 21;
    v[137][6] = 91;
    v[138][6] = 109;
    v[139][6] = 107;
    v[140][6] = 45;
    v[141][6] = 15;
    v[142][6] = 93;
    v[143][6] = 127;
    v[144][6] = 3;
    v[145][6] = 53;
    v[146][6] = 81;
    v[147][6] = 79;
    v[148][6] = 107;
    v[149][6] = 79;
    v[150][6] = 87;
    v[151][6] = 35;
    v[152][6] = 109;
    v[153][6] = 73;
    v[154][6] = 35;
    v[155][6] = 83;
    v[156][6] = 107;
    v[157][6] = 1;
    v[158][6] = 51;
    v[159][6] = 7;
    v[160][6] = 59;
    v[161][6] = 33;
    v[162][6] = 115;
    v[163][6] = 43;
    v[164][6] = 111;
    v[165][6] = 45;
    v[166][6] = 121;
    v[167][6] = 105;
    v[168][6] = 125;
    v[169][6] = 87;
    v[170][6] = 101;
    v[171][6] = 41;
    v[172][6] = 95;
    v[173][6] = 75;
    v[174][6] = 1;
    v[175][6] = 57;
    v[176][6] = 117;
    v[177][6] = 21;
    v[178][6] = 27;
    v[179][6] = 67;
    v[180][6] = 29;
    v[181][6] = 53;
    v[182][6] = 117;
    v[183][6] = 63;
    v[184][6] = 1;
    v[185][6] = 77;
    v[186][6] = 89;
    v[187][6] = 115;
    v[188][6] = 49;
    v[189][6] = 127;
    v[190][6] = 15;
    v[191][6] = 79;
    v[192][6] = 81;
    v[193][6] = 29;
    v[194][6] = 65;
    v[195][6] = 103;
    v[196][6] = 33;
    v[197][6] = 73;
    v[198][6] = 79;
    v[199][6] = 29;
    v[200][6] = 21;
    v[201][6] = 113;
    v[202][6] = 31;
    v[203][6] = 33;
    v[204][6] = 107;
    v[205][6] = 95;
    v[206][6] = 111;
    v[207][6] = 59;
    v[208][6] = 99;
    v[209][6] = 117;
    v[210][6] = 63;
    v[211][6] = 63;
    v[212][6] = 99;
    v[213][6] = 39;
    v[214][6] = 9;
    v[215][6] = 35;
    v[216][6] = 63;
    v[217][6] = 125;
    v[218][6] = 99;
    v[219][6] = 45;
    v[220][6] = 93;
    v[221][6] = 33;
    v[222][6] = 93;
    v[223][6] = 9;
    v[224][6] = 105;
    v[225][6] = 75;
    v[226][6] = 51;
    v[227][6] = 115;
    v[228][6] = 11;
    v[229][6] = 37;
    v[230][6] = 17;
    v[231][6] = 41;
    v[232][6] = 21;
    v[233][6] = 43;
    v[234][6] = 73;
    v[235][6] = 19;
    v[236][6] = 93;
    v[237][6] = 7;
    v[238][6] = 95;
    v[239][6] = 81;
    v[240][6] = 93;
    v[241][6] = 79;
    v[242][6] = 81;
    v[243][6] = 55;
    v[244][6] = 9;
    v[245][6] = 51;
    v[246][6] = 63;
    v[247][6] = 45;
    v[248][6] = 89;
    v[249][6] = 73;
    v[250][6] = 19;
    v[251][6] = 115;
    v[252][6] = 39;
    v[253][6] = 47;
    v[254][6] = 81;
    v[255][6] = 39;
    v[256][6] = 5;
    v[257][6] = 5;
    v[258][6] = 45;
    v[259][6] = 53;
    v[260][6] = 65;
    v[261][6] = 49;
    v[262][6] = 17;
    v[263][6] = 105;
    v[264][6] = 13;
    v[265][6] = 107;
    v[266][6] = 5;
    v[267][6] = 5;
    v[268][6] = 19;
    v[269][6] = 73;
    v[270][6] = 59;
    v[271][6] = 43;
    v[272][6] = 83;
    v[273][6] = 97;
    v[274][6] = 115;
    v[275][6] = 27;
    v[276][6] = 1;
    v[277][6] = 69;
    v[278][6] = 103;
    v[279][6] = 3;
    v[280][6] = 99;
    v[281][6] = 103;
    v[282][6] = 63;
    v[283][6] = 67;
    v[284][6] = 25;
    v[285][6] = 121;
    v[286][6] = 97;
    v[287][6] = 77;
    v[288][6] = 13;
    v[289][6] = 83;
    v[290][6] = 103;
    v[291][6] = 41;
    v[292][6] = 11;
    v[293][6] = 27;
    v[294][6] = 81;
    v[295][6] = 37;
    v[296][6] = 33;
    v[297][6] = 125;
    v[298][6] = 71;
    v[299][6] = 41;
    v[300][6] = 41;
    v[301][6] = 59;
    v[302][6] = 41;
    v[303][6] = 87;
    v[304][6] = 123;
    v[305][6] = 43;
    v[306][6] = 101;
    v[307][6] = 63;
    v[308][6] = 45;
    v[309][6] = 39;
    v[310][6] = 21;
    v[311][6] = 97;
    v[312][6] = 15;
    v[313][6] = 97;
    v[314][6] = 111;
    v[315][6] = 21;
    v[316][6] = 49;
    v[317][6] = 13;
    v[318][6] = 17;
    v[319][6] = 79;
    v[320][6] = 91;
    v[321][6] = 65;
    v[322][6] = 105;
    v[323][6] = 75;
    v[324][6] = 1;
    v[325][6] = 45;
    v[326][6] = 67;
    v[327][6] = 83;
    v[328][6] = 107;
    v[329][6] = 125;
    v[330][6] = 87;
    v[331][6] = 15;
    v[332][6] = 81;
    v[333][6] = 95;
    v[334][6] = 105;
    v[335][6] = 65;
    v[336][6] = 45;
    v[337][6] = 59;
    v[338][6] = 103;
    v[339][6] = 23;
    v[340][6] = 103;
    v[341][6] = 99;
    v[342][6] = 67;
    v[343][6] = 99;
    v[344][6] = 47;
    v[345][6] = 117;
    v[346][6] = 71;
    v[347][6] = 89;
    v[348][6] = 35;
    v[349][6] = 53;
    v[350][6] = 73;
    v[351][6] = 9;
    v[352][6] = 115;
    v[353][6] = 49;
    v[354][6] = 37;
    v[355][6] = 1;
    v[356][6] = 35;
    v[357][6] = 9;
    v[358][6] = 45;
    v[359][6] = 81;
    v[360][6] = 19;
    v[361][6] = 127;
    v[362][6] = 17;
    v[363][6] = 17;
    v[364][6] = 105;
    v[365][6] = 89;
    v[366][6] = 49;
    v[367][6] = 101;
    v[368][6] = 7;
    v[369][6] = 37;
    v[370][6] = 33;
    v[371][6] = 11;
    v[372][6] = 95;
    v[373][6] = 95;
    v[374][6] = 17;
    v[375][6] = 111;
    v[376][6] = 105;
    v[377][6] = 41;
    v[378][6] = 115;
    v[379][6] = 5;
    v[380][6] = 69;
    v[381][6] = 101;
    v[382][6] = 27;
    v[383][6] = 27;
    v[384][6] = 101;
    v[385][6] = 103;
    v[386][6] = 53;
    v[387][6] = 9;
    v[388][6] = 21;
    v[389][6] = 43;
    v[390][6] = 79;
    v[391][6] = 91;
    v[392][6] = 65;
    v[393][6] = 117;
    v[394][6] = 87;
    v[395][6] = 125;
    v[396][6] = 55;
    v[397][6] = 45;
    v[398][6] = 63;
    v[399][6] = 85;
    v[400][6] = 83;
    v[401][6] = 97;
    v[402][6] = 45;
    v[403][6] = 83;
    v[404][6] = 87;
    v[405][6] = 113;
    v[406][6] = 93;
    v[407][6] = 95;
    v[408][6] = 5;
    v[409][6] = 17;
    v[410][6] = 77;
    v[411][6] = 77;
    v[412][6] = 127;
    v[413][6] = 123;
    v[414][6] = 45;
    v[415][6] = 81;
    v[416][6] = 85;
    v[417][6] = 121;
    v[418][6] = 119;
    v[419][6] = 27;
    v[420][6] = 85;
    v[421][6] = 41;
    v[422][6] = 49;
    v[423][6] = 15;
    v[424][6] = 107;
    v[425][6] = 21;
    v[426][6] = 51;
    v[427][6] = 119;
    v[428][6] = 11;
    v[429][6] = 87;
    v[430][6] = 101;
    v[431][6] = 115;
    v[432][6] = 63;
    v[433][6] = 63;
    v[434][6] = 37;
    v[435][6] = 121;
    v[436][6] = 109;
    v[437][6] = 7;
    v[438][6] = 43;
    v[439][6] = 69;
    v[440][6] = 19;
    v[441][6] = 77;
    v[442][6] = 49;
    v[443][6] = 71;
    v[444][6] = 59;
    v[445][6] = 35;
    v[446][6] = 7;
    v[447][6] = 13;
    v[448][6] = 55;
    v[449][6] = 101;
    v[450][6] = 127;
    v[451][6] = 103;
    v[452][6] = 85;
    v[453][6] = 109;
    v[454][6] = 29;
    v[455][6] = 61;
    v[456][6] = 67;
    v[457][6] = 21;
    v[458][6] = 111;
    v[459][6] = 67;
    v[460][6] = 23;
    v[461][6] = 57;
    v[462][6] = 75;
    v[463][6] = 71;
    v[464][6] = 101;
    v[465][6] = 123;
    v[466][6] = 41;
    v[467][6] = 107;
    v[468][6] = 101;
    v[469][6] = 107;
    v[470][6] = 125;
    v[471][6] = 27;
    v[472][6] = 47;
    v[473][6] = 119;
    v[474][6] = 41;
    v[475][6] = 19;
    v[476][6] = 127;
    v[477][6] = 33;
    v[478][6] = 31;
    v[479][6] = 109;
    v[480][6] = 7;
    v[481][6] = 91;
    v[482][6] = 91;
    v[483][6] = 39;
    v[484][6] = 125;
    v[485][6] = 105;
    v[486][6] = 47;
    v[487][6] = 125;
    v[488][6] = 123;
    v[489][6] = 91;
    v[490][6] = 9;
    v[491][6] = 103;
    v[492][6] = 45;
    v[493][6] = 23;
    v[494][6] = 117;
    v[495][6] = 9;
    v[496][6] = 125;
    v[497][6] = 73;
    v[498][6] = 11;
    v[499][6] = 37;
    v[500][6] = 61;
    v[501][6] = 79;
    v[502][6] = 21;
    v[503][6] = 5;
    v[504][6] = 47;
    v[505][6] = 117;
    v[506][6] = 67;
    v[507][6] = 53;
    v[508][6] = 85;
    v[509][6] = 33;
    v[510][6] = 81;
    v[511][6] = 121;
    v[512][6] = 47;
    v[513][6] = 61;
    v[514][6] = 51;
    v[515][6] = 127;
    v[516][6] = 29;
    v[517][6] = 65;
    v[518][6] = 45;
    v[519][6] = 41;
    v[520][6] = 95;
    v[521][6] = 57;
    v[522][6] = 73;
    v[523][6] = 33;
    v[524][6] = 117;
    v[525][6] = 61;
    v[526][6] = 111;
    v[527][6] = 59;
    v[528][6] = 123;
    v[529][6] = 65;
    v[530][6] = 47;
    v[531][6] = 105;
    v[532][6] = 23;
    v[533][6] = 29;
    v[534][6] = 107;
    v[535][6] = 37;
    v[536][6] = 81;
    v[537][6] = 67;
    v[538][6] = 29;
    v[539][6] = 115;
    v[540][6] = 119;
    v[541][6] = 75;
    v[542][6] = 73;
    v[543][6] = 99;
    v[544][6] = 103;
    v[545][6] = 7;
    v[546][6] = 57;
    v[547][6] = 45;
    v[548][6] = 61;
    v[549][6] = 95;
    v[550][6] = 49;
    v[551][6] = 101;
    v[552][6] = 101;
    v[553][6] = 35;
    v[554][6] = 47;
    v[555][6] = 119;
    v[556][6] = 39;
    v[557][6] = 67;
    v[558][6] = 31;
    v[559][6] = 103;
    v[560][6] = 7;
    v[561][6] = 61;
    v[562][6] = 127;
    v[563][6] = 87;
    v[564][6] = 3;
    v[565][6] = 35;
    v[566][6] = 29;
    v[567][6] = 73;
    v[568][6] = 95;
    v[569][6] = 103;
    v[570][6] = 71;
    v[571][6] = 75;
    v[572][6] = 51;
    v[573][6] = 87;
    v[574][6] = 57;
    v[575][6] = 97;
    v[576][6] = 11;
    v[577][6] = 105;
    v[578][6] = 87;
    v[579][6] = 41;
    v[580][6] = 73;
    v[581][6] = 109;
    v[582][6] = 69;
    v[583][6] = 35;
    v[584][6] = 121;
    v[585][6] = 39;
    v[586][6] = 111;
    v[587][6] = 1;
    v[588][6] = 77;
    v[589][6] = 39;
    v[590][6] = 47;
    v[591][6] = 53;
    v[592][6] = 91;
    v[593][6] = 3;
    v[594][6] = 17;
    v[595][6] = 51;
    v[596][6] = 83;
    v[597][6] = 39;
    v[598][6] = 125;
    v[599][6] = 85;
    v[600][6] = 111;
    v[601][6] = 21;
    v[602][6] = 69;
    v[603][6] = 85;
    v[604][6] = 29;
    v[605][6] = 55;
    v[606][6] = 11;
    v[607][6] = 117;
    v[608][6] = 1;
    v[609][6] = 47;
    v[610][6] = 17;
    v[611][6] = 65;
    v[612][6] = 63;
    v[613][6] = 47;
    v[614][6] = 117;
    v[615][6] = 17;
    v[616][6] = 115;
    v[617][6] = 51;
    v[618][6] = 25;
    v[619][6] = 33;
    v[620][6] = 123;
    v[621][6] = 123;
    v[622][6] = 83;
    v[623][6] = 51;
    v[624][6] = 113;
    v[625][6] = 95;
    v[626][6] = 121;
    v[627][6] = 51;
    v[628][6] = 91;
    v[629][6] = 109;
    v[630][6] = 43;
    v[631][6] = 55;
    v[632][6] = 35;
    v[633][6] = 55;
    v[634][6] = 87;
    v[635][6] = 33;
    v[636][6] = 37;
    v[637][6] = 5;
    v[638][6] = 3;
    v[639][6] = 45;
    v[640][6] = 21;
    v[641][6] = 105;
    v[642][6] = 127;
    v[643][6] = 35;
    v[644][6] = 17;
    v[645][6] = 35;
    v[646][6] = 37;
    v[647][6] = 97;
    v[648][6] = 97;
    v[649][6] = 21;
    v[650][6] = 77;
    v[651][6] = 123;
    v[652][6] = 17;
    v[653][6] = 89;
    v[654][6] = 53;
    v[655][6] = 105;
    v[656][6] = 75;
    v[657][6] = 25;
    v[658][6] = 125;
    v[659][6] = 13;
    v[660][6] = 47;
    v[661][6] = 21;
    v[662][6] = 125;
    v[663][6] = 23;
    v[664][6] = 55;
    v[665][6] = 63;
    v[666][6] = 61;
    v[667][6] = 5;
    v[668][6] = 17;
    v[669][6] = 93;
    v[670][6] = 57;
    v[671][6] = 121;
    v[672][6] = 69;
    v[673][6] = 73;
    v[674][6] = 93;
    v[675][6] = 121;
    v[676][6] = 105;
    v[677][6] = 75;
    v[678][6] = 91;
    v[679][6] = 67;
    v[680][6] = 95;
    v[681][6] = 75;
    v[682][6] = 9;
    v[683][6] = 69;
    v[684][6] = 97;
    v[685][6] = 99;
    v[686][6] = 93;
    v[687][6] = 11;
    v[688][6] = 53;
    v[689][6] = 19;
    v[690][6] = 73;
    v[691][6] = 5;
    v[692][6] = 33;
    v[693][6] = 79;
    v[694][6] = 107;
    v[695][6] = 65;
    v[696][6] = 69;
    v[697][6] = 79;
    v[698][6] = 125;
    v[699][6] = 25;
    v[700][6] = 93;
    v[701][6] = 55;
    v[702][6] = 61;
    v[703][6] = 17;
    v[704][6] = 117;
    v[705][6] = 69;
    v[706][6] = 97;
    v[707][6] = 87;
    v[708][6] = 111;
    v[709][6] = 37;
    v[710][6] = 93;
    v[711][6] = 59;
    v[712][6] = 79;
    v[713][6] = 95;
    v[714][6] = 53;
    v[715][6] = 115;
    v[716][6] = 53;
    v[717][6] = 85;
    v[718][6] = 85;
    v[719][6] = 65;
    v[720][6] = 59;
    v[721][6] = 23;
    v[722][6] = 75;
    v[723][6] = 21;
    v[724][6] = 67;
    v[725][6] = 27;
    v[726][6] = 99;
    v[727][6] = 79;
    v[728][6] = 27;
    v[729][6] = 3;
    v[730][6] = 95;
    v[731][6] = 27;
    v[732][6] = 69;
    v[733][6] = 19;
    v[734][6] = 75;
    v[735][6] = 47;
    v[736][6] = 59;
    v[737][6] = 41;
    v[738][6] = 85;
    v[739][6] = 77;
    v[740][6] = 99;
    v[741][6] = 55;
    v[742][6] = 49;
    v[743][6] = 93;
    v[744][6] = 93;
    v[745][6] = 119;
    v[746][6] = 51;
    v[747][6] = 125;
    v[748][6] = 63;
    v[749][6] = 13;
    v[750][6] = 15;
    v[751][6] = 45;
    v[752][6] = 61;
    v[753][6] = 19;
    v[754][6] = 105;
    v[755][6] = 115;
    v[756][6] = 17;
    v[757][6] = 83;
    v[758][6] = 7;
    v[759][6] = 7;
    v[760][6] = 11;
    v[761][6] = 61;
    v[762][6] = 37;
    v[763][6] = 63;
    v[764][6] = 89;
    v[765][6] = 95;
    v[766][6] = 119;
    v[767][6] = 113;
    v[768][6] = 67;
    v[769][6] = 123;
    v[770][6] = 91;
    v[771][6] = 33;
    v[772][6] = 37;
    v[773][6] = 99;
    v[774][6] = 43;
    v[775][6] = 11;
    v[776][6] = 33;
    v[777][6] = 65;
    v[778][6] = 81;
    v[779][6] = 79;
    v[780][6] = 81;
    v[781][6] = 107;
    v[782][6] = 63;
    v[783][6] = 63;
    v[784][6] = 55;
    v[785][6] = 89;
    v[786][6] = 91;
    v[787][6] = 25;
    v[788][6] = 93;
    v[789][6] = 101;
    v[790][6] = 27;
    v[791][6] = 55;
    v[792][6] = 75;
    v[793][6] = 121;
    v[794][6] = 79;
    v[795][6] = 43;
    v[796][6] = 125;
    v[797][6] = 73;
    v[798][6] = 27;
    v[799][6] = 109;
    v[800][6] = 35;
    v[801][6] = 21;
    v[802][6] = 71;
    v[803][6] = 113;
    v[804][6] = 89;
    v[805][6] = 59;
    v[806][6] = 95;
    v[807][6] = 41;
    v[808][6] = 45;
    v[809][6] = 113;
    v[810][6] = 119;
    v[811][6] = 113;
    v[812][6] = 39;
    v[813][6] = 59;
    v[814][6] = 73;
    v[815][6] = 15;
    v[816][6] = 13;
    v[817][6] = 59;
    v[818][6] = 67;
    v[819][6] = 121;
    v[820][6] = 27;
    v[821][6] = 7;
    v[822][6] = 105;
    v[823][6] = 15;
    v[824][6] = 59;
    v[825][6] = 59;
    v[826][6] = 35;
    v[827][6] = 91;
    v[828][6] = 89;
    v[829][6] = 23;
    v[830][6] = 125;
    v[831][6] = 97;
    v[832][6] = 53;
    v[833][6] = 41;
    v[834][6] = 91;
    v[835][6] = 111;
    v[836][6] = 29;
    v[837][6] = 31;
    v[838][6] = 3;
    v[839][6] = 103;
    v[840][6] = 61;
    v[841][6] = 71;
    v[842][6] = 35;
    v[843][6] = 7;
    v[844][6] = 119;
    v[845][6] = 29;
    v[846][6] = 45;
    v[847][6] = 49;
    v[848][6] = 111;
    v[849][6] = 41;
    v[850][6] = 109;
    v[851][6] = 59;
    v[852][6] = 125;
    v[853][6] = 13;
    v[854][6] = 27;
    v[855][6] = 19;
    v[856][6] = 79;
    v[857][6] = 9;
    v[858][6] = 75;
    v[859][6] = 83;
    v[860][6] = 81;
    v[861][6] = 33;
    v[862][6] = 91;
    v[863][6] = 109;
    v[864][6] = 33;
    v[865][6] = 29;
    v[866][6] = 107;
    v[867][6] = 111;
    v[868][6] = 101;
    v[869][6] = 107;
    v[870][6] = 109;
    v[871][6] = 65;
    v[872][6] = 59;
    v[873][6] = 43;
    v[874][6] = 37;
    v[875][6] = 1;
    v[876][6] = 9;
    v[877][6] = 15;
    v[878][6] = 109;
    v[879][6] = 37;
    v[880][6] = 111;
    v[881][6] = 113;
    v[882][6] = 119;
    v[883][6] = 79;
    v[884][6] = 73;
    v[885][6] = 65;
    v[886][6] = 71;
    v[887][6] = 93;
    v[888][6] = 17;
    v[889][6] = 101;
    v[890][6] = 87;
    v[891][6] = 97;
    v[892][6] = 43;
    v[893][6] = 23;
    v[894][6] = 75;
    v[895][6] = 109;
    v[896][6] = 41;
    v[897][6] = 49;
    v[898][6] = 53;
    v[899][6] = 31;
    v[900][6] = 97;
    v[901][6] = 105;
    v[902][6] = 109;
    v[903][6] = 119;
    v[904][6] = 51;
    v[905][6] = 9;
    v[906][6] = 53;
    v[907][6] = 113;
    v[908][6] = 97;
    v[909][6] = 73;
    v[910][6] = 89;
    v[911][6] = 79;
    v[912][6] = 49;
    v[913][6] = 61;
    v[914][6] = 105;
    v[915][6] = 13;
    v[916][6] = 99;
    v[917][6] = 53;
    v[918][6] = 71;
    v[919][6] = 7;
    v[920][6] = 87;
    v[921][6] = 21;
    v[922][6] = 101;
    v[923][6] = 5;
    v[924][6] = 71;
    v[925][6] = 31;
    v[926][6] = 123;
    v[927][6] = 121;
    v[928][6] = 121;
    v[929][6] = 73;
    v[930][6] = 79;
    v[931][6] = 115;
    v[932][6] = 13;
    v[933][6] = 39;
    v[934][6] = 101;
    v[935][6] = 19;
    v[936][6] = 37;
    v[937][6] = 51;
    v[938][6] = 83;
    v[939][6] = 97;
    v[940][6] = 55;
    v[941][6] = 81;
    v[942][6] = 91;
    v[943][6] = 127;
    v[944][6] = 105;
    v[945][6] = 89;
    v[946][6] = 63;
    v[947][6] = 47;
    v[948][6] = 49;
    v[949][6] = 75;
    v[950][6] = 37;
    v[951][6] = 77;
    v[952][6] = 15;
    v[953][6] = 49;
    v[954][6] = 107;
    v[955][6] = 23;
    v[956][6] = 23;
    v[957][6] = 35;
    v[958][6] = 19;
    v[959][6] = 69;
    v[960][6] = 17;
    v[961][6] = 59;
    v[962][6] = 63;
    v[963][6] = 73;
    v[964][6] = 29;
    v[965][6] = 125;
    v[966][6] = 61;
    v[967][6] = 65;
    v[968][6] = 95;
    v[969][6] = 101;
    v[970][6] = 81;
    v[971][6] = 57;
    v[972][6] = 69;
    v[973][6] = 83;
    v[974][6] = 37;
    v[975][6] = 11;
    v[976][6] = 37;
    v[977][6] = 95;
    v[978][6] = 1;
    v[979][6] = 73;
    v[980][6] = 27;
    v[981][6] = 29;
    v[982][6] = 57;
    v[983][6] = 7;
    v[984][6] = 65;
    v[985][6] = 83;
    v[986][6] = 99;
    v[987][6] = 69;
    v[988][6] = 19;
    v[989][6] = 103;
    v[990][6] = 43;
    v[991][6] = 95;
    v[992][6] = 25;
    v[993][6] = 19;
    v[994][6] = 103;
    v[995][6] = 41;
    v[996][6] = 125;
    v[997][6] = 97;
    v[998][6] = 71;
    v[999][6] = 105;
    v[1000][6] = 83;
    v[1001][6] = 83;
    v[1002][6] = 61;
    v[1003][6] = 39;
    v[1004][6] = 9;
    v[1005][6] = 45;
    v[1006][6] = 117;
    v[1007][6] = 63;
    v[1008][6] = 31;
    v[1009][6] = 5;
    v[1010][6] = 117;
    v[1011][6] = 67;
    v[1012][6] = 125;
    v[1013][6] = 41;
    v[1014][6] = 117;
    v[1015][6] = 43;
    v[1016][6] = 77;
    v[1017][6] = 97;
    v[1018][6] = 15;
    v[1019][6] = 29;
    v[1020][6] = 5;
    v[1021][6] = 59;
    v[1022][6] = 25;
    v[1023][6] = 63;
    v[1024][6] = 87;
    v[1025][6] = 39;
    v[1026][6] = 39;
    v[1027][6] = 77;
    v[1028][6] = 85;
    v[1029][6] = 37;
    v[1030][6] = 81;
    v[1031][6] = 73;
    v[1032][6] = 89;
    v[1033][6] = 29;
    v[1034][6] = 125;
    v[1035][6] = 109;
    v[1036][6] = 21;
    v[1037][6] = 23;
    v[1038][6] = 119;
    v[1039][6] = 105;
    v[1040][6] = 43;
    v[1041][6] = 93;
    v[1042][6] = 97;
    v[1043][6] = 15;
    v[1044][6] = 125;
    v[1045][6] = 29;
    v[1046][6] = 51;
    v[1047][6] = 69;
    v[1048][6] = 37;
    v[1049][6] = 45;
    v[1050][6] = 31;
    v[1051][6] = 75;
    v[1052][6] = 109;
    v[1053][6] = 119;
    v[1054][6] = 53;
    v[1055][6] = 5;
    v[1056][6] = 101;
    v[1057][6] = 125;
    v[1058][6] = 121;
    v[1059][6] = 35;
    v[1060][6] = 29;
    v[1061][6] = 7;
    v[1062][6] = 63;
    v[1063][6] = 17;
    v[1064][6] = 63;
    v[1065][6] = 13;
    v[1066][6] = 69;
    v[1067][6] = 15;
    v[1068][6] = 105;
    v[1069][6] = 51;
    v[1070][6] = 127;
    v[1071][6] = 105;
    v[1072][6] = 9;
    v[1073][6] = 57;
    v[1074][6] = 95;
    v[1075][6] = 59;
    v[1076][6] = 109;
    v[1077][6] = 35;
    v[1078][6] = 49;
    v[1079][6] = 23;
    v[1080][6] = 33;
    v[1081][6] = 107;
    v[1082][6] = 55;
    v[1083][6] = 33;
    v[1084][6] = 57;
    v[1085][6] = 79;
    v[1086][6] = 73;
    v[1087][6] = 69;
    v[1088][6] = 59;
    v[1089][6] = 107;
    v[1090][6] = 55;
    v[1091][6] = 11;
    v[1092][6] = 63;
    v[1093][6] = 95;
    v[1094][6] = 103;
    v[1095][6] = 23;
    v[1096][6] = 125;
    v[1097][6] = 91;
    v[1098][6] = 31;
    v[1099][6] = 91;
    v[1100][6] = 51;
    v[1101][6] = 65;
    v[1102][6] = 61;
    v[1103][6] = 75;
    v[1104][6] = 69;
    v[1105][6] = 107;
    v[1106][6] = 65;
    v[1107][6] = 101;
    v[1108][6] = 59;
    v[1109][6] = 35;
    v[1110][6] = 15;

    v[37][7] = 7;
    v[38][7] = 23;
    v[39][7] = 39;
    v[40][7] = 217;
    v[41][7] = 141;
    v[42][7] = 27;
    v[43][7] = 53;
    v[44][7] = 181;
    v[45][7] = 169;
    v[46][7] = 35;
    v[47][7] = 15;
    v[48][7] = 207;
    v[49][7] = 45;
    v[50][7] = 247;
    v[51][7] = 185;
    v[52][7] = 117;
    v[53][7] = 41;
    v[54][7] = 81;
    v[55][7] = 223;
    v[56][7] = 151;
    v[57][7] = 81;
    v[58][7] = 189;
    v[59][7] = 61;
    v[60][7] = 95;
    v[61][7] = 185;
    v[62][7] = 23;
    v[63][7] = 73;
    v[64][7] = 113;
    v[65][7] = 239;
    v[66][7] = 85;
    v[67][7] = 9;
    v[68][7] = 201;
    v[69][7] = 83;
    v[70][7] = 53;
    v[71][7] = 183;
    v[72][7] = 203;
    v[73][7] = 91;
    v[74][7] = 149;
    v[75][7] = 101;
    v[76][7] = 13;
    v[77][7] = 111;
    v[78][7] = 239;
    v[79][7] = 3;
    v[80][7] = 205;
    v[81][7] = 253;
    v[82][7] = 247;
    v[83][7] = 121;
    v[84][7] = 189;
    v[85][7] = 169;
    v[86][7] = 179;
    v[87][7] = 197;
    v[88][7] = 175;
    v[89][7] = 217;
    v[90][7] = 249;
    v[91][7] = 195;
    v[92][7] = 95;
    v[93][7] = 63;
    v[94][7] = 19;
    v[95][7] = 7;
    v[96][7] = 5;
    v[97][7] = 75;
    v[98][7] = 217;
    v[99][7] = 245;
    v[100][7] = 111;
    v[101][7] = 189;
    v[102][7] = 165;
    v[103][7] = 169;
    v[104][7] = 141;
    v[105][7] = 221;
    v[106][7] = 249;
    v[107][7] = 159;
    v[108][7] = 253;
    v[109][7] = 207;
    v[110][7] = 249;
    v[111][7] = 219;
    v[112][7] = 23;
    v[113][7] = 49;
    v[114][7] = 127;
    v[115][7] = 237;
    v[116][7] = 5;
    v[117][7] = 25;
    v[118][7] = 177;
    v[119][7] = 37;
    v[120][7] = 103;
    v[121][7] = 65;
    v[122][7] = 167;
    v[123][7] = 81;
    v[124][7] = 87;
    v[125][7] = 119;
    v[126][7] = 45;
    v[127][7] = 79;
    v[128][7] = 143;
    v[129][7] = 57;
    v[130][7] = 79;
    v[131][7] = 187;
    v[132][7] = 143;
    v[133][7] = 183;
    v[134][7] = 75;
    v[135][7] = 97;
    v[136][7] = 211;
    v[137][7] = 149;
    v[138][7] = 175;
    v[139][7] = 37;
    v[140][7] = 135;
    v[141][7] = 189;
    v[142][7] = 225;
    v[143][7] = 241;
    v[144][7] = 63;
    v[145][7] = 33;
    v[146][7] = 43;
    v[147][7] = 13;
    v[148][7] = 73;
    v[149][7] = 213;
    v[150][7] = 57;
    v[151][7] = 239;
    v[152][7] = 183;
    v[153][7] = 117;
    v[154][7] = 21;
    v[155][7] = 29;
    v[156][7] = 115;
    v[157][7] = 43;
    v[158][7] = 205;
    v[159][7] = 223;
    v[160][7] = 15;
    v[161][7] = 3;
    v[162][7] = 159;
    v[163][7] = 51;
    v[164][7] = 101;
    v[165][7] = 127;
    v[166][7] = 99;
    v[167][7] = 239;
    v[168][7] = 171;
    v[169][7] = 113;
    v[170][7] = 171;
    v[171][7] = 119;
    v[172][7] = 189;
    v[173][7] = 245;
    v[174][7] = 201;
    v[175][7] = 27;
    v[176][7] = 185;
    v[177][7] = 229;
    v[178][7] = 105;
    v[179][7] = 153;
    v[180][7] = 189;
    v[181][7] = 33;
    v[182][7] = 35;
    v[183][7] = 137;
    v[184][7] = 77;
    v[185][7] = 97;
    v[186][7] = 17;
    v[187][7] = 181;
    v[188][7] = 55;
    v[189][7] = 197;
    v[190][7] = 201;
    v[191][7] = 155;
    v[192][7] = 37;
    v[193][7] = 197;
    v[194][7] = 137;
    v[195][7] = 223;
    v[196][7] = 25;
    v[197][7] = 179;
    v[198][7] = 91;
    v[199][7] = 23;
    v[200][7] = 235;
    v[201][7] = 53;
    v[202][7] = 253;
    v[203][7] = 49;
    v[204][7] = 181;
    v[205][7] = 249;
    v[206][7] = 53;
    v[207][7] = 173;
    v[208][7] = 97;
    v[209][7] = 247;
    v[210][7] = 67;
    v[211][7] = 115;
    v[212][7] = 103;
    v[213][7] = 159;
    v[214][7] = 239;
    v[215][7] = 69;
    v[216][7] = 173;
    v[217][7] = 217;
    v[218][7] = 95;
    v[219][7] = 221;
    v[220][7] = 247;
    v[221][7] = 97;
    v[222][7] = 91;
    v[223][7] = 123;
    v[224][7] = 223;
    v[225][7] = 213;
    v[226][7] = 129;
    v[227][7] = 181;
    v[228][7] = 87;
    v[229][7] = 239;
    v[230][7] = 85;
    v[231][7] = 89;
    v[232][7] = 249;
    v[233][7] = 141;
    v[234][7] = 39;
    v[235][7] = 57;
    v[236][7] = 249;
    v[237][7] = 71;
    v[238][7] = 101;
    v[239][7] = 159;
    v[240][7] = 33;
    v[241][7] = 137;
    v[242][7] = 189;
    v[243][7] = 71;
    v[244][7] = 253;
    v[245][7] = 205;
    v[246][7] = 171;
    v[247][7] = 13;
    v[248][7] = 249;
    v[249][7] = 109;
    v[250][7] = 131;
    v[251][7] = 199;
    v[252][7] = 189;
    v[253][7] = 179;
    v[254][7] = 31;
    v[255][7] = 99;
    v[256][7] = 113;
    v[257][7] = 41;
    v[258][7] = 173;
    v[259][7] = 23;
    v[260][7] = 189;
    v[261][7] = 197;
    v[262][7] = 3;
    v[263][7] = 135;
    v[264][7] = 9;
    v[265][7] = 95;
    v[266][7] = 195;
    v[267][7] = 27;
    v[268][7] = 183;
    v[269][7] = 1;
    v[270][7] = 123;
    v[271][7] = 73;
    v[272][7] = 53;
    v[273][7] = 99;
    v[274][7] = 197;
    v[275][7] = 59;
    v[276][7] = 27;
    v[277][7] = 101;
    v[278][7] = 55;
    v[279][7] = 193;
    v[280][7] = 31;
    v[281][7] = 61;
    v[282][7] = 119;
    v[283][7] = 11;
    v[284][7] = 7;
    v[285][7] = 255;
    v[286][7] = 233;
    v[287][7] = 53;
    v[288][7] = 157;
    v[289][7] = 193;
    v[290][7] = 97;
    v[291][7] = 83;
    v[292][7] = 65;
    v[293][7] = 81;
    v[294][7] = 239;
    v[295][7] = 167;
    v[296][7] = 69;
    v[297][7] = 71;
    v[298][7] = 109;
    v[299][7] = 97;
    v[300][7] = 137;
    v[301][7] = 71;
    v[302][7] = 193;
    v[303][7] = 189;
    v[304][7] = 115;
    v[305][7] = 79;
    v[306][7] = 205;
    v[307][7] = 37;
    v[308][7] = 227;
    v[309][7] = 53;
    v[310][7] = 33;
    v[311][7] = 91;
    v[312][7] = 229;
    v[313][7] = 245;
    v[314][7] = 105;
    v[315][7] = 77;
    v[316][7] = 229;
    v[317][7] = 161;
    v[318][7] = 103;
    v[319][7] = 93;
    v[320][7] = 13;
    v[321][7] = 161;
    v[322][7] = 229;
    v[323][7] = 223;
    v[324][7] = 69;
    v[325][7] = 15;
    v[326][7] = 25;
    v[327][7] = 23;
    v[328][7] = 233;
    v[329][7] = 93;
    v[330][7] = 25;
    v[331][7] = 217;
    v[332][7] = 247;
    v[333][7] = 61;
    v[334][7] = 75;
    v[335][7] = 27;
    v[336][7] = 9;
    v[337][7] = 223;
    v[338][7] = 213;
    v[339][7] = 55;
    v[340][7] = 197;
    v[341][7] = 145;
    v[342][7] = 89;
    v[343][7] = 199;
    v[344][7] = 41;
    v[345][7] = 201;
    v[346][7] = 5;
    v[347][7] = 149;
    v[348][7] = 35;
    v[349][7] = 119;
    v[350][7] = 183;
    v[351][7] = 53;
    v[352][7] = 11;
    v[353][7] = 13;
    v[354][7] = 3;
    v[355][7] = 179;
    v[356][7] = 229;
    v[357][7] = 43;
    v[358][7] = 55;
    v[359][7] = 187;
    v[360][7] = 233;
    v[361][7] = 47;
    v[362][7] = 133;
    v[363][7] = 91;
    v[364][7] = 47;
    v[365][7] = 71;
    v[366][7] = 93;
    v[367][7] = 105;
    v[368][7] = 145;
    v[369][7] = 45;
    v[370][7] = 255;
    v[371][7] = 221;
    v[372][7] = 115;
    v[373][7] = 175;
    v[374][7] = 19;
    v[375][7] = 129;
    v[376][7] = 5;
    v[377][7] = 209;
    v[378][7] = 197;
    v[379][7] = 57;
    v[380][7] = 177;
    v[381][7] = 115;
    v[382][7] = 187;
    v[383][7] = 119;
    v[384][7] = 77;
    v[385][7] = 211;
    v[386][7] = 111;
    v[387][7] = 33;
    v[388][7] = 113;
    v[389][7] = 23;
    v[390][7] = 87;
    v[391][7] = 137;
    v[392][7] = 41;
    v[393][7] = 7;
    v[394][7] = 83;
    v[395][7] = 43;
    v[396][7] = 121;
    v[397][7] = 145;
    v[398][7] = 5;
    v[399][7] = 219;
    v[400][7] = 27;
    v[401][7] = 11;
    v[402][7] = 111;
    v[403][7] = 207;
    v[404][7] = 55;
    v[405][7] = 97;
    v[406][7] = 63;
    v[407][7] = 229;
    v[408][7] = 53;
    v[409][7] = 33;
    v[410][7] = 149;
    v[411][7] = 23;
    v[412][7] = 187;
    v[413][7] = 153;
    v[414][7] = 91;
    v[415][7] = 193;
    v[416][7] = 183;
    v[417][7] = 59;
    v[418][7] = 211;
    v[419][7] = 93;
    v[420][7] = 139;
    v[421][7] = 59;
    v[422][7] = 179;
    v[423][7] = 163;
    v[424][7] = 209;
    v[425][7] = 77;
    v[426][7] = 39;
    v[427][7] = 111;
    v[428][7] = 79;
    v[429][7] = 229;
    v[430][7] = 85;
    v[431][7] = 237;
    v[432][7] = 199;
    v[433][7] = 137;
    v[434][7] = 147;
    v[435][7] = 25;
    v[436][7] = 73;
    v[437][7] = 121;
    v[438][7] = 129;
    v[439][7] = 83;
    v[440][7] = 87;
    v[441][7] = 93;
    v[442][7] = 205;
    v[443][7] = 167;
    v[444][7] = 53;
    v[445][7] = 107;
    v[446][7] = 229;
    v[447][7] = 213;
    v[448][7] = 95;
    v[449][7] = 219;
    v[450][7] = 109;
    v[451][7] = 175;
    v[452][7] = 13;
    v[453][7] = 209;
    v[454][7] = 97;
    v[455][7] = 61;
    v[456][7] = 147;
    v[457][7] = 19;
    v[458][7] = 13;
    v[459][7] = 123;
    v[460][7] = 73;
    v[461][7] = 35;
    v[462][7] = 141;
    v[463][7] = 81;
    v[464][7] = 19;
    v[465][7] = 171;
    v[466][7] = 255;
    v[467][7] = 111;
    v[468][7] = 107;
    v[469][7] = 233;
    v[470][7] = 113;
    v[471][7] = 133;
    v[472][7] = 89;
    v[473][7] = 9;
    v[474][7] = 231;
    v[475][7] = 95;
    v[476][7] = 69;
    v[477][7] = 33;
    v[478][7] = 1;
    v[479][7] = 253;
    v[480][7] = 219;
    v[481][7] = 253;
    v[482][7] = 247;
    v[483][7] = 129;
    v[484][7] = 11;
    v[485][7] = 251;
    v[486][7] = 221;
    v[487][7] = 153;
    v[488][7] = 35;
    v[489][7] = 103;
    v[490][7] = 239;
    v[491][7] = 7;
    v[492][7] = 27;
    v[493][7] = 235;
    v[494][7] = 181;
    v[495][7] = 5;
    v[496][7] = 207;
    v[497][7] = 53;
    v[498][7] = 149;
    v[499][7] = 155;
    v[500][7] = 225;
    v[501][7] = 165;
    v[502][7] = 137;
    v[503][7] = 155;
    v[504][7] = 201;
    v[505][7] = 97;
    v[506][7] = 245;
    v[507][7] = 203;
    v[508][7] = 47;
    v[509][7] = 39;
    v[510][7] = 35;
    v[511][7] = 105;
    v[512][7] = 239;
    v[513][7] = 49;
    v[514][7] = 15;
    v[515][7] = 253;
    v[516][7] = 7;
    v[517][7] = 237;
    v[518][7] = 213;
    v[519][7] = 55;
    v[520][7] = 87;
    v[521][7] = 199;
    v[522][7] = 27;
    v[523][7] = 175;
    v[524][7] = 49;
    v[525][7] = 41;
    v[526][7] = 229;
    v[527][7] = 85;
    v[528][7] = 3;
    v[529][7] = 149;
    v[530][7] = 179;
    v[531][7] = 129;
    v[532][7] = 185;
    v[533][7] = 249;
    v[534][7] = 197;
    v[535][7] = 15;
    v[536][7] = 97;
    v[537][7] = 197;
    v[538][7] = 139;
    v[539][7] = 203;
    v[540][7] = 63;
    v[541][7] = 33;
    v[542][7] = 251;
    v[543][7] = 217;
    v[544][7] = 199;
    v[545][7] = 199;
    v[546][7] = 99;
    v[547][7] = 249;
    v[548][7] = 33;
    v[549][7] = 229;
    v[550][7] = 177;
    v[551][7] = 13;
    v[552][7] = 209;
    v[553][7] = 147;
    v[554][7] = 97;
    v[555][7] = 31;
    v[556][7] = 125;
    v[557][7] = 177;
    v[558][7] = 137;
    v[559][7] = 187;
    v[560][7] = 11;
    v[561][7] = 91;
    v[562][7] = 223;
    v[563][7] = 29;
    v[564][7] = 169;
    v[565][7] = 231;
    v[566][7] = 59;
    v[567][7] = 31;
    v[568][7] = 163;
    v[569][7] = 41;
    v[570][7] = 57;
    v[571][7] = 87;
    v[572][7] = 247;
    v[573][7] = 25;
    v[574][7] = 127;
    v[575][7] = 101;
    v[576][7] = 207;
    v[577][7] = 187;
    v[578][7] = 73;
    v[579][7] = 61;
    v[580][7] = 105;
    v[581][7] = 27;
    v[582][7] = 91;
    v[583][7] = 171;
    v[584][7] = 243;
    v[585][7] = 33;
    v[586][7] = 3;
    v[587][7] = 1;
    v[588][7] = 21;
    v[589][7] = 229;
    v[590][7] = 93;
    v[591][7] = 71;
    v[592][7] = 61;
    v[593][7] = 37;
    v[594][7] = 183;
    v[595][7] = 65;
    v[596][7] = 211;
    v[597][7] = 53;
    v[598][7] = 11;
    v[599][7] = 151;
    v[600][7] = 165;
    v[601][7] = 47;
    v[602][7] = 5;
    v[603][7] = 129;
    v[604][7] = 79;
    v[605][7] = 101;
    v[606][7] = 147;
    v[607][7] = 169;
    v[608][7] = 181;
    v[609][7] = 19;
    v[610][7] = 95;
    v[611][7] = 77;
    v[612][7] = 139;
    v[613][7] = 197;
    v[614][7] = 219;
    v[615][7] = 97;
    v[616][7] = 239;
    v[617][7] = 183;
    v[618][7] = 143;
    v[619][7] = 9;
    v[620][7] = 13;
    v[621][7] = 209;
    v[622][7] = 23;
    v[623][7] = 215;
    v[624][7] = 53;
    v[625][7] = 137;
    v[626][7] = 203;
    v[627][7] = 19;
    v[628][7] = 151;
    v[629][7] = 171;
    v[630][7] = 133;
    v[631][7] = 219;
    v[632][7] = 231;
    v[633][7] = 3;
    v[634][7] = 15;
    v[635][7] = 253;
    v[636][7] = 225;
    v[637][7] = 33;
    v[638][7] = 111;
    v[639][7] = 183;
    v[640][7] = 213;
    v[641][7] = 169;
    v[642][7] = 119;
    v[643][7] = 111;
    v[644][7] = 15;
    v[645][7] = 201;
    v[646][7] = 123;
    v[647][7] = 121;
    v[648][7] = 225;
    v[649][7] = 113;
    v[650][7] = 113;
    v[651][7] = 225;
    v[652][7] = 161;
    v[653][7] = 165;
    v[654][7] = 1;
    v[655][7] = 139;
    v[656][7] = 55;
    v[657][7] = 3;
    v[658][7] = 93;
    v[659][7] = 217;
    v[660][7] = 193;
    v[661][7] = 97;
    v[662][7] = 29;
    v[663][7] = 69;
    v[664][7] = 231;
    v[665][7] = 161;
    v[666][7] = 93;
    v[667][7] = 69;
    v[668][7] = 143;
    v[669][7] = 137;
    v[670][7] = 9;
    v[671][7] = 87;
    v[672][7] = 183;
    v[673][7] = 113;
    v[674][7] = 183;
    v[675][7] = 73;
    v[676][7] = 215;
    v[677][7] = 137;
    v[678][7] = 89;
    v[679][7] = 251;
    v[680][7] = 163;
    v[681][7] = 41;
    v[682][7] = 227;
    v[683][7] = 145;
    v[684][7] = 57;
    v[685][7] = 81;
    v[686][7] = 57;
    v[687][7] = 11;
    v[688][7] = 135;
    v[689][7] = 145;
    v[690][7] = 161;
    v[691][7] = 175;
    v[692][7] = 159;
    v[693][7] = 25;
    v[694][7] = 55;
    v[695][7] = 167;
    v[696][7] = 157;
    v[697][7] = 211;
    v[698][7] = 97;
    v[699][7] = 247;
    v[700][7] = 249;
    v[701][7] = 23;
    v[702][7] = 129;
    v[703][7] = 159;
    v[704][7] = 71;
    v[705][7] = 197;
    v[706][7] = 127;
    v[707][7] = 141;
    v[708][7] = 219;
    v[709][7] = 5;
    v[710][7] = 233;
    v[711][7] = 131;
    v[712][7] = 217;
    v[713][7] = 101;
    v[714][7] = 131;
    v[715][7] = 33;
    v[716][7] = 157;
    v[717][7] = 173;
    v[718][7] = 69;
    v[719][7] = 207;
    v[720][7] = 239;
    v[721][7] = 81;
    v[722][7] = 205;
    v[723][7] = 11;
    v[724][7] = 41;
    v[725][7] = 169;
    v[726][7] = 65;
    v[727][7] = 193;
    v[728][7] = 77;
    v[729][7] = 201;
    v[730][7] = 173;
    v[731][7] = 1;
    v[732][7] = 221;
    v[733][7] = 157;
    v[734][7] = 1;
    v[735][7] = 15;
    v[736][7] = 113;
    v[737][7] = 147;
    v[738][7] = 137;
    v[739][7] = 205;
    v[740][7] = 225;
    v[741][7] = 73;
    v[742][7] = 45;
    v[743][7] = 49;
    v[744][7] = 149;
    v[745][7] = 113;
    v[746][7] = 253;
    v[747][7] = 99;
    v[748][7] = 17;
    v[749][7] = 119;
    v[750][7] = 105;
    v[751][7] = 117;
    v[752][7] = 129;
    v[753][7] = 243;
    v[754][7] = 75;
    v[755][7] = 203;
    v[756][7] = 53;
    v[757][7] = 29;
    v[758][7] = 247;
    v[759][7] = 35;
    v[760][7] = 247;
    v[761][7] = 171;
    v[762][7] = 31;
    v[763][7] = 199;
    v[764][7] = 213;
    v[765][7] = 29;
    v[766][7] = 251;
    v[767][7] = 7;
    v[768][7] = 251;
    v[769][7] = 187;
    v[770][7] = 91;
    v[771][7] = 11;
    v[772][7] = 149;
    v[773][7] = 13;
    v[774][7] = 205;
    v[775][7] = 37;
    v[776][7] = 249;
    v[777][7] = 137;
    v[778][7] = 139;
    v[779][7] = 9;
    v[780][7] = 7;
    v[781][7] = 113;
    v[782][7] = 183;
    v[783][7] = 205;
    v[784][7] = 187;
    v[785][7] = 39;
    v[786][7] = 3;
    v[787][7] = 79;
    v[788][7] = 155;
    v[789][7] = 227;
    v[790][7] = 89;
    v[791][7] = 185;
    v[792][7] = 51;
    v[793][7] = 127;
    v[794][7] = 63;
    v[795][7] = 83;
    v[796][7] = 41;
    v[797][7] = 133;
    v[798][7] = 183;
    v[799][7] = 181;
    v[800][7] = 127;
    v[801][7] = 19;
    v[802][7] = 255;
    v[803][7] = 219;
    v[804][7] = 59;
    v[805][7] = 251;
    v[806][7] = 3;
    v[807][7] = 187;
    v[808][7] = 57;
    v[809][7] = 217;
    v[810][7] = 115;
    v[811][7] = 217;
    v[812][7] = 229;
    v[813][7] = 181;
    v[814][7] = 185;
    v[815][7] = 149;
    v[816][7] = 83;
    v[817][7] = 115;
    v[818][7] = 11;
    v[819][7] = 123;
    v[820][7] = 19;
    v[821][7] = 109;
    v[822][7] = 165;
    v[823][7] = 103;
    v[824][7] = 123;
    v[825][7] = 219;
    v[826][7] = 129;
    v[827][7] = 155;
    v[828][7] = 207;
    v[829][7] = 177;
    v[830][7] = 9;
    v[831][7] = 49;
    v[832][7] = 181;
    v[833][7] = 231;
    v[834][7] = 33;
    v[835][7] = 233;
    v[836][7] = 67;
    v[837][7] = 155;
    v[838][7] = 41;
    v[839][7] = 9;
    v[840][7] = 95;
    v[841][7] = 123;
    v[842][7] = 65;
    v[843][7] = 117;
    v[844][7] = 249;
    v[845][7] = 85;
    v[846][7] = 169;
    v[847][7] = 129;
    v[848][7] = 241;
    v[849][7] = 173;
    v[850][7] = 251;
    v[851][7] = 225;
    v[852][7] = 147;
    v[853][7] = 165;
    v[854][7] = 69;
    v[855][7] = 81;
    v[856][7] = 239;
    v[857][7] = 95;
    v[858][7] = 23;
    v[859][7] = 83;
    v[860][7] = 227;
    v[861][7] = 249;
    v[862][7] = 143;
    v[863][7] = 171;
    v[864][7] = 193;
    v[865][7] = 9;
    v[866][7] = 21;
    v[867][7] = 57;
    v[868][7] = 73;
    v[869][7] = 97;
    v[870][7] = 57;
    v[871][7] = 29;
    v[872][7] = 239;
    v[873][7] = 151;
    v[874][7] = 159;
    v[875][7] = 191;
    v[876][7] = 47;
    v[877][7] = 51;
    v[878][7] = 1;
    v[879][7] = 223;
    v[880][7] = 251;
    v[881][7] = 251;
    v[882][7] = 151;
    v[883][7] = 41;
    v[884][7] = 119;
    v[885][7] = 127;
    v[886][7] = 131;
    v[887][7] = 33;
    v[888][7] = 209;
    v[889][7] = 123;
    v[890][7] = 53;
    v[891][7] = 241;
    v[892][7] = 25;
    v[893][7] = 31;
    v[894][7] = 183;
    v[895][7] = 107;
    v[896][7] = 25;
    v[897][7] = 115;
    v[898][7] = 39;
    v[899][7] = 11;
    v[900][7] = 213;
    v[901][7] = 239;
    v[902][7] = 219;
    v[903][7] = 109;
    v[904][7] = 185;
    v[905][7] = 35;
    v[906][7] = 133;
    v[907][7] = 123;
    v[908][7] = 185;
    v[909][7] = 27;
    v[910][7] = 55;
    v[911][7] = 245;
    v[912][7] = 61;
    v[913][7] = 75;
    v[914][7] = 205;
    v[915][7] = 213;
    v[916][7] = 169;
    v[917][7] = 163;
    v[918][7] = 63;
    v[919][7] = 55;
    v[920][7] = 49;
    v[921][7] = 83;
    v[922][7] = 195;
    v[923][7] = 51;
    v[924][7] = 31;
    v[925][7] = 41;
    v[926][7] = 15;
    v[927][7] = 203;
    v[928][7] = 41;
    v[929][7] = 63;
    v[930][7] = 127;
    v[931][7] = 161;
    v[932][7] = 5;
    v[933][7] = 143;
    v[934][7] = 7;
    v[935][7] = 199;
    v[936][7] = 251;
    v[937][7] = 95;
    v[938][7] = 75;
    v[939][7] = 101;
    v[940][7] = 15;
    v[941][7] = 43;
    v[942][7] = 237;
    v[943][7] = 197;
    v[944][7] = 117;
    v[945][7] = 167;
    v[946][7] = 155;
    v[947][7] = 21;
    v[948][7] = 83;
    v[949][7] = 205;
    v[950][7] = 255;
    v[951][7] = 49;
    v[952][7] = 101;
    v[953][7] = 213;
    v[954][7] = 237;
    v[955][7] = 135;
    v[956][7] = 135;
    v[957][7] = 21;
    v[958][7] = 73;
    v[959][7] = 93;
    v[960][7] = 115;
    v[961][7] = 7;
    v[962][7] = 85;
    v[963][7] = 223;
    v[964][7] = 237;
    v[965][7] = 79;
    v[966][7] = 89;
    v[967][7] = 5;
    v[968][7] = 57;
    v[969][7] = 239;
    v[970][7] = 67;
    v[971][7] = 65;
    v[972][7] = 201;
    v[973][7] = 155;
    v[974][7] = 71;
    v[975][7] = 85;
    v[976][7] = 195;
    v[977][7] = 89;
    v[978][7] = 181;
    v[979][7] = 119;
    v[980][7] = 135;
    v[981][7] = 147;
    v[982][7] = 237;
    v[983][7] = 173;
    v[984][7] = 41;
    v[985][7] = 155;
    v[986][7] = 67;
    v[987][7] = 113;
    v[988][7] = 111;
    v[989][7] = 21;
    v[990][7] = 183;
    v[991][7] = 23;
    v[992][7] = 103;
    v[993][7] = 207;
    v[994][7] = 253;
    v[995][7] = 69;
    v[996][7] = 219;
    v[997][7] = 205;
    v[998][7] = 195;
    v[999][7] = 43;
    v[1000][7] = 197;
    v[1001][7] = 229;
    v[1002][7] = 139;
    v[1003][7] = 177;
    v[1004][7] = 129;
    v[1005][7] = 69;
    v[1006][7] = 97;
    v[1007][7] = 201;
    v[1008][7] = 163;
    v[1009][7] = 189;
    v[1010][7] = 11;
    v[1011][7] = 99;
    v[1012][7] = 91;
    v[1013][7] = 253;
    v[1014][7] = 239;
    v[1015][7] = 91;
    v[1016][7] = 145;
    v[1017][7] = 19;
    v[1018][7] = 179;
    v[1019][7] = 231;
    v[1020][7] = 121;
    v[1021][7] = 7;
    v[1022][7] = 225;
    v[1023][7] = 237;
    v[1024][7] = 125;
    v[1025][7] = 191;
    v[1026][7] = 119;
    v[1027][7] = 59;
    v[1028][7] = 175;
    v[1029][7] = 237;
    v[1030][7] = 131;
    v[1031][7] = 79;
    v[1032][7] = 43;
    v[1033][7] = 45;
    v[1034][7] = 205;
    v[1035][7] = 199;
    v[1036][7] = 251;
    v[1037][7] = 153;
    v[1038][7] = 207;
    v[1039][7] = 37;
    v[1040][7] = 179;
    v[1041][7] = 113;
    v[1042][7] = 255;
    v[1043][7] = 107;
    v[1044][7] = 217;
    v[1045][7] = 61;
    v[1046][7] = 7;
    v[1047][7] = 181;
    v[1048][7] = 247;
    v[1049][7] = 31;
    v[1050][7] = 13;
    v[1051][7] = 113;
    v[1052][7] = 145;
    v[1053][7] = 107;
    v[1054][7] = 233;
    v[1055][7] = 233;
    v[1056][7] = 43;
    v[1057][7] = 79;
    v[1058][7] = 23;
    v[1059][7] = 169;
    v[1060][7] = 137;
    v[1061][7] = 129;
    v[1062][7] = 183;
    v[1063][7] = 53;
    v[1064][7] = 91;
    v[1065][7] = 55;
    v[1066][7] = 103;
    v[1067][7] = 223;
    v[1068][7] = 87;
    v[1069][7] = 177;
    v[1070][7] = 157;
    v[1071][7] = 79;
    v[1072][7] = 213;
    v[1073][7] = 139;
    v[1074][7] = 183;
    v[1075][7] = 231;
    v[1076][7] = 205;
    v[1077][7] = 143;
    v[1078][7] = 129;
    v[1079][7] = 243;
    v[1080][7] = 205;
    v[1081][7] = 93;
    v[1082][7] = 59;
    v[1083][7] = 15;
    v[1084][7] = 89;
    v[1085][7] = 9;
    v[1086][7] = 11;
    v[1087][7] = 47;
    v[1088][7] = 133;
    v[1089][7] = 227;
    v[1090][7] = 75;
    v[1091][7] = 9;
    v[1092][7] = 91;
    v[1093][7] = 19;
    v[1094][7] = 171;
    v[1095][7] = 163;
    v[1096][7] = 79;
    v[1097][7] = 7;
    v[1098][7] = 103;
    v[1099][7] = 5;
    v[1100][7] = 119;
    v[1101][7] = 155;
    v[1102][7] = 75;
    v[1103][7] = 11;
    v[1104][7] = 71;
    v[1105][7] = 95;
    v[1106][7] = 17;
    v[1107][7] = 13;
    v[1108][7] = 243;
    v[1109][7] = 207;
    v[1110][7] = 187;

    v[53][8] = 235;
    v[54][8] = 307;
    v[55][8] = 495;
    v[56][8] = 417;
    v[57][8] = 57;
    v[58][8] = 151;
    v[59][8] = 19;
    v[60][8] = 119;
    v[61][8] = 375;
    v[62][8] = 451;
    v[63][8] = 55;
    v[64][8] = 449;
    v[65][8] = 501;
    v[66][8] = 53;
    v[67][8] = 185;
    v[68][8] = 317;
    v[69][8] = 17;
    v[70][8] = 21;
    v[71][8] = 487;
    v[72][8] = 13;
    v[73][8] = 347;
    v[74][8] = 393;
    v[75][8] = 15;
    v[76][8] = 391;
    v[77][8] = 307;
    v[78][8] = 189;
    v[79][8] = 381;
    v[80][8] = 71;
    v[81][8] = 163;
    v[82][8] = 99;
    v[83][8] = 467;
    v[84][8] = 167;
    v[85][8] = 433;
    v[86][8] = 337;
    v[87][8] = 257;
    v[88][8] = 179;
    v[89][8] = 47;
    v[90][8] = 385;
    v[91][8] = 23;
    v[92][8] = 117;
    v[93][8] = 369;
    v[94][8] = 425;
    v[95][8] = 207;
    v[96][8] = 433;
    v[97][8] = 301;
    v[98][8] = 147;
    v[99][8] = 333;
    v[100][8] = 85;
    v[101][8] = 221;
    v[102][8] = 423;
    v[103][8] = 49;
    v[104][8] = 3;
    v[105][8] = 43;
    v[106][8] = 229;
    v[107][8] = 227;
    v[108][8] = 201;
    v[109][8] = 383;
    v[110][8] = 281;
    v[111][8] = 229;
    v[112][8] = 207;
    v[113][8] = 21;
    v[114][8] = 343;
    v[115][8] = 251;
    v[116][8] = 397;
    v[117][8] = 173;
    v[118][8] = 507;
    v[119][8] = 421;
    v[120][8] = 443;
    v[121][8] = 399;
    v[122][8] = 53;
    v[123][8] = 345;
    v[124][8] = 77;
    v[125][8] = 385;
    v[126][8] = 317;
    v[127][8] = 155;
    v[128][8] = 187;
    v[129][8] = 269;
    v[130][8] = 501;
    v[131][8] = 19;
    v[132][8] = 169;
    v[133][8] = 235;
    v[134][8] = 415;
    v[135][8] = 61;
    v[136][8] = 247;
    v[137][8] = 183;
    v[138][8] = 5;
    v[139][8] = 257;
    v[140][8] = 401;
    v[141][8] = 451;
    v[142][8] = 95;
    v[143][8] = 455;
    v[144][8] = 49;
    v[145][8] = 489;
    v[146][8] = 75;
    v[147][8] = 459;
    v[148][8] = 377;
    v[149][8] = 87;
    v[150][8] = 463;
    v[151][8] = 155;
    v[152][8] = 233;
    v[153][8] = 115;
    v[154][8] = 429;
    v[155][8] = 211;
    v[156][8] = 419;
    v[157][8] = 143;
    v[158][8] = 487;
    v[159][8] = 195;
    v[160][8] = 209;
    v[161][8] = 461;
    v[162][8] = 193;
    v[163][8] = 157;
    v[164][8] = 193;
    v[165][8] = 363;
    v[166][8] = 181;
    v[167][8] = 271;
    v[168][8] = 445;
    v[169][8] = 381;
    v[170][8] = 231;
    v[171][8] = 135;
    v[172][8] = 327;
    v[173][8] = 403;
    v[174][8] = 171;
    v[175][8] = 197;
    v[176][8] = 181;
    v[177][8] = 343;
    v[178][8] = 113;
    v[179][8] = 313;
    v[180][8] = 393;
    v[181][8] = 311;
    v[182][8] = 415;
    v[183][8] = 267;
    v[184][8] = 247;
    v[185][8] = 425;
    v[186][8] = 233;
    v[187][8] = 289;
    v[188][8] = 55;
    v[189][8] = 39;
    v[190][8] = 247;
    v[191][8] = 327;
    v[192][8] = 141;
    v[193][8] = 5;
    v[194][8] = 189;
    v[195][8] = 183;
    v[196][8] = 27;
    v[197][8] = 337;
    v[198][8] = 341;
    v[199][8] = 327;
    v[200][8] = 87;
    v[201][8] = 429;
    v[202][8] = 357;
    v[203][8] = 265;
    v[204][8] = 251;
    v[205][8] = 437;
    v[206][8] = 201;
    v[207][8] = 29;
    v[208][8] = 339;
    v[209][8] = 257;
    v[210][8] = 377;
    v[211][8] = 17;
    v[212][8] = 53;
    v[213][8] = 327;
    v[214][8] = 47;
    v[215][8] = 375;
    v[216][8] = 393;
    v[217][8] = 369;
    v[218][8] = 403;
    v[219][8] = 125;
    v[220][8] = 429;
    v[221][8] = 257;
    v[222][8] = 157;
    v[223][8] = 217;
    v[224][8] = 85;
    v[225][8] = 267;
    v[226][8] = 117;
    v[227][8] = 337;
    v[228][8] = 447;
    v[229][8] = 219;
    v[230][8] = 501;
    v[231][8] = 41;
    v[232][8] = 41;
    v[233][8] = 193;
    v[234][8] = 509;
    v[235][8] = 131;
    v[236][8] = 207;
    v[237][8] = 505;
    v[238][8] = 421;
    v[239][8] = 149;
    v[240][8] = 111;
    v[241][8] = 177;
    v[242][8] = 167;
    v[243][8] = 223;
    v[244][8] = 291;
    v[245][8] = 91;
    v[246][8] = 29;
    v[247][8] = 305;
    v[248][8] = 151;
    v[249][8] = 177;
    v[250][8] = 337;
    v[251][8] = 183;
    v[252][8] = 361;
    v[253][8] = 435;
    v[254][8] = 307;
    v[255][8] = 507;
    v[256][8] = 77;
    v[257][8] = 181;
    v[258][8] = 507;
    v[259][8] = 315;
    v[260][8] = 145;
    v[261][8] = 423;
    v[262][8] = 71;
    v[263][8] = 103;
    v[264][8] = 493;
    v[265][8] = 271;
    v[266][8] = 469;
    v[267][8] = 339;
    v[268][8] = 237;
    v[269][8] = 437;
    v[270][8] = 483;
    v[271][8] = 31;
    v[272][8] = 219;
    v[273][8] = 61;
    v[274][8] = 131;
    v[275][8] = 391;
    v[276][8] = 233;
    v[277][8] = 219;
    v[278][8] = 69;
    v[279][8] = 57;
    v[280][8] = 459;
    v[281][8] = 225;
    v[282][8] = 421;
    v[283][8] = 7;
    v[284][8] = 461;
    v[285][8] = 111;
    v[286][8] = 451;
    v[287][8] = 277;
    v[288][8] = 185;
    v[289][8] = 193;
    v[290][8] = 125;
    v[291][8] = 251;
    v[292][8] = 199;
    v[293][8] = 73;
    v[294][8] = 71;
    v[295][8] = 7;
    v[296][8] = 409;
    v[297][8] = 417;
    v[298][8] = 149;
    v[299][8] = 193;
    v[300][8] = 53;
    v[301][8] = 437;
    v[302][8] = 29;
    v[303][8] = 467;
    v[304][8] = 229;
    v[305][8] = 31;
    v[306][8] = 35;
    v[307][8] = 75;
    v[308][8] = 105;
    v[309][8] = 503;
    v[310][8] = 75;
    v[311][8] = 317;
    v[312][8] = 401;
    v[313][8] = 367;
    v[314][8] = 131;
    v[315][8] = 365;
    v[316][8] = 441;
    v[317][8] = 433;
    v[318][8] = 93;
    v[319][8] = 377;
    v[320][8] = 405;
    v[321][8] = 465;
    v[322][8] = 259;
    v[323][8] = 283;
    v[324][8] = 443;
    v[325][8] = 143;
    v[326][8] = 445;
    v[327][8] = 3;
    v[328][8] = 461;
    v[329][8] = 329;
    v[330][8] = 309;
    v[331][8] = 77;
    v[332][8] = 323;
    v[333][8] = 155;
    v[334][8] = 347;
    v[335][8] = 45;
    v[336][8] = 381;
    v[337][8] = 315;
    v[338][8] = 463;
    v[339][8] = 207;
    v[340][8] = 321;
    v[341][8] = 157;
    v[342][8] = 109;
    v[343][8] = 479;
    v[344][8] = 313;
    v[345][8] = 345;
    v[346][8] = 167;
    v[347][8] = 439;
    v[348][8] = 307;
    v[349][8] = 235;
    v[350][8] = 473;
    v[351][8] = 79;
    v[352][8] = 101;
    v[353][8] = 245;
    v[354][8] = 19;
    v[355][8] = 381;
    v[356][8] = 251;
    v[357][8] = 35;
    v[358][8] = 25;
    v[359][8] = 107;
    v[360][8] = 187;
    v[361][8] = 115;
    v[362][8] = 113;
    v[363][8] = 321;
    v[364][8] = 115;
    v[365][8] = 445;
    v[366][8] = 61;
    v[367][8] = 77;
    v[368][8] = 293;
    v[369][8] = 405;
    v[370][8] = 13;
    v[371][8] = 53;
    v[372][8] = 17;
    v[373][8] = 171;
    v[374][8] = 299;
    v[375][8] = 41;
    v[376][8] = 79;
    v[377][8] = 3;
    v[378][8] = 485;
    v[379][8] = 331;
    v[380][8] = 13;
    v[381][8] = 257;
    v[382][8] = 59;
    v[383][8] = 201;
    v[384][8] = 497;
    v[385][8] = 81;
    v[386][8] = 451;
    v[387][8] = 199;
    v[388][8] = 171;
    v[389][8] = 81;
    v[390][8] = 253;
    v[391][8] = 365;
    v[392][8] = 75;
    v[393][8] = 451;
    v[394][8] = 149;
    v[395][8] = 483;
    v[396][8] = 81;
    v[397][8] = 453;
    v[398][8] = 469;
    v[399][8] = 485;
    v[400][8] = 305;
    v[401][8] = 163;
    v[402][8] = 401;
    v[403][8] = 15;
    v[404][8] = 91;
    v[405][8] = 3;
    v[406][8] = 129;
    v[407][8] = 35;
    v[408][8] = 239;
    v[409][8] = 355;
    v[410][8] = 211;
    v[411][8] = 387;
    v[412][8] = 101;
    v[413][8] = 299;
    v[414][8] = 67;
    v[415][8] = 375;
    v[416][8] = 405;
    v[417][8] = 357;
    v[418][8] = 267;
    v[419][8] = 363;
    v[420][8] = 79;
    v[421][8] = 83;
    v[422][8] = 437;
    v[423][8] = 457;
    v[424][8] = 39;
    v[425][8] = 97;
    v[426][8] = 473;
    v[427][8] = 289;
    v[428][8] = 179;
    v[429][8] = 57;
    v[430][8] = 23;
    v[431][8] = 49;
    v[432][8] = 79;
    v[433][8] = 71;
    v[434][8] = 341;
    v[435][8] = 287;
    v[436][8] = 95;
    v[437][8] = 229;
    v[438][8] = 271;
    v[439][8] = 475;
    v[440][8] = 49;
    v[441][8] = 241;
    v[442][8] = 261;
    v[443][8] = 495;
    v[444][8] = 353;
    v[445][8] = 381;
    v[446][8] = 13;
    v[447][8] = 291;
    v[448][8] = 37;
    v[449][8] = 251;
    v[450][8] = 105;
    v[451][8] = 399;
    v[452][8] = 81;
    v[453][8] = 89;
    v[454][8] = 265;
    v[455][8] = 507;
    v[456][8] = 205;
    v[457][8] = 145;
    v[458][8] = 331;
    v[459][8] = 129;
    v[460][8] = 119;
    v[461][8] = 503;
    v[462][8] = 249;
    v[463][8] = 1;
    v[464][8] = 289;
    v[465][8] = 463;
    v[466][8] = 163;
    v[467][8] = 443;
    v[468][8] = 63;
    v[469][8] = 123;
    v[470][8] = 361;
    v[471][8] = 261;
    v[472][8] = 49;
    v[473][8] = 429;
    v[474][8] = 137;
    v[475][8] = 355;
    v[476][8] = 175;
    v[477][8] = 507;
    v[478][8] = 59;
    v[479][8] = 277;
    v[480][8] = 391;
    v[481][8] = 25;
    v[482][8] = 185;
    v[483][8] = 381;
    v[484][8] = 197;
    v[485][8] = 39;
    v[486][8] = 5;
    v[487][8] = 429;
    v[488][8] = 119;
    v[489][8] = 247;
    v[490][8] = 177;
    v[491][8] = 329;
    v[492][8] = 465;
    v[493][8] = 421;
    v[494][8] = 271;
    v[495][8] = 467;
    v[496][8] = 151;
    v[497][8] = 45;
    v[498][8] = 429;
    v[499][8] = 137;
    v[500][8] = 471;
    v[501][8] = 11;
    v[502][8] = 17;
    v[503][8] = 409;
    v[504][8] = 347;
    v[505][8] = 199;
    v[506][8] = 463;
    v[507][8] = 177;
    v[508][8] = 11;
    v[509][8] = 51;
    v[510][8] = 361;
    v[511][8] = 95;
    v[512][8] = 497;
    v[513][8] = 163;
    v[514][8] = 351;
    v[515][8] = 127;
    v[516][8] = 395;
    v[517][8] = 511;
    v[518][8] = 327;
    v[519][8] = 353;
    v[520][8] = 49;
    v[521][8] = 105;
    v[522][8] = 151;
    v[523][8] = 321;
    v[524][8] = 331;
    v[525][8] = 329;
    v[526][8] = 509;
    v[527][8] = 107;
    v[528][8] = 109;
    v[529][8] = 303;
    v[530][8] = 467;
    v[531][8] = 287;
    v[532][8] = 161;
    v[533][8] = 45;
    v[534][8] = 385;
    v[535][8] = 289;
    v[536][8] = 363;
    v[537][8] = 331;
    v[538][8] = 265;
    v[539][8] = 407;
    v[540][8] = 37;
    v[541][8] = 433;
    v[542][8] = 315;
    v[543][8] = 343;
    v[544][8] = 63;
    v[545][8] = 51;
    v[546][8] = 185;
    v[547][8] = 71;
    v[548][8] = 27;
    v[549][8] = 267;
    v[550][8] = 503;
    v[551][8] = 239;
    v[552][8] = 293;
    v[553][8] = 245;
    v[554][8] = 281;
    v[555][8] = 297;
    v[556][8] = 75;
    v[557][8] = 461;
    v[558][8] = 371;
    v[559][8] = 129;
    v[560][8] = 189;
    v[561][8] = 189;
    v[562][8] = 339;
    v[563][8] = 287;
    v[564][8] = 111;
    v[565][8] = 111;
    v[566][8] = 379;
    v[567][8] = 93;
    v[568][8] = 27;
    v[569][8] = 185;
    v[570][8] = 347;
    v[571][8] = 337;
    v[572][8] = 247;
    v[573][8] = 507;
    v[574][8] = 161;
    v[575][8] = 231;
    v[576][8] = 43;
    v[577][8] = 499;
    v[578][8] = 73;
    v[579][8] = 327;
    v[580][8] = 263;
    v[581][8] = 331;
    v[582][8] = 249;
    v[583][8] = 493;
    v[584][8] = 37;
    v[585][8] = 25;
    v[586][8] = 115;
    v[587][8] = 3;
    v[588][8] = 167;
    v[589][8] = 197;
    v[590][8] = 127;
    v[591][8] = 357;
    v[592][8] = 497;
    v[593][8] = 103;
    v[594][8] = 125;
    v[595][8] = 191;
    v[596][8] = 165;
    v[597][8] = 55;
    v[598][8] = 101;
    v[599][8] = 95;
    v[600][8] = 79;
    v[601][8] = 351;
    v[602][8] = 341;
    v[603][8] = 43;
    v[604][8] = 125;
    v[605][8] = 135;
    v[606][8] = 173;
    v[607][8] = 289;
    v[608][8] = 373;
    v[609][8] = 133;
    v[610][8] = 421;
    v[611][8] = 241;
    v[612][8] = 281;
    v[613][8] = 213;
    v[614][8] = 177;
    v[615][8] = 363;
    v[616][8] = 151;
    v[617][8] = 227;
    v[618][8] = 145;
    v[619][8] = 363;
    v[620][8] = 239;
    v[621][8] = 431;
    v[622][8] = 81;
    v[623][8] = 397;
    v[624][8] = 241;
    v[625][8] = 67;
    v[626][8] = 291;
    v[627][8] = 255;
    v[628][8] = 405;
    v[629][8] = 421;
    v[630][8] = 399;
    v[631][8] = 75;
    v[632][8] = 399;
    v[633][8] = 105;
    v[634][8] = 329;
    v[635][8] = 41;
    v[636][8] = 425;
    v[637][8] = 7;
    v[638][8] = 283;
    v[639][8] = 375;
    v[640][8] = 475;
    v[641][8] = 427;
    v[642][8] = 277;
    v[643][8] = 209;
    v[644][8] = 411;
    v[645][8] = 3;
    v[646][8] = 137;
    v[647][8] = 195;
    v[648][8] = 289;
    v[649][8] = 509;
    v[650][8] = 121;
    v[651][8] = 55;
    v[652][8] = 147;
    v[653][8] = 275;
    v[654][8] = 251;
    v[655][8] = 19;
    v[656][8] = 129;
    v[657][8] = 285;
    v[658][8] = 415;
    v[659][8] = 487;
    v[660][8] = 491;
    v[661][8] = 193;
    v[662][8] = 219;
    v[663][8] = 403;
    v[664][8] = 23;
    v[665][8] = 97;
    v[666][8] = 65;
    v[667][8] = 285;
    v[668][8] = 75;
    v[669][8] = 21;
    v[670][8] = 373;
    v[671][8] = 261;
    v[672][8] = 339;
    v[673][8] = 239;
    v[674][8] = 495;
    v[675][8] = 415;
    v[676][8] = 333;
    v[677][8] = 107;
    v[678][8] = 435;
    v[679][8] = 297;
    v[680][8] = 213;
    v[681][8] = 149;
    v[682][8] = 463;
    v[683][8] = 199;
    v[684][8] = 323;
    v[685][8] = 45;
    v[686][8] = 19;
    v[687][8] = 301;
    v[688][8] = 121;
    v[689][8] = 499;
    v[690][8] = 187;
    v[691][8] = 229;
    v[692][8] = 63;
    v[693][8] = 425;
    v[694][8] = 99;
    v[695][8] = 281;
    v[696][8] = 35;
    v[697][8] = 125;
    v[698][8] = 349;
    v[699][8] = 87;
    v[700][8] = 101;
    v[701][8] = 59;
    v[702][8] = 195;
    v[703][8] = 511;
    v[704][8] = 355;
    v[705][8] = 73;
    v[706][8] = 263;
    v[707][8] = 243;
    v[708][8] = 101;
    v[709][8] = 165;
    v[710][8] = 141;
    v[711][8] = 11;
    v[712][8] = 389;
    v[713][8] = 219;
    v[714][8] = 187;
    v[715][8] = 449;
    v[716][8] = 447;
    v[717][8] = 393;
    v[718][8] = 477;
    v[719][8] = 305;
    v[720][8] = 221;
    v[721][8] = 51;
    v[722][8] = 355;
    v[723][8] = 209;
    v[724][8] = 499;
    v[725][8] = 479;
    v[726][8] = 265;
    v[727][8] = 377;
    v[728][8] = 145;
    v[729][8] = 411;
    v[730][8] = 173;
    v[731][8] = 11;
    v[732][8] = 433;
    v[733][8] = 483;
    v[734][8] = 135;
    v[735][8] = 385;
    v[736][8] = 341;
    v[737][8] = 89;
    v[738][8] = 209;
    v[739][8] = 391;
    v[740][8] = 33;
    v[741][8] = 395;
    v[742][8] = 319;
    v[743][8] = 451;
    v[744][8] = 119;
    v[745][8] = 341;
    v[746][8] = 227;
    v[747][8] = 375;
    v[748][8] = 61;
    v[749][8] = 331;
    v[750][8] = 493;
    v[751][8] = 411;
    v[752][8] = 293;
    v[753][8] = 47;
    v[754][8] = 203;
    v[755][8] = 375;
    v[756][8] = 167;
    v[757][8] = 395;
    v[758][8] = 155;
    v[759][8] = 5;
    v[760][8] = 237;
    v[761][8] = 361;
    v[762][8] = 489;
    v[763][8] = 127;
    v[764][8] = 21;
    v[765][8] = 345;
    v[766][8] = 101;
    v[767][8] = 371;
    v[768][8] = 233;
    v[769][8] = 431;
    v[770][8] = 109;
    v[771][8] = 119;
    v[772][8] = 277;
    v[773][8] = 125;
    v[774][8] = 263;
    v[775][8] = 73;
    v[776][8] = 135;
    v[777][8] = 123;
    v[778][8] = 83;
    v[779][8] = 123;
    v[780][8] = 405;
    v[781][8] = 69;
    v[782][8] = 75;
    v[783][8] = 287;
    v[784][8] = 401;
    v[785][8] = 23;
    v[786][8] = 283;
    v[787][8] = 393;
    v[788][8] = 41;
    v[789][8] = 379;
    v[790][8] = 431;
    v[791][8] = 11;
    v[792][8] = 475;
    v[793][8] = 505;
    v[794][8] = 19;
    v[795][8] = 365;
    v[796][8] = 265;
    v[797][8] = 271;
    v[798][8] = 499;
    v[799][8] = 489;
    v[800][8] = 443;
    v[801][8] = 165;
    v[802][8] = 91;
    v[803][8] = 83;
    v[804][8] = 291;
    v[805][8] = 319;
    v[806][8] = 199;
    v[807][8] = 107;
    v[808][8] = 245;
    v[809][8] = 389;
    v[810][8] = 143;
    v[811][8] = 137;
    v[812][8] = 89;
    v[813][8] = 125;
    v[814][8] = 281;
    v[815][8] = 381;
    v[816][8] = 215;
    v[817][8] = 131;
    v[818][8] = 299;
    v[819][8] = 249;
    v[820][8] = 375;
    v[821][8] = 455;
    v[822][8] = 43;
    v[823][8] = 73;
    v[824][8] = 281;
    v[825][8] = 217;
    v[826][8] = 297;
    v[827][8] = 229;
    v[828][8] = 431;
    v[829][8] = 357;
    v[830][8] = 81;
    v[831][8] = 357;
    v[832][8] = 171;
    v[833][8] = 451;
    v[834][8] = 481;
    v[835][8] = 13;
    v[836][8] = 387;
    v[837][8] = 491;
    v[838][8] = 489;
    v[839][8] = 439;
    v[840][8] = 385;
    v[841][8] = 487;
    v[842][8] = 177;
    v[843][8] = 393;
    v[844][8] = 33;
    v[845][8] = 71;
    v[846][8] = 375;
    v[847][8] = 443;
    v[848][8] = 129;
    v[849][8] = 407;
    v[850][8] = 395;
    v[851][8] = 127;
    v[852][8] = 65;
    v[853][8] = 333;
    v[854][8] = 309;
    v[855][8] = 119;
    v[856][8] = 197;
    v[857][8] = 435;
    v[858][8] = 497;
    v[859][8] = 373;
    v[860][8] = 71;
    v[861][8] = 379;
    v[862][8] = 509;
    v[863][8] = 387;
    v[864][8] = 159;
    v[865][8] = 265;
    v[866][8] = 477;
    v[867][8] = 463;
    v[868][8] = 449;
    v[869][8] = 47;
    v[870][8] = 353;
    v[871][8] = 249;
    v[872][8] = 335;
    v[873][8] = 505;
    v[874][8] = 89;
    v[875][8] = 141;
    v[876][8] = 55;
    v[877][8] = 235;
    v[878][8] = 187;
    v[879][8] = 87;
    v[880][8] = 363;
    v[881][8] = 93;
    v[882][8] = 363;
    v[883][8] = 101;
    v[884][8] = 67;
    v[885][8] = 215;
    v[886][8] = 321;
    v[887][8] = 331;
    v[888][8] = 305;
    v[889][8] = 261;
    v[890][8] = 411;
    v[891][8] = 491;
    v[892][8] = 479;
    v[893][8] = 65;
    v[894][8] = 307;
    v[895][8] = 469;
    v[896][8] = 415;
    v[897][8] = 131;
    v[898][8] = 315;
    v[899][8] = 487;
    v[900][8] = 83;
    v[901][8] = 455;
    v[902][8] = 19;
    v[903][8] = 113;
    v[904][8] = 163;
    v[905][8] = 503;
    v[906][8] = 99;
    v[907][8] = 499;
    v[908][8] = 251;
    v[909][8] = 239;
    v[910][8] = 81;
    v[911][8] = 167;
    v[912][8] = 391;
    v[913][8] = 255;
    v[914][8] = 317;
    v[915][8] = 363;
    v[916][8] = 359;
    v[917][8] = 395;
    v[918][8] = 419;
    v[919][8] = 307;
    v[920][8] = 251;
    v[921][8] = 267;
    v[922][8] = 171;
    v[923][8] = 461;
    v[924][8] = 183;
    v[925][8] = 465;
    v[926][8] = 165;
    v[927][8] = 163;
    v[928][8] = 293;
    v[929][8] = 477;
    v[930][8] = 223;
    v[931][8] = 403;
    v[932][8] = 389;
    v[933][8] = 97;
    v[934][8] = 335;
    v[935][8] = 357;
    v[936][8] = 297;
    v[937][8] = 19;
    v[938][8] = 469;
    v[939][8] = 501;
    v[940][8] = 249;
    v[941][8] = 85;
    v[942][8] = 213;
    v[943][8] = 311;
    v[944][8] = 265;
    v[945][8] = 379;
    v[946][8] = 297;
    v[947][8] = 283;
    v[948][8] = 393;
    v[949][8] = 449;
    v[950][8] = 463;
    v[951][8] = 289;
    v[952][8] = 159;
    v[953][8] = 289;
    v[954][8] = 499;
    v[955][8] = 407;
    v[956][8] = 129;
    v[957][8] = 137;
    v[958][8] = 221;
    v[959][8] = 43;
    v[960][8] = 89;
    v[961][8] = 403;
    v[962][8] = 271;
    v[963][8] = 75;
    v[964][8] = 83;
    v[965][8] = 445;
    v[966][8] = 453;
    v[967][8] = 389;
    v[968][8] = 149;
    v[969][8] = 143;
    v[970][8] = 423;
    v[971][8] = 499;
    v[972][8] = 317;
    v[973][8] = 445;
    v[974][8] = 157;
    v[975][8] = 137;
    v[976][8] = 453;
    v[977][8] = 163;
    v[978][8] = 87;
    v[979][8] = 23;
    v[980][8] = 391;
    v[981][8] = 119;
    v[982][8] = 427;
    v[983][8] = 323;
    v[984][8] = 173;
    v[985][8] = 89;
    v[986][8] = 259;
    v[987][8] = 377;
    v[988][8] = 511;
    v[989][8] = 249;
    v[990][8] = 31;
    v[991][8] = 363;
    v[992][8] = 229;
    v[993][8] = 353;
    v[994][8] = 329;
    v[995][8] = 493;
    v[996][8] = 427;
    v[997][8] = 57;
    v[998][8] = 205;
    v[999][8] = 389;
    v[1000][8] = 91;
    v[1001][8] = 83;
    v[1002][8] = 13;
    v[1003][8] = 219;
    v[1004][8] = 439;
    v[1005][8] = 45;
    v[1006][8] = 35;
    v[1007][8] = 371;
    v[1008][8] = 441;
    v[1009][8] = 17;
    v[1010][8] = 267;
    v[1011][8] = 501;
    v[1012][8] = 53;
    v[1013][8] = 25;
    v[1014][8] = 333;
    v[1015][8] = 17;
    v[1016][8] = 201;
    v[1017][8] = 475;
    v[1018][8] = 257;
    v[1019][8] = 417;
    v[1020][8] = 345;
    v[1021][8] = 381;
    v[1022][8] = 377;
    v[1023][8] = 55;
    v[1024][8] = 403;
    v[1025][8] = 77;
    v[1026][8] = 389;
    v[1027][8] = 347;
    v[1028][8] = 363;
    v[1029][8] = 211;
    v[1030][8] = 413;
    v[1031][8] = 419;
    v[1032][8] = 5;
    v[1033][8] = 167;
    v[1034][8] = 219;
    v[1035][8] = 201;
    v[1036][8] = 285;
    v[1037][8] = 425;
    v[1038][8] = 11;
    v[1039][8] = 77;
    v[1040][8] = 269;
    v[1041][8] = 489;
    v[1042][8] = 281;
    v[1043][8] = 403;
    v[1044][8] = 79;
    v[1045][8] = 425;
    v[1046][8] = 125;
    v[1047][8] = 81;
    v[1048][8] = 331;
    v[1049][8] = 437;
    v[1050][8] = 271;
    v[1051][8] = 397;
    v[1052][8] = 299;
    v[1053][8] = 475;
    v[1054][8] = 271;
    v[1055][8] = 249;
    v[1056][8] = 413;
    v[1057][8] = 233;
    v[1058][8] = 261;
    v[1059][8] = 495;
    v[1060][8] = 171;
    v[1061][8] = 69;
    v[1062][8] = 27;
    v[1063][8] = 409;
    v[1064][8] = 21;
    v[1065][8] = 421;
    v[1066][8] = 367;
    v[1067][8] = 81;
    v[1068][8] = 483;
    v[1069][8] = 255;
    v[1070][8] = 15;
    v[1071][8] = 219;
    v[1072][8] = 365;
    v[1073][8] = 497;
    v[1074][8] = 181;
    v[1075][8] = 75;
    v[1076][8] = 431;
    v[1077][8] = 99;
    v[1078][8] = 325;
    v[1079][8] = 407;
    v[1080][8] = 229;
    v[1081][8] = 281;
    v[1082][8] = 63;
    v[1083][8] = 83;
    v[1084][8] = 493;
    v[1085][8] = 5;
    v[1086][8] = 113;
    v[1087][8] = 15;
    v[1088][8] = 271;
    v[1089][8] = 37;
    v[1090][8] = 87;
    v[1091][8] = 451;
    v[1092][8] = 299;
    v[1093][8] = 83;
    v[1094][8] = 451;
    v[1095][8] = 311;
    v[1096][8] = 441;
    v[1097][8] = 47;
    v[1098][8] = 455;
    v[1099][8] = 47;
    v[1100][8] = 253;
    v[1101][8] = 13;
    v[1102][8] = 109;
    v[1103][8] = 369;
    v[1104][8] = 347;
    v[1105][8] = 11;
    v[1106][8] = 409;
    v[1107][8] = 275;
    v[1108][8] = 63;
    v[1109][8] = 441;
    v[1110][8] = 15;

    v[101][9] = 519;
    v[102][9] = 307;
    v[103][9] = 931;
    v[104][9] = 1023;
    v[105][9] = 517;
    v[106][9] = 771;
    v[107][9] = 151;
    v[108][9] = 1023;
    v[109][9] = 539;
    v[110][9] = 725;
    v[111][9] = 45;
    v[112][9] = 927;
    v[113][9] = 707;
    v[114][9] = 29;
    v[115][9] = 125;
    v[116][9] = 371;
    v[117][9] = 275;
    v[118][9] = 279;
    v[119][9] = 817;
    v[120][9] = 389;
    v[121][9] = 453;
    v[122][9] = 989;
    v[123][9] = 1015;
    v[124][9] = 29;
    v[125][9] = 169;
    v[126][9] = 743;
    v[127][9] = 99;
    v[128][9] = 923;
    v[129][9] = 981;
    v[130][9] = 181;
    v[131][9] = 693;
    v[132][9] = 309;
    v[133][9] = 227;
    v[134][9] = 111;
    v[135][9] = 219;
    v[136][9] = 897;
    v[137][9] = 377;
    v[138][9] = 425;
    v[139][9] = 609;
    v[140][9] = 227;
    v[141][9] = 19;
    v[142][9] = 221;
    v[143][9] = 143;
    v[144][9] = 581;
    v[145][9] = 147;
    v[146][9] = 919;
    v[147][9] = 127;
    v[148][9] = 725;
    v[149][9] = 793;
    v[150][9] = 289;
    v[151][9] = 411;
    v[152][9] = 835;
    v[153][9] = 921;
    v[154][9] = 957;
    v[155][9] = 443;
    v[156][9] = 349;
    v[157][9] = 813;
    v[158][9] = 5;
    v[159][9] = 105;
    v[160][9] = 457;
    v[161][9] = 393;
    v[162][9] = 539;
    v[163][9] = 101;
    v[164][9] = 197;
    v[165][9] = 697;
    v[166][9] = 27;
    v[167][9] = 343;
    v[168][9] = 515;
    v[169][9] = 69;
    v[170][9] = 485;
    v[171][9] = 383;
    v[172][9] = 855;
    v[173][9] = 693;
    v[174][9] = 133;
    v[175][9] = 87;
    v[176][9] = 743;
    v[177][9] = 747;
    v[178][9] = 475;
    v[179][9] = 87;
    v[180][9] = 469;
    v[181][9] = 763;
    v[182][9] = 721;
    v[183][9] = 345;
    v[184][9] = 479;
    v[185][9] = 965;
    v[186][9] = 527;
    v[187][9] = 121;
    v[188][9] = 271;
    v[189][9] = 353;
    v[190][9] = 467;
    v[191][9] = 177;
    v[192][9] = 245;
    v[193][9] = 627;
    v[194][9] = 113;
    v[195][9] = 357;
    v[196][9] = 7;
    v[197][9] = 691;
    v[198][9] = 725;
    v[199][9] = 355;
    v[200][9] = 889;
    v[201][9] = 635;
    v[202][9] = 737;
    v[203][9] = 429;
    v[204][9] = 545;
    v[205][9] = 925;
    v[206][9] = 357;
    v[207][9] = 873;
    v[208][9] = 187;
    v[209][9] = 351;
    v[210][9] = 677;
    v[211][9] = 999;
    v[212][9] = 921;
    v[213][9] = 477;
    v[214][9] = 233;
    v[215][9] = 765;
    v[216][9] = 495;
    v[217][9] = 81;
    v[218][9] = 953;
    v[219][9] = 479;
    v[220][9] = 89;
    v[221][9] = 173;
    v[222][9] = 473;
    v[223][9] = 131;
    v[224][9] = 961;
    v[225][9] = 411;
    v[226][9] = 291;
    v[227][9] = 967;
    v[228][9] = 65;
    v[229][9] = 511;
    v[230][9] = 13;
    v[231][9] = 805;
    v[232][9] = 945;
    v[233][9] = 369;
    v[234][9] = 827;
    v[235][9] = 295;
    v[236][9] = 163;
    v[237][9] = 835;
    v[238][9] = 259;
    v[239][9] = 207;
    v[240][9] = 331;
    v[241][9] = 29;
    v[242][9] = 315;
    v[243][9] = 999;
    v[244][9] = 133;
    v[245][9] = 967;
    v[246][9] = 41;
    v[247][9] = 117;
    v[248][9] = 677;
    v[249][9] = 471;
    v[250][9] = 717;
    v[251][9] = 881;
    v[252][9] = 755;
    v[253][9] = 351;
    v[254][9] = 723;
    v[255][9] = 259;
    v[256][9] = 879;
    v[257][9] = 455;
    v[258][9] = 721;
    v[259][9] = 289;
    v[260][9] = 149;
    v[261][9] = 199;
    v[262][9] = 805;
    v[263][9] = 987;
    v[264][9] = 851;
    v[265][9] = 423;
    v[266][9] = 597;
    v[267][9] = 129;
    v[268][9] = 11;
    v[269][9] = 733;
    v[270][9] = 549;
    v[271][9] = 153;
    v[272][9] = 285;
    v[273][9] = 451;
    v[274][9] = 559;
    v[275][9] = 377;
    v[276][9] = 109;
    v[277][9] = 357;
    v[278][9] = 143;
    v[279][9] = 693;
    v[280][9] = 615;
    v[281][9] = 677;
    v[282][9] = 701;
    v[283][9] = 475;
    v[284][9] = 767;
    v[285][9] = 85;
    v[286][9] = 229;
    v[287][9] = 509;
    v[288][9] = 547;
    v[289][9] = 151;
    v[290][9] = 389;
    v[291][9] = 711;
    v[292][9] = 785;
    v[293][9] = 657;
    v[294][9] = 319;
    v[295][9] = 509;
    v[296][9] = 99;
    v[297][9] = 1007;
    v[298][9] = 775;
    v[299][9] = 359;
    v[300][9] = 697;
    v[301][9] = 677;
    v[302][9] = 85;
    v[303][9] = 497;
    v[304][9] = 105;
    v[305][9] = 615;
    v[306][9] = 891;
    v[307][9] = 71;
    v[308][9] = 449;
    v[309][9] = 835;
    v[310][9] = 609;
    v[311][9] = 377;
    v[312][9] = 693;
    v[313][9] = 665;
    v[314][9] = 627;
    v[315][9] = 215;
    v[316][9] = 911;
    v[317][9] = 503;
    v[318][9] = 729;
    v[319][9] = 131;
    v[320][9] = 19;
    v[321][9] = 895;
    v[322][9] = 199;
    v[323][9] = 161;
    v[324][9] = 239;
    v[325][9] = 633;
    v[326][9] = 1013;
    v[327][9] = 537;
    v[328][9] = 255;
    v[329][9] = 23;
    v[330][9] = 149;
    v[331][9] = 679;
    v[332][9] = 1021;
    v[333][9] = 595;
    v[334][9] = 199;
    v[335][9] = 557;
    v[336][9] = 659;
    v[337][9] = 251;
    v[338][9] = 829;
    v[339][9] = 727;
    v[340][9] = 439;
    v[341][9] = 495;
    v[342][9] = 647;
    v[343][9] = 223;
    v[344][9] = 949;
    v[345][9] = 625;
    v[346][9] = 87;
    v[347][9] = 481;
    v[348][9] = 85;
    v[349][9] = 799;
    v[350][9] = 917;
    v[351][9] = 769;
    v[352][9] = 949;
    v[353][9] = 739;
    v[354][9] = 115;
    v[355][9] = 499;
    v[356][9] = 945;
    v[357][9] = 547;
    v[358][9] = 225;
    v[359][9] = 1015;
    v[360][9] = 469;
    v[361][9] = 737;
    v[362][9] = 495;
    v[363][9] = 353;
    v[364][9] = 103;
    v[365][9] = 17;
    v[366][9] = 665;
    v[367][9] = 639;
    v[368][9] = 525;
    v[369][9] = 75;
    v[370][9] = 447;
    v[371][9] = 185;
    v[372][9] = 43;
    v[373][9] = 729;
    v[374][9] = 577;
    v[375][9] = 863;
    v[376][9] = 735;
    v[377][9] = 317;
    v[378][9] = 99;
    v[379][9] = 17;
    v[380][9] = 477;
    v[381][9] = 893;
    v[382][9] = 537;
    v[383][9] = 519;
    v[384][9] = 1017;
    v[385][9] = 375;
    v[386][9] = 297;
    v[387][9] = 325;
    v[388][9] = 999;
    v[389][9] = 353;
    v[390][9] = 343;
    v[391][9] = 729;
    v[392][9] = 135;
    v[393][9] = 489;
    v[394][9] = 859;
    v[395][9] = 267;
    v[396][9] = 141;
    v[397][9] = 831;
    v[398][9] = 141;
    v[399][9] = 893;
    v[400][9] = 249;
    v[401][9] = 807;
    v[402][9] = 53;
    v[403][9] = 613;
    v[404][9] = 131;
    v[405][9] = 547;
    v[406][9] = 977;
    v[407][9] = 131;
    v[408][9] = 999;
    v[409][9] = 175;
    v[410][9] = 31;
    v[411][9] = 341;
    v[412][9] = 739;
    v[413][9] = 467;
    v[414][9] = 675;
    v[415][9] = 241;
    v[416][9] = 645;
    v[417][9] = 247;
    v[418][9] = 391;
    v[419][9] = 583;
    v[420][9] = 183;
    v[421][9] = 973;
    v[422][9] = 433;
    v[423][9] = 367;
    v[424][9] = 131;
    v[425][9] = 467;
    v[426][9] = 571;
    v[427][9] = 309;
    v[428][9] = 385;
    v[429][9] = 977;
    v[430][9] = 111;
    v[431][9] = 917;
    v[432][9] = 935;
    v[433][9] = 473;
    v[434][9] = 345;
    v[435][9] = 411;
    v[436][9] = 313;
    v[437][9] = 97;
    v[438][9] = 149;
    v[439][9] = 959;
    v[440][9] = 841;
    v[441][9] = 839;
    v[442][9] = 669;
    v[443][9] = 431;
    v[444][9] = 51;
    v[445][9] = 41;
    v[446][9] = 301;
    v[447][9] = 247;
    v[448][9] = 1015;
    v[449][9] = 377;
    v[450][9] = 329;
    v[451][9] = 945;
    v[452][9] = 269;
    v[453][9] = 67;
    v[454][9] = 979;
    v[455][9] = 581;
    v[456][9] = 643;
    v[457][9] = 823;
    v[458][9] = 557;
    v[459][9] = 91;
    v[460][9] = 405;
    v[461][9] = 117;
    v[462][9] = 801;
    v[463][9] = 509;
    v[464][9] = 347;
    v[465][9] = 893;
    v[466][9] = 303;
    v[467][9] = 227;
    v[468][9] = 783;
    v[469][9] = 555;
    v[470][9] = 867;
    v[471][9] = 99;
    v[472][9] = 703;
    v[473][9] = 111;
    v[474][9] = 797;
    v[475][9] = 873;
    v[476][9] = 541;
    v[477][9] = 919;
    v[478][9] = 513;
    v[479][9] = 343;
    v[480][9] = 319;
    v[481][9] = 517;
    v[482][9] = 135;
    v[483][9] = 871;
    v[484][9] = 917;
    v[485][9] = 285;
    v[486][9] = 663;
    v[487][9] = 301;
    v[488][9] = 15;
    v[489][9] = 763;
    v[490][9] = 89;
    v[491][9] = 323;
    v[492][9] = 757;
    v[493][9] = 317;
    v[494][9] = 807;
    v[495][9] = 309;
    v[496][9] = 1013;
    v[497][9] = 345;
    v[498][9] = 499;
    v[499][9] = 279;
    v[500][9] = 711;
    v[501][9] = 915;
    v[502][9] = 411;
    v[503][9] = 281;
    v[504][9] = 193;
    v[505][9] = 739;
    v[506][9] = 365;
    v[507][9] = 315;
    v[508][9] = 375;
    v[509][9] = 809;
    v[510][9] = 469;
    v[511][9] = 487;
    v[512][9] = 621;
    v[513][9] = 857;
    v[514][9] = 975;
    v[515][9] = 537;
    v[516][9] = 939;
    v[517][9] = 585;
    v[518][9] = 129;
    v[519][9] = 625;
    v[520][9] = 447;
    v[521][9] = 129;
    v[522][9] = 1017;
    v[523][9] = 133;
    v[524][9] = 83;
    v[525][9] = 3;
    v[526][9] = 415;
    v[527][9] = 661;
    v[528][9] = 53;
    v[529][9] = 115;
    v[530][9] = 903;
    v[531][9] = 49;
    v[532][9] = 79;
    v[533][9] = 55;
    v[534][9] = 385;
    v[535][9] = 261;
    v[536][9] = 345;
    v[537][9] = 297;
    v[538][9] = 199;
    v[539][9] = 385;
    v[540][9] = 617;
    v[541][9] = 25;
    v[542][9] = 515;
    v[543][9] = 275;
    v[544][9] = 849;
    v[545][9] = 401;
    v[546][9] = 471;
    v[547][9] = 377;
    v[548][9] = 661;
    v[549][9] = 535;
    v[550][9] = 505;
    v[551][9] = 939;
    v[552][9] = 465;
    v[553][9] = 225;
    v[554][9] = 929;
    v[555][9] = 219;
    v[556][9] = 955;
    v[557][9] = 659;
    v[558][9] = 441;
    v[559][9] = 117;
    v[560][9] = 527;
    v[561][9] = 427;
    v[562][9] = 515;
    v[563][9] = 287;
    v[564][9] = 191;
    v[565][9] = 33;
    v[566][9] = 389;
    v[567][9] = 197;
    v[568][9] = 825;
    v[569][9] = 63;
    v[570][9] = 417;
    v[571][9] = 949;
    v[572][9] = 35;
    v[573][9] = 571;
    v[574][9] = 9;
    v[575][9] = 131;
    v[576][9] = 609;
    v[577][9] = 439;
    v[578][9] = 95;
    v[579][9] = 19;
    v[580][9] = 569;
    v[581][9] = 893;
    v[582][9] = 451;
    v[583][9] = 397;
    v[584][9] = 971;
    v[585][9] = 801;
    v[586][9] = 125;
    v[587][9] = 471;
    v[588][9] = 187;
    v[589][9] = 257;
    v[590][9] = 67;
    v[591][9] = 949;
    v[592][9] = 621;
    v[593][9] = 453;
    v[594][9] = 411;
    v[595][9] = 621;
    v[596][9] = 955;
    v[597][9] = 309;
    v[598][9] = 783;
    v[599][9] = 893;
    v[600][9] = 597;
    v[601][9] = 377;
    v[602][9] = 753;
    v[603][9] = 145;
    v[604][9] = 637;
    v[605][9] = 941;
    v[606][9] = 593;
    v[607][9] = 317;
    v[608][9] = 555;
    v[609][9] = 375;
    v[610][9] = 575;
    v[611][9] = 175;
    v[612][9] = 403;
    v[613][9] = 571;
    v[614][9] = 555;
    v[615][9] = 109;
    v[616][9] = 377;
    v[617][9] = 931;
    v[618][9] = 499;
    v[619][9] = 649;
    v[620][9] = 653;
    v[621][9] = 329;
    v[622][9] = 279;
    v[623][9] = 271;
    v[624][9] = 647;
    v[625][9] = 721;
    v[626][9] = 665;
    v[627][9] = 429;
    v[628][9] = 957;
    v[629][9] = 803;
    v[630][9] = 767;
    v[631][9] = 425;
    v[632][9] = 477;
    v[633][9] = 995;
    v[634][9] = 105;
    v[635][9] = 495;
    v[636][9] = 575;
    v[637][9] = 687;
    v[638][9] = 385;
    v[639][9] = 227;
    v[640][9] = 923;
    v[641][9] = 563;
    v[642][9] = 723;
    v[643][9] = 481;
    v[644][9] = 717;
    v[645][9] = 111;
    v[646][9] = 633;
    v[647][9] = 113;
    v[648][9] = 369;
    v[649][9] = 955;
    v[650][9] = 253;
    v[651][9] = 321;
    v[652][9] = 409;
    v[653][9] = 909;
    v[654][9] = 367;
    v[655][9] = 33;
    v[656][9] = 967;
    v[657][9] = 453;
    v[658][9] = 863;
    v[659][9] = 449;
    v[660][9] = 539;
    v[661][9] = 781;
    v[662][9] = 911;
    v[663][9] = 113;
    v[664][9] = 7;
    v[665][9] = 219;
    v[666][9] = 725;
    v[667][9] = 1015;
    v[668][9] = 971;
    v[669][9] = 1021;
    v[670][9] = 525;
    v[671][9] = 785;
    v[672][9] = 873;
    v[673][9] = 191;
    v[674][9] = 893;
    v[675][9] = 297;
    v[676][9] = 507;
    v[677][9] = 215;
    v[678][9] = 21;
    v[679][9] = 153;
    v[680][9] = 645;
    v[681][9] = 913;
    v[682][9] = 755;
    v[683][9] = 371;
    v[684][9] = 881;
    v[685][9] = 113;
    v[686][9] = 903;
    v[687][9] = 225;
    v[688][9] = 49;
    v[689][9] = 587;
    v[690][9] = 201;
    v[691][9] = 927;
    v[692][9] = 429;
    v[693][9] = 599;
    v[694][9] = 513;
    v[695][9] = 97;
    v[696][9] = 319;
    v[697][9] = 331;
    v[698][9] = 833;
    v[699][9] = 325;
    v[700][9] = 887;
    v[701][9] = 139;
    v[702][9] = 927;
    v[703][9] = 399;
    v[704][9] = 163;
    v[705][9] = 307;
    v[706][9] = 803;
    v[707][9] = 169;
    v[708][9] = 1019;
    v[709][9] = 869;
    v[710][9] = 537;
    v[711][9] = 907;
    v[712][9] = 479;
    v[713][9] = 335;
    v[714][9] = 697;
    v[715][9] = 479;
    v[716][9] = 353;
    v[717][9] = 769;
    v[718][9] = 787;
    v[719][9] = 1023;
    v[720][9] = 855;
    v[721][9] = 493;
    v[722][9] = 883;
    v[723][9] = 521;
    v[724][9] = 735;
    v[725][9] = 297;
    v[726][9] = 1011;
    v[727][9] = 991;
    v[728][9] = 879;
    v[729][9] = 855;
    v[730][9] = 591;
    v[731][9] = 415;
    v[732][9] = 917;
    v[733][9] = 375;
    v[734][9] = 453;
    v[735][9] = 553;
    v[736][9] = 189;
    v[737][9] = 841;
    v[738][9] = 339;
    v[739][9] = 211;
    v[740][9] = 601;
    v[741][9] = 57;
    v[742][9] = 765;
    v[743][9] = 745;
    v[744][9] = 621;
    v[745][9] = 209;
    v[746][9] = 875;
    v[747][9] = 639;
    v[748][9] = 7;
    v[749][9] = 595;
    v[750][9] = 971;
    v[751][9] = 263;
    v[752][9] = 1009;
    v[753][9] = 201;
    v[754][9] = 23;
    v[755][9] = 77;
    v[756][9] = 621;
    v[757][9] = 33;
    v[758][9] = 535;
    v[759][9] = 963;
    v[760][9] = 661;
    v[761][9] = 523;
    v[762][9] = 263;
    v[763][9] = 917;
    v[764][9] = 103;
    v[765][9] = 623;
    v[766][9] = 231;
    v[767][9] = 47;
    v[768][9] = 301;
    v[769][9] = 549;
    v[770][9] = 337;
    v[771][9] = 675;
    v[772][9] = 189;
    v[773][9] = 357;
    v[774][9] = 1005;
    v[775][9] = 789;
    v[776][9] = 189;
    v[777][9] = 319;
    v[778][9] = 721;
    v[779][9] = 1005;
    v[780][9] = 525;
    v[781][9] = 675;
    v[782][9] = 539;
    v[783][9] = 191;
    v[784][9] = 813;
    v[785][9] = 917;
    v[786][9] = 51;
    v[787][9] = 167;
    v[788][9] = 415;
    v[789][9] = 579;
    v[790][9] = 755;
    v[791][9] = 605;
    v[792][9] = 721;
    v[793][9] = 837;
    v[794][9] = 529;
    v[795][9] = 31;
    v[796][9] = 327;
    v[797][9] = 799;
    v[798][9] = 961;
    v[799][9] = 279;
    v[800][9] = 409;
    v[801][9] = 847;
    v[802][9] = 649;
    v[803][9] = 241;
    v[804][9] = 285;
    v[805][9] = 545;
    v[806][9] = 407;
    v[807][9] = 161;
    v[808][9] = 591;
    v[809][9] = 73;
    v[810][9] = 313;
    v[811][9] = 811;
    v[812][9] = 17;
    v[813][9] = 663;
    v[814][9] = 269;
    v[815][9] = 261;
    v[816][9] = 37;
    v[817][9] = 783;
    v[818][9] = 127;
    v[819][9] = 917;
    v[820][9] = 231;
    v[821][9] = 577;
    v[822][9] = 975;
    v[823][9] = 793;
    v[824][9] = 921;
    v[825][9] = 343;
    v[826][9] = 751;
    v[827][9] = 139;
    v[828][9] = 221;
    v[829][9] = 79;
    v[830][9] = 817;
    v[831][9] = 393;
    v[832][9] = 545;
    v[833][9] = 11;
    v[834][9] = 781;
    v[835][9] = 71;
    v[836][9] = 1;
    v[837][9] = 699;
    v[838][9] = 767;
    v[839][9] = 917;
    v[840][9] = 9;
    v[841][9] = 107;
    v[842][9] = 341;
    v[843][9] = 587;
    v[844][9] = 903;
    v[845][9] = 965;
    v[846][9] = 599;
    v[847][9] = 507;
    v[848][9] = 843;
    v[849][9] = 739;
    v[850][9] = 579;
    v[851][9] = 397;
    v[852][9] = 397;
    v[853][9] = 325;
    v[854][9] = 775;
    v[855][9] = 565;
    v[856][9] = 925;
    v[857][9] = 75;
    v[858][9] = 55;
    v[859][9] = 979;
    v[860][9] = 931;
    v[861][9] = 93;
    v[862][9] = 957;
    v[863][9] = 857;
    v[864][9] = 753;
    v[865][9] = 965;
    v[866][9] = 795;
    v[867][9] = 67;
    v[868][9] = 5;
    v[869][9] = 87;
    v[870][9] = 909;
    v[871][9] = 97;
    v[872][9] = 995;
    v[873][9] = 271;
    v[874][9] = 875;
    v[875][9] = 671;
    v[876][9] = 613;
    v[877][9] = 33;
    v[878][9] = 351;
    v[879][9] = 69;
    v[880][9] = 811;
    v[881][9] = 669;
    v[882][9] = 729;
    v[883][9] = 401;
    v[884][9] = 647;
    v[885][9] = 241;
    v[886][9] = 435;
    v[887][9] = 447;
    v[888][9] = 721;
    v[889][9] = 271;
    v[890][9] = 745;
    v[891][9] = 53;
    v[892][9] = 775;
    v[893][9] = 99;
    v[894][9] = 343;
    v[895][9] = 451;
    v[896][9] = 427;
    v[897][9] = 593;
    v[898][9] = 339;
    v[899][9] = 845;
    v[900][9] = 243;
    v[901][9] = 345;
    v[902][9] = 17;
    v[903][9] = 573;
    v[904][9] = 421;
    v[905][9] = 517;
    v[906][9] = 971;
    v[907][9] = 499;
    v[908][9] = 435;
    v[909][9] = 769;
    v[910][9] = 75;
    v[911][9] = 203;
    v[912][9] = 793;
    v[913][9] = 985;
    v[914][9] = 343;
    v[915][9] = 955;
    v[916][9] = 735;
    v[917][9] = 523;
    v[918][9] = 659;
    v[919][9] = 703;
    v[920][9] = 303;
    v[921][9] = 421;
    v[922][9] = 951;
    v[923][9] = 405;
    v[924][9] = 631;
    v[925][9] = 825;
    v[926][9] = 735;
    v[927][9] = 433;
    v[928][9] = 841;
    v[929][9] = 485;
    v[930][9] = 49;
    v[931][9] = 749;
    v[932][9] = 107;
    v[933][9] = 669;
    v[934][9] = 211;
    v[935][9] = 497;
    v[936][9] = 143;
    v[937][9] = 99;
    v[938][9] = 57;
    v[939][9] = 277;
    v[940][9] = 969;
    v[941][9] = 107;
    v[942][9] = 397;
    v[943][9] = 563;
    v[944][9] = 551;
    v[945][9] = 447;
    v[946][9] = 381;
    v[947][9] = 187;
    v[948][9] = 57;
    v[949][9] = 405;
    v[950][9] = 731;
    v[951][9] = 769;
    v[952][9] = 923;
    v[953][9] = 955;
    v[954][9] = 915;
    v[955][9] = 737;
    v[956][9] = 595;
    v[957][9] = 341;
    v[958][9] = 253;
    v[959][9] = 823;
    v[960][9] = 197;
    v[961][9] = 321;
    v[962][9] = 315;
    v[963][9] = 181;
    v[964][9] = 885;
    v[965][9] = 497;
    v[966][9] = 159;
    v[967][9] = 571;
    v[968][9] = 981;
    v[969][9] = 899;
    v[970][9] = 785;
    v[971][9] = 947;
    v[972][9] = 217;
    v[973][9] = 217;
    v[974][9] = 135;
    v[975][9] = 753;
    v[976][9] = 623;
    v[977][9] = 565;
    v[978][9] = 717;
    v[979][9] = 903;
    v[980][9] = 581;
    v[981][9] = 955;
    v[982][9] = 621;
    v[983][9] = 361;
    v[984][9] = 869;
    v[985][9] = 87;
    v[986][9] = 943;
    v[987][9] = 907;
    v[988][9] = 853;
    v[989][9] = 353;
    v[990][9] = 335;
    v[991][9] = 197;
    v[992][9] = 771;
    v[993][9] = 433;
    v[994][9] = 743;
    v[995][9] = 195;
    v[996][9] = 91;
    v[997][9] = 1023;
    v[998][9] = 63;
    v[999][9] = 301;
    v[1000][9] = 647;
    v[1001][9] = 205;
    v[1002][9] = 485;
    v[1003][9] = 927;
    v[1004][9] = 1003;
    v[1005][9] = 987;
    v[1006][9] = 359;
    v[1007][9] = 577;
    v[1008][9] = 147;
    v[1009][9] = 141;
    v[1010][9] = 1017;
    v[1011][9] = 701;
    v[1012][9] = 273;
    v[1013][9] = 89;
    v[1014][9] = 589;
    v[1015][9] = 487;
    v[1016][9] = 859;
    v[1017][9] = 343;
    v[1018][9] = 91;
    v[1019][9] = 847;
    v[1020][9] = 341;
    v[1021][9] = 173;
    v[1022][9] = 287;
    v[1023][9] = 1003;
    v[1024][9] = 289;
    v[1025][9] = 639;
    v[1026][9] = 983;
    v[1027][9] = 685;
    v[1028][9] = 697;
    v[1029][9] = 35;
    v[1030][9] = 701;
    v[1031][9] = 645;
    v[1032][9] = 911;
    v[1033][9] = 501;
    v[1034][9] = 705;
    v[1035][9] = 873;
    v[1036][9] = 763;
    v[1037][9] = 745;
    v[1038][9] = 657;
    v[1039][9] = 559;
    v[1040][9] = 699;
    v[1041][9] = 315;
    v[1042][9] = 347;
    v[1043][9] = 429;
    v[1044][9] = 197;
    v[1045][9] = 165;
    v[1046][9] = 955;
    v[1047][9] = 859;
    v[1048][9] = 167;
    v[1049][9] = 303;
    v[1050][9] = 833;
    v[1051][9] = 531;
    v[1052][9] = 473;
    v[1053][9] = 635;
    v[1054][9] = 641;
    v[1055][9] = 195;
    v[1056][9] = 589;
    v[1057][9] = 821;
    v[1058][9] = 205;
    v[1059][9] = 3;
    v[1060][9] = 635;
    v[1061][9] = 371;
    v[1062][9] = 891;
    v[1063][9] = 249;
    v[1064][9] = 123;
    v[1065][9] = 77;
    v[1066][9] = 623;
    v[1067][9] = 993;
    v[1068][9] = 401;
    v[1069][9] = 525;
    v[1070][9] = 427;
    v[1071][9] = 71;
    v[1072][9] = 655;
    v[1073][9] = 951;
    v[1074][9] = 357;
    v[1075][9] = 851;
    v[1076][9] = 899;
    v[1077][9] = 535;
    v[1078][9] = 493;
    v[1079][9] = 323;
    v[1080][9] = 1003;
    v[1081][9] = 343;
    v[1082][9] = 515;
    v[1083][9] = 859;
    v[1084][9] = 1017;
    v[1085][9] = 5;
    v[1086][9] = 423;
    v[1087][9] = 315;
    v[1088][9] = 1011;
    v[1089][9] = 703;
    v[1090][9] = 41;
    v[1091][9] = 777;
    v[1092][9] = 163;
    v[1093][9] = 95;
    v[1094][9] = 831;
    v[1095][9] = 79;
    v[1096][9] = 975;
    v[1097][9] = 235;
    v[1098][9] = 633;
    v[1099][9] = 723;
    v[1100][9] = 297;
    v[1101][9] = 589;
    v[1102][9] = 317;
    v[1103][9] = 679;
    v[1104][9] = 981;
    v[1105][9] = 195;
    v[1106][9] = 399;
    v[1107][9] = 1003;
    v[1108][9] = 121;
    v[1109][9] = 501;
    v[1110][9] = 155;

    v[161][10] = 7;
    v[162][10] = 2011;
    v[163][10] = 1001;
    v[164][10] = 49;
    v[165][10] = 825;
    v[166][10] = 415;
    v[167][10] = 1441;
    v[168][10] = 383;
    v[169][10] = 1581;
    v[170][10] = 623;
    v[171][10] = 1621;
    v[172][10] = 1319;
    v[173][10] = 1387;
    v[174][10] = 619;
    v[175][10] = 839;
    v[176][10] = 217;
    v[177][10] = 75;
    v[178][10] = 1955;
    v[179][10] = 505;
    v[180][10] = 281;
    v[181][10] = 1629;
    v[182][10] = 1379;
    v[183][10] = 53;
    v[184][10] = 1111;
    v[185][10] = 1399;
    v[186][10] = 301;
    v[187][10] = 209;
    v[188][10] = 49;
    v[189][10] = 155;
    v[190][10] = 1647;
    v[191][10] = 631;
    v[192][10] = 129;
    v[193][10] = 1569;
    v[194][10] = 335;
    v[195][10] = 67;
    v[196][10] = 1955;
    v[197][10] = 1611;
    v[198][10] = 2021;
    v[199][10] = 1305;
    v[200][10] = 121;
    v[201][10] = 37;
    v[202][10] = 877;
    v[203][10] = 835;
    v[204][10] = 1457;
    v[205][10] = 669;
    v[206][10] = 1405;
    v[207][10] = 935;
    v[208][10] = 1735;
    v[209][10] = 665;
    v[210][10] = 551;
    v[211][10] = 789;
    v[212][10] = 1543;
    v[213][10] = 1267;
    v[214][10] = 1027;
    v[215][10] = 1;
    v[216][10] = 1911;
    v[217][10] = 163;
    v[218][10] = 1929;
    v[219][10] = 67;
    v[220][10] = 1975;
    v[221][10] = 1681;
    v[222][10] = 1413;
    v[223][10] = 191;
    v[224][10] = 1711;
    v[225][10] = 1307;
    v[226][10] = 401;
    v[227][10] = 725;
    v[228][10] = 1229;
    v[229][10] = 1403;
    v[230][10] = 1609;
    v[231][10] = 2035;
    v[232][10] = 917;
    v[233][10] = 921;
    v[234][10] = 1789;
    v[235][10] = 41;
    v[236][10] = 2003;
    v[237][10] = 187;
    v[238][10] = 67;
    v[239][10] = 1635;
    v[240][10] = 717;
    v[241][10] = 1449;
    v[242][10] = 277;
    v[243][10] = 1903;
    v[244][10] = 1179;
    v[245][10] = 363;
    v[246][10] = 1211;
    v[247][10] = 1231;
    v[248][10] = 647;
    v[249][10] = 1261;
    v[250][10] = 1029;
    v[251][10] = 1485;
    v[252][10] = 1309;
    v[253][10] = 1149;
    v[254][10] = 317;
    v[255][10] = 1335;
    v[256][10] = 171;
    v[257][10] = 243;
    v[258][10] = 271;
    v[259][10] = 1055;
    v[260][10] = 1601;
    v[261][10] = 1129;
    v[262][10] = 1653;
    v[263][10] = 205;
    v[264][10] = 1463;
    v[265][10] = 1681;
    v[266][10] = 1621;
    v[267][10] = 197;
    v[268][10] = 951;
    v[269][10] = 573;
    v[270][10] = 1697;
    v[271][10] = 1265;
    v[272][10] = 1321;
    v[273][10] = 1805;
    v[274][10] = 1235;
    v[275][10] = 1853;
    v[276][10] = 1307;
    v[277][10] = 945;
    v[278][10] = 1197;
    v[279][10] = 1411;
    v[280][10] = 833;
    v[281][10] = 273;
    v[282][10] = 1517;
    v[283][10] = 1747;
    v[284][10] = 1095;
    v[285][10] = 1345;
    v[286][10] = 869;
    v[287][10] = 57;
    v[288][10] = 1383;
    v[289][10] = 221;
    v[290][10] = 1713;
    v[291][10] = 335;
    v[292][10] = 1751;
    v[293][10] = 1141;
    v[294][10] = 839;
    v[295][10] = 523;
    v[296][10] = 1861;
    v[297][10] = 1105;
    v[298][10] = 389;
    v[299][10] = 1177;
    v[300][10] = 1877;
    v[301][10] = 805;
    v[302][10] = 93;
    v[303][10] = 1591;
    v[304][10] = 423;
    v[305][10] = 1835;
    v[306][10] = 99;
    v[307][10] = 1781;
    v[308][10] = 1515;
    v[309][10] = 1909;
    v[310][10] = 1011;
    v[311][10] = 303;
    v[312][10] = 385;
    v[313][10] = 1635;
    v[314][10] = 357;
    v[315][10] = 973;
    v[316][10] = 1781;
    v[317][10] = 1707;
    v[318][10] = 1363;
    v[319][10] = 1053;
    v[320][10] = 649;
    v[321][10] = 1469;
    v[322][10] = 623;
    v[323][10] = 1429;
    v[324][10] = 1241;
    v[325][10] = 1151;
    v[326][10] = 1055;
    v[327][10] = 503;
    v[328][10] = 921;
    v[329][10] = 3;
    v[330][10] = 349;
    v[331][10] = 1149;
    v[332][10] = 293;
    v[333][10] = 45;
    v[334][10] = 303;
    v[335][10] = 877;
    v[336][10] = 1565;
    v[337][10] = 1583;
    v[338][10] = 1001;
    v[339][10] = 663;
    v[340][10] = 1535;
    v[341][10] = 395;
    v[342][10] = 1141;
    v[343][10] = 1481;
    v[344][10] = 1797;
    v[345][10] = 643;
    v[346][10] = 1507;
    v[347][10] = 465;
    v[348][10] = 2027;
    v[349][10] = 1695;
    v[350][10] = 367;
    v[351][10] = 937;
    v[352][10] = 719;
    v[353][10] = 545;
    v[354][10] = 1991;
    v[355][10] = 83;
    v[356][10] = 819;
    v[357][10] = 239;
    v[358][10] = 1791;
    v[359][10] = 1461;
    v[360][10] = 1647;
    v[361][10] = 1501;
    v[362][10] = 1161;
    v[363][10] = 1629;
    v[364][10] = 139;
    v[365][10] = 1595;
    v[366][10] = 1921;
    v[367][10] = 1267;
    v[368][10] = 1415;
    v[369][10] = 509;
    v[370][10] = 347;
    v[371][10] = 777;
    v[372][10] = 1083;
    v[373][10] = 363;
    v[374][10] = 269;
    v[375][10] = 1015;
    v[376][10] = 1809;
    v[377][10] = 1105;
    v[378][10] = 1429;
    v[379][10] = 1471;
    v[380][10] = 2019;
    v[381][10] = 381;
    v[382][10] = 2025;
    v[383][10] = 1223;
    v[384][10] = 827;
    v[385][10] = 1733;
    v[386][10] = 887;
    v[387][10] = 1321;
    v[388][10] = 803;
    v[389][10] = 1951;
    v[390][10] = 1297;
    v[391][10] = 1995;
    v[392][10] = 833;
    v[393][10] = 1107;
    v[394][10] = 1135;
    v[395][10] = 1181;
    v[396][10] = 1251;
    v[397][10] = 983;
    v[398][10] = 1389;
    v[399][10] = 1565;
    v[400][10] = 273;
    v[401][10] = 137;
    v[402][10] = 71;
    v[403][10] = 735;
    v[404][10] = 1005;
    v[405][10] = 933;
    v[406][10] = 67;
    v[407][10] = 1471;
    v[408][10] = 551;
    v[409][10] = 457;
    v[410][10] = 1667;
    v[411][10] = 1729;
    v[412][10] = 919;
    v[413][10] = 285;
    v[414][10] = 1629;
    v[415][10] = 1815;
    v[416][10] = 653;
    v[417][10] = 1919;
    v[418][10] = 1039;
    v[419][10] = 531;
    v[420][10] = 393;
    v[421][10] = 1411;
    v[422][10] = 359;
    v[423][10] = 221;
    v[424][10] = 699;
    v[425][10] = 1485;
    v[426][10] = 471;
    v[427][10] = 1357;
    v[428][10] = 1715;
    v[429][10] = 595;
    v[430][10] = 1677;
    v[431][10] = 153;
    v[432][10] = 1903;
    v[433][10] = 1281;
    v[434][10] = 215;
    v[435][10] = 781;
    v[436][10] = 543;
    v[437][10] = 293;
    v[438][10] = 1807;
    v[439][10] = 965;
    v[440][10] = 1695;
    v[441][10] = 443;
    v[442][10] = 1985;
    v[443][10] = 321;
    v[444][10] = 879;
    v[445][10] = 1227;
    v[446][10] = 1915;
    v[447][10] = 839;
    v[448][10] = 1945;
    v[449][10] = 1993;
    v[450][10] = 1165;
    v[451][10] = 51;
    v[452][10] = 557;
    v[453][10] = 723;
    v[454][10] = 1491;
    v[455][10] = 817;
    v[456][10] = 1237;
    v[457][10] = 947;
    v[458][10] = 1215;
    v[459][10] = 1911;
    v[460][10] = 1225;
    v[461][10] = 1965;
    v[462][10] = 1889;
    v[463][10] = 1503;
    v[464][10] = 1177;
    v[465][10] = 73;
    v[466][10] = 1767;
    v[467][10] = 303;
    v[468][10] = 177;
    v[469][10] = 1897;
    v[470][10] = 1401;
    v[471][10] = 321;
    v[472][10] = 921;
    v[473][10] = 217;
    v[474][10] = 1779;
    v[475][10] = 327;
    v[476][10] = 1889;
    v[477][10] = 333;
    v[478][10] = 615;
    v[479][10] = 1665;
    v[480][10] = 1825;
    v[481][10] = 1639;
    v[482][10] = 237;
    v[483][10] = 1205;
    v[484][10] = 361;
    v[485][10] = 129;
    v[486][10] = 1655;
    v[487][10] = 983;
    v[488][10] = 1089;
    v[489][10] = 1171;
    v[490][10] = 401;
    v[491][10] = 677;
    v[492][10] = 643;
    v[493][10] = 749;
    v[494][10] = 303;
    v[495][10] = 1407;
    v[496][10] = 1873;
    v[497][10] = 1579;
    v[498][10] = 1491;
    v[499][10] = 1393;
    v[500][10] = 1247;
    v[501][10] = 789;
    v[502][10] = 763;
    v[503][10] = 49;
    v[504][10] = 5;
    v[505][10] = 1607;
    v[506][10] = 1891;
    v[507][10] = 735;
    v[508][10] = 1557;
    v[509][10] = 1909;
    v[510][10] = 1765;
    v[511][10] = 1777;
    v[512][10] = 1127;
    v[513][10] = 813;
    v[514][10] = 695;
    v[515][10] = 97;
    v[516][10] = 731;
    v[517][10] = 1503;
    v[518][10] = 1751;
    v[519][10] = 333;
    v[520][10] = 769;
    v[521][10] = 865;
    v[522][10] = 693;
    v[523][10] = 377;
    v[524][10] = 1919;
    v[525][10] = 957;
    v[526][10] = 1359;
    v[527][10] = 1627;
    v[528][10] = 1039;
    v[529][10] = 1783;
    v[530][10] = 1065;
    v[531][10] = 1665;
    v[532][10] = 1917;
    v[533][10] = 1947;
    v[534][10] = 991;
    v[535][10] = 1997;
    v[536][10] = 841;
    v[537][10] = 459;
    v[538][10] = 221;
    v[539][10] = 327;
    v[540][10] = 1595;
    v[541][10] = 1881;
    v[542][10] = 1269;
    v[543][10] = 1007;
    v[544][10] = 129;
    v[545][10] = 1413;
    v[546][10] = 475;
    v[547][10] = 1105;
    v[548][10] = 791;
    v[549][10] = 1983;
    v[550][10] = 1359;
    v[551][10] = 503;
    v[552][10] = 691;
    v[553][10] = 659;
    v[554][10] = 691;
    v[555][10] = 343;
    v[556][10] = 1375;
    v[557][10] = 1919;
    v[558][10] = 263;
    v[559][10] = 1373;
    v[560][10] = 603;
    v[561][10] = 1383;
    v[562][10] = 297;
    v[563][10] = 781;
    v[564][10] = 145;
    v[565][10] = 285;
    v[566][10] = 767;
    v[567][10] = 1739;
    v[568][10] = 1715;
    v[569][10] = 715;
    v[570][10] = 317;
    v[571][10] = 1333;
    v[572][10] = 85;
    v[573][10] = 831;
    v[574][10] = 1615;
    v[575][10] = 81;
    v[576][10] = 1667;
    v[577][10] = 1467;
    v[578][10] = 1457;
    v[579][10] = 1453;
    v[580][10] = 1825;
    v[581][10] = 109;
    v[582][10] = 387;
    v[583][10] = 1207;
    v[584][10] = 2039;
    v[585][10] = 213;
    v[586][10] = 1351;
    v[587][10] = 1329;
    v[588][10] = 1173;
    v[589][10] = 57;
    v[590][10] = 1769;
    v[591][10] = 951;
    v[592][10] = 183;
    v[593][10] = 23;
    v[594][10] = 451;
    v[595][10] = 1155;
    v[596][10] = 1551;
    v[597][10] = 2037;
    v[598][10] = 811;
    v[599][10] = 635;
    v[600][10] = 1671;
    v[601][10] = 1451;
    v[602][10] = 863;
    v[603][10] = 1499;
    v[604][10] = 1673;
    v[605][10] = 363;
    v[606][10] = 1029;
    v[607][10] = 1077;
    v[608][10] = 1525;
    v[609][10] = 277;
    v[610][10] = 1023;
    v[611][10] = 655;
    v[612][10] = 665;
    v[613][10] = 1869;
    v[614][10] = 1255;
    v[615][10] = 965;
    v[616][10] = 277;
    v[617][10] = 1601;
    v[618][10] = 329;
    v[619][10] = 1603;
    v[620][10] = 1901;
    v[621][10] = 395;
    v[622][10] = 65;
    v[623][10] = 1307;
    v[624][10] = 2029;
    v[625][10] = 21;
    v[626][10] = 1321;
    v[627][10] = 543;
    v[628][10] = 1569;
    v[629][10] = 1185;
    v[630][10] = 1905;
    v[631][10] = 1701;
    v[632][10] = 413;
    v[633][10] = 2041;
    v[634][10] = 1697;
    v[635][10] = 725;
    v[636][10] = 1417;
    v[637][10] = 1847;
    v[638][10] = 411;
    v[639][10] = 211;
    v[640][10] = 915;
    v[641][10] = 1891;
    v[642][10] = 17;
    v[643][10] = 1877;
    v[644][10] = 1699;
    v[645][10] = 687;
    v[646][10] = 1089;
    v[647][10] = 1973;
    v[648][10] = 1809;
    v[649][10] = 851;
    v[650][10] = 1495;
    v[651][10] = 1257;
    v[652][10] = 63;
    v[653][10] = 1323;
    v[654][10] = 1307;
    v[655][10] = 609;
    v[656][10] = 881;
    v[657][10] = 1543;
    v[658][10] = 177;
    v[659][10] = 617;
    v[660][10] = 1505;
    v[661][10] = 1747;
    v[662][10] = 1537;
    v[663][10] = 925;
    v[664][10] = 183;
    v[665][10] = 77;
    v[666][10] = 1723;
    v[667][10] = 1877;
    v[668][10] = 1703;
    v[669][10] = 397;
    v[670][10] = 459;
    v[671][10] = 521;
    v[672][10] = 257;
    v[673][10] = 1177;
    v[674][10] = 389;
    v[675][10] = 1947;
    v[676][10] = 1553;
    v[677][10] = 1583;
    v[678][10] = 1831;
    v[679][10] = 261;
    v[680][10] = 485;
    v[681][10] = 289;
    v[682][10] = 1281;
    v[683][10] = 1543;
    v[684][10] = 1591;
    v[685][10] = 1123;
    v[686][10] = 573;
    v[687][10] = 821;
    v[688][10] = 1065;
    v[689][10] = 1933;
    v[690][10] = 1373;
    v[691][10] = 2005;
    v[692][10] = 905;
    v[693][10] = 207;
    v[694][10] = 173;
    v[695][10] = 1573;
    v[696][10] = 1597;
    v[697][10] = 573;
    v[698][10] = 1883;
    v[699][10] = 1795;
    v[700][10] = 1499;
    v[701][10] = 1743;
    v[702][10] = 553;
    v[703][10] = 335;
    v[704][10] = 333;
    v[705][10] = 1645;
    v[706][10] = 791;
    v[707][10] = 871;
    v[708][10] = 1157;
    v[709][10] = 969;
    v[710][10] = 557;
    v[711][10] = 141;
    v[712][10] = 223;
    v[713][10] = 1129;
    v[714][10] = 1685;
    v[715][10] = 423;
    v[716][10] = 1069;
    v[717][10] = 391;
    v[718][10] = 99;
    v[719][10] = 95;
    v[720][10] = 1847;
    v[721][10] = 531;
    v[722][10] = 1859;
    v[723][10] = 1833;
    v[724][10] = 1833;
    v[725][10] = 341;
    v[726][10] = 237;
    v[727][10] = 1997;
    v[728][10] = 1799;
    v[729][10] = 409;
    v[730][10] = 431;
    v[731][10] = 1917;
    v[732][10] = 363;
    v[733][10] = 335;
    v[734][10] = 1039;
    v[735][10] = 1085;
    v[736][10] = 1657;
    v[737][10] = 1975;
    v[738][10] = 1527;
    v[739][10] = 1111;
    v[740][10] = 659;
    v[741][10] = 389;
    v[742][10] = 899;
    v[743][10] = 595;
    v[744][10] = 1439;
    v[745][10] = 1861;
    v[746][10] = 1979;
    v[747][10] = 1569;
    v[748][10] = 1087;
    v[749][10] = 1009;
    v[750][10] = 165;
    v[751][10] = 1895;
    v[752][10] = 1481;
    v[753][10] = 1583;
    v[754][10] = 29;
    v[755][10] = 1193;
    v[756][10] = 1673;
    v[757][10] = 1075;
    v[758][10] = 301;
    v[759][10] = 1081;
    v[760][10] = 1377;
    v[761][10] = 1747;
    v[762][10] = 1497;
    v[763][10] = 1103;
    v[764][10] = 1789;
    v[765][10] = 887;
    v[766][10] = 739;
    v[767][10] = 1577;
    v[768][10] = 313;
    v[769][10] = 1367;
    v[770][10] = 1299;
    v[771][10] = 1801;
    v[772][10] = 1131;
    v[773][10] = 1837;
    v[774][10] = 73;
    v[775][10] = 1865;
    v[776][10] = 1065;
    v[777][10] = 843;
    v[778][10] = 635;
    v[779][10] = 55;
    v[780][10] = 1655;
    v[781][10] = 913;
    v[782][10] = 1037;
    v[783][10] = 223;
    v[784][10] = 1871;
    v[785][10] = 1161;
    v[786][10] = 461;
    v[787][10] = 479;
    v[788][10] = 511;
    v[789][10] = 1721;
    v[790][10] = 1107;
    v[791][10] = 389;
    v[792][10] = 151;
    v[793][10] = 35;
    v[794][10] = 375;
    v[795][10] = 1099;
    v[796][10] = 937;
    v[797][10] = 1185;
    v[798][10] = 1701;
    v[799][10] = 769;
    v[800][10] = 639;
    v[801][10] = 1633;
    v[802][10] = 1609;
    v[803][10] = 379;
    v[804][10] = 1613;
    v[805][10] = 2031;
    v[806][10] = 685;
    v[807][10] = 289;
    v[808][10] = 975;
    v[809][10] = 671;
    v[810][10] = 1599;
    v[811][10] = 1447;
    v[812][10] = 871;
    v[813][10] = 647;
    v[814][10] = 99;
    v[815][10] = 139;
    v[816][10] = 1427;
    v[817][10] = 959;
    v[818][10] = 89;
    v[819][10] = 117;
    v[820][10] = 841;
    v[821][10] = 891;
    v[822][10] = 1959;
    v[823][10] = 223;
    v[824][10] = 1697;
    v[825][10] = 1145;
    v[826][10] = 499;
    v[827][10] = 1435;
    v[828][10] = 1809;
    v[829][10] = 1413;
    v[830][10] = 1445;
    v[831][10] = 1675;
    v[832][10] = 171;
    v[833][10] = 1073;
    v[834][10] = 1349;
    v[835][10] = 1545;
    v[836][10] = 2039;
    v[837][10] = 1027;
    v[838][10] = 1563;
    v[839][10] = 859;
    v[840][10] = 215;
    v[841][10] = 1673;
    v[842][10] = 1919;
    v[843][10] = 1633;
    v[844][10] = 779;
    v[845][10] = 411;
    v[846][10] = 1845;
    v[847][10] = 1477;
    v[848][10] = 1489;
    v[849][10] = 447;
    v[850][10] = 1545;
    v[851][10] = 351;
    v[852][10] = 1989;
    v[853][10] = 495;
    v[854][10] = 183;
    v[855][10] = 1639;
    v[856][10] = 1385;
    v[857][10] = 1805;
    v[858][10] = 1097;
    v[859][10] = 1249;
    v[860][10] = 1431;
    v[861][10] = 1571;
    v[862][10] = 591;
    v[863][10] = 697;
    v[864][10] = 1509;
    v[865][10] = 709;
    v[866][10] = 31;
    v[867][10] = 1563;
    v[868][10] = 165;
    v[869][10] = 513;
    v[870][10] = 1425;
    v[871][10] = 1299;
    v[872][10] = 1081;
    v[873][10] = 145;
    v[874][10] = 1841;
    v[875][10] = 1211;
    v[876][10] = 941;
    v[877][10] = 609;
    v[878][10] = 845;
    v[879][10] = 1169;
    v[880][10] = 1865;
    v[881][10] = 1593;
    v[882][10] = 347;
    v[883][10] = 293;
    v[884][10] = 1277;
    v[885][10] = 157;
    v[886][10] = 211;
    v[887][10] = 93;
    v[888][10] = 1679;
    v[889][10] = 1799;
    v[890][10] = 527;
    v[891][10] = 41;
    v[892][10] = 473;
    v[893][10] = 563;
    v[894][10] = 187;
    v[895][10] = 1525;
    v[896][10] = 575;
    v[897][10] = 1579;
    v[898][10] = 857;
    v[899][10] = 703;
    v[900][10] = 1211;
    v[901][10] = 647;
    v[902][10] = 709;
    v[903][10] = 981;
    v[904][10] = 285;
    v[905][10] = 697;
    v[906][10] = 163;
    v[907][10] = 981;
    v[908][10] = 153;
    v[909][10] = 1515;
    v[910][10] = 47;
    v[911][10] = 1553;
    v[912][10] = 599;
    v[913][10] = 225;
    v[914][10] = 1147;
    v[915][10] = 381;
    v[916][10] = 135;
    v[917][10] = 821;
    v[918][10] = 1965;
    v[919][10] = 609;
    v[920][10] = 1033;
    v[921][10] = 983;
    v[922][10] = 503;
    v[923][10] = 1117;
    v[924][10] = 327;
    v[925][10] = 453;
    v[926][10] = 2005;
    v[927][10] = 1257;
    v[928][10] = 343;
    v[929][10] = 1649;
    v[930][10] = 1199;
    v[931][10] = 599;
    v[932][10] = 1877;
    v[933][10] = 569;
    v[934][10] = 695;
    v[935][10] = 1587;
    v[936][10] = 1475;
    v[937][10] = 187;
    v[938][10] = 973;
    v[939][10] = 233;
    v[940][10] = 511;
    v[941][10] = 51;
    v[942][10] = 1083;
    v[943][10] = 665;
    v[944][10] = 1321;
    v[945][10] = 531;
    v[946][10] = 1875;
    v[947][10] = 1939;
    v[948][10] = 859;
    v[949][10] = 1507;
    v[950][10] = 1979;
    v[951][10] = 1203;
    v[952][10] = 1965;
    v[953][10] = 737;
    v[954][10] = 921;
    v[955][10] = 1565;
    v[956][10] = 1943;
    v[957][10] = 819;
    v[958][10] = 223;
    v[959][10] = 365;
    v[960][10] = 167;
    v[961][10] = 1705;
    v[962][10] = 413;
    v[963][10] = 1577;
    v[964][10] = 745;
    v[965][10] = 1573;
    v[966][10] = 655;
    v[967][10] = 1633;
    v[968][10] = 1003;
    v[969][10] = 91;
    v[970][10] = 1123;
    v[971][10] = 477;
    v[972][10] = 1741;
    v[973][10] = 1663;
    v[974][10] = 35;
    v[975][10] = 715;
    v[976][10] = 37;
    v[977][10] = 1513;
    v[978][10] = 815;
    v[979][10] = 941;
    v[980][10] = 1379;
    v[981][10] = 263;
    v[982][10] = 1831;
    v[983][10] = 1735;
    v[984][10] = 1111;
    v[985][10] = 1449;
    v[986][10] = 353;
    v[987][10] = 1941;
    v[988][10] = 1655;
    v[989][10] = 1349;
    v[990][10] = 877;
    v[991][10] = 285;
    v[992][10] = 1723;
    v[993][10] = 125;
    v[994][10] = 1753;
    v[995][10] = 985;
    v[996][10] = 723;
    v[997][10] = 175;
    v[998][10] = 439;
    v[999][10] = 791;
    v[1000][10] = 1051;
    v[1001][10] = 1261;
    v[1002][10] = 717;
    v[1003][10] = 1555;
    v[1004][10] = 1757;
    v[1005][10] = 1777;
    v[1006][10] = 577;
    v[1007][10] = 1583;
    v[1008][10] = 1957;
    v[1009][10] = 873;
    v[1010][10] = 331;
    v[1011][10] = 1163;
    v[1012][10] = 313;
    v[1013][10] = 1;
    v[1014][10] = 1963;
    v[1015][10] = 963;
    v[1016][10] = 1905;
    v[1017][10] = 821;
    v[1018][10] = 1677;
    v[1019][10] = 185;
    v[1020][10] = 709;
    v[1021][10] = 545;
    v[1022][10] = 1723;
    v[1023][10] = 215;
    v[1024][10] = 1885;
    v[1025][10] = 1249;
    v[1026][10] = 583;
    v[1027][10] = 1803;
    v[1028][10] = 839;
    v[1029][10] = 885;
    v[1030][10] = 485;
    v[1031][10] = 413;
    v[1032][10] = 1767;
    v[1033][10] = 425;
    v[1034][10] = 129;
    v[1035][10] = 1035;
    v[1036][10] = 329;
    v[1037][10] = 1263;
    v[1038][10] = 1881;
    v[1039][10] = 1779;
    v[1040][10] = 1565;
    v[1041][10] = 359;
    v[1042][10] = 367;
    v[1043][10] = 453;
    v[1044][10] = 707;
    v[1045][10] = 1419;
    v[1046][10] = 831;
    v[1047][10] = 1889;
    v[1048][10] = 887;
    v[1049][10] = 1871;
    v[1050][10] = 1869;
    v[1051][10] = 747;
    v[1052][10] = 223;
    v[1053][10] = 1547;
    v[1054][10] = 1799;
    v[1055][10] = 433;
    v[1056][10] = 1441;
    v[1057][10] = 553;
    v[1058][10] = 2021;
    v[1059][10] = 1303;
    v[1060][10] = 1505;
    v[1061][10] = 1735;
    v[1062][10] = 1619;
    v[1063][10] = 1065;
    v[1064][10] = 1161;
    v[1065][10] = 2047;
    v[1066][10] = 347;
    v[1067][10] = 867;
    v[1068][10] = 881;
    v[1069][10] = 1447;
    v[1070][10] = 329;
    v[1071][10] = 781;
    v[1072][10] = 1065;
    v[1073][10] = 219;
    v[1074][10] = 589;
    v[1075][10] = 645;
    v[1076][10] = 1257;
    v[1077][10] = 1833;
    v[1078][10] = 749;
    v[1079][10] = 1841;
    v[1080][10] = 1733;
    v[1081][10] = 1179;
    v[1082][10] = 1191;
    v[1083][10] = 1025;
    v[1084][10] = 1639;
    v[1085][10] = 1955;
    v[1086][10] = 1423;
    v[1087][10] = 1685;
    v[1088][10] = 1711;
    v[1089][10] = 493;
    v[1090][10] = 549;
    v[1091][10] = 783;
    v[1092][10] = 1653;
    v[1093][10] = 397;
    v[1094][10] = 895;
    v[1095][10] = 233;
    v[1096][10] = 759;
    v[1097][10] = 1505;
    v[1098][10] = 677;
    v[1099][10] = 1449;
    v[1100][10] = 1573;
    v[1101][10] = 1297;
    v[1102][10] = 1821;
    v[1103][10] = 1691;
    v[1104][10] = 791;
    v[1105][10] = 289;
    v[1106][10] = 1187;
    v[1107][10] = 867;
    v[1108][10] = 1535;
    v[1109][10] = 575;
    v[1110][10] = 183;

    v[337][11] = 3915;
    v[338][11] = 97;
    v[339][11] = 3047;
    v[340][11] = 937;
    v[341][11] = 2897;
    v[342][11] = 953;
    v[343][11] = 127;
    v[344][11] = 1201;
    v[345][11] = 3819;
    v[346][11] = 193;
    v[347][11] = 2053;
    v[348][11] = 3061;
    v[349][11] = 3759;
    v[350][11] = 1553;
    v[351][11] = 2007;
    v[352][11] = 2493;
    v[353][11] = 603;
    v[354][11] = 3343;
    v[355][11] = 3751;
    v[356][11] = 1059;
    v[357][11] = 783;
    v[358][11] = 1789;
    v[359][11] = 1589;
    v[360][11] = 283;
    v[361][11] = 1093;
    v[362][11] = 3919;
    v[363][11] = 2747;
    v[364][11] = 277;
    v[365][11] = 2605;
    v[366][11] = 2169;
    v[367][11] = 2905;
    v[368][11] = 721;
    v[369][11] = 4069;
    v[370][11] = 233;
    v[371][11] = 261;
    v[372][11] = 1137;
    v[373][11] = 3993;
    v[374][11] = 3619;
    v[375][11] = 2881;
    v[376][11] = 1275;
    v[377][11] = 3865;
    v[378][11] = 1299;
    v[379][11] = 3757;
    v[380][11] = 1193;
    v[381][11] = 733;
    v[382][11] = 993;
    v[383][11] = 1153;
    v[384][11] = 2945;
    v[385][11] = 3163;
    v[386][11] = 3179;
    v[387][11] = 437;
    v[388][11] = 271;
    v[389][11] = 3493;
    v[390][11] = 3971;
    v[391][11] = 1005;
    v[392][11] = 2615;
    v[393][11] = 2253;
    v[394][11] = 1131;
    v[395][11] = 585;
    v[396][11] = 2775;
    v[397][11] = 2171;
    v[398][11] = 2383;
    v[399][11] = 2937;
    v[400][11] = 2447;
    v[401][11] = 1745;
    v[402][11] = 663;
    v[403][11] = 1515;
    v[404][11] = 3767;
    v[405][11] = 2709;
    v[406][11] = 1767;
    v[407][11] = 3185;
    v[408][11] = 3017;
    v[409][11] = 2815;
    v[410][11] = 1829;
    v[411][11] = 87;
    v[412][11] = 3341;
    v[413][11] = 793;
    v[414][11] = 2627;
    v[415][11] = 2169;
    v[416][11] = 1875;
    v[417][11] = 3745;
    v[418][11] = 367;
    v[419][11] = 3783;
    v[420][11] = 783;
    v[421][11] = 827;
    v[422][11] = 3253;
    v[423][11] = 2639;
    v[424][11] = 2955;
    v[425][11] = 3539;
    v[426][11] = 1579;
    v[427][11] = 2109;
    v[428][11] = 379;
    v[429][11] = 2939;
    v[430][11] = 3019;
    v[431][11] = 1999;
    v[432][11] = 2253;
    v[433][11] = 2911;
    v[434][11] = 3733;
    v[435][11] = 481;
    v[436][11] = 1767;
    v[437][11] = 1055;
    v[438][11] = 4019;
    v[439][11] = 4085;
    v[440][11] = 105;
    v[441][11] = 1829;
    v[442][11] = 2097;
    v[443][11] = 2379;
    v[444][11] = 1567;
    v[445][11] = 2713;
    v[446][11] = 737;
    v[447][11] = 3423;
    v[448][11] = 3941;
    v[449][11] = 2659;
    v[450][11] = 3961;
    v[451][11] = 1755;
    v[452][11] = 3613;
    v[453][11] = 1937;
    v[454][11] = 1559;
    v[455][11] = 2287;
    v[456][11] = 2743;
    v[457][11] = 67;
    v[458][11] = 2859;
    v[459][11] = 325;
    v[460][11] = 2601;
    v[461][11] = 1149;
    v[462][11] = 3259;
    v[463][11] = 2403;
    v[464][11] = 3947;
    v[465][11] = 2011;
    v[466][11] = 175;
    v[467][11] = 3389;
    v[468][11] = 3915;
    v[469][11] = 1315;
    v[470][11] = 2447;
    v[471][11] = 141;
    v[472][11] = 359;
    v[473][11] = 3609;
    v[474][11] = 3933;
    v[475][11] = 729;
    v[476][11] = 2051;
    v[477][11] = 1755;
    v[478][11] = 2149;
    v[479][11] = 2107;
    v[480][11] = 1741;
    v[481][11] = 1051;
    v[482][11] = 3681;
    v[483][11] = 471;
    v[484][11] = 1055;
    v[485][11] = 845;
    v[486][11] = 257;
    v[487][11] = 1559;
    v[488][11] = 1061;
    v[489][11] = 2803;
    v[490][11] = 2219;
    v[491][11] = 1315;
    v[492][11] = 1369;
    v[493][11] = 3211;
    v[494][11] = 4027;
    v[495][11] = 105;
    v[496][11] = 11;
    v[497][11] = 1077;
    v[498][11] = 2857;
    v[499][11] = 337;
    v[500][11] = 3553;
    v[501][11] = 3503;
    v[502][11] = 3917;
    v[503][11] = 2665;
    v[504][11] = 3823;
    v[505][11] = 3403;
    v[506][11] = 3711;
    v[507][11] = 2085;
    v[508][11] = 1103;
    v[509][11] = 1641;
    v[510][11] = 701;
    v[511][11] = 4095;
    v[512][11] = 2883;
    v[513][11] = 1435;
    v[514][11] = 653;
    v[515][11] = 2363;
    v[516][11] = 1597;
    v[517][11] = 767;
    v[518][11] = 869;
    v[519][11] = 1825;
    v[520][11] = 1117;
    v[521][11] = 1297;
    v[522][11] = 501;
    v[523][11] = 505;
    v[524][11] = 149;
    v[525][11] = 873;
    v[526][11] = 2673;
    v[527][11] = 551;
    v[528][11] = 1499;
    v[529][11] = 2793;
    v[530][11] = 3277;
    v[531][11] = 2143;
    v[532][11] = 3663;
    v[533][11] = 533;
    v[534][11] = 3991;
    v[535][11] = 575;
    v[536][11] = 1877;
    v[537][11] = 1009;
    v[538][11] = 3929;
    v[539][11] = 473;
    v[540][11] = 3009;
    v[541][11] = 2595;
    v[542][11] = 3249;
    v[543][11] = 675;
    v[544][11] = 3593;
    v[545][11] = 2453;
    v[546][11] = 1567;
    v[547][11] = 973;
    v[548][11] = 595;
    v[549][11] = 1335;
    v[550][11] = 1715;
    v[551][11] = 589;
    v[552][11] = 85;
    v[553][11] = 2265;
    v[554][11] = 3069;
    v[555][11] = 461;
    v[556][11] = 1659;
    v[557][11] = 2627;
    v[558][11] = 1307;
    v[559][11] = 1731;
    v[560][11] = 1501;
    v[561][11] = 1699;
    v[562][11] = 3545;
    v[563][11] = 3803;
    v[564][11] = 2157;
    v[565][11] = 453;
    v[566][11] = 2813;
    v[567][11] = 2047;
    v[568][11] = 2999;
    v[569][11] = 3841;
    v[570][11] = 2361;
    v[571][11] = 1079;
    v[572][11] = 573;
    v[573][11] = 69;
    v[574][11] = 1363;
    v[575][11] = 1597;
    v[576][11] = 3427;
    v[577][11] = 2899;
    v[578][11] = 2771;
    v[579][11] = 1327;
    v[580][11] = 1117;
    v[581][11] = 1523;
    v[582][11] = 3521;
    v[583][11] = 2393;
    v[584][11] = 2537;
    v[585][11] = 1979;
    v[586][11] = 3179;
    v[587][11] = 683;
    v[588][11] = 2453;
    v[589][11] = 453;
    v[590][11] = 1227;
    v[591][11] = 779;
    v[592][11] = 671;
    v[593][11] = 3483;
    v[594][11] = 2135;
    v[595][11] = 3139;
    v[596][11] = 3381;
    v[597][11] = 3945;
    v[598][11] = 57;
    v[599][11] = 1541;
    v[600][11] = 3405;
    v[601][11] = 3381;
    v[602][11] = 2371;
    v[603][11] = 2879;
    v[604][11] = 1985;
    v[605][11] = 987;
    v[606][11] = 3017;
    v[607][11] = 3031;
    v[608][11] = 3839;
    v[609][11] = 1401;
    v[610][11] = 3749;
    v[611][11] = 2977;
    v[612][11] = 681;
    v[613][11] = 1175;
    v[614][11] = 1519;
    v[615][11] = 3355;
    v[616][11] = 907;
    v[617][11] = 117;
    v[618][11] = 771;
    v[619][11] = 3741;
    v[620][11] = 3337;
    v[621][11] = 1743;
    v[622][11] = 1227;
    v[623][11] = 3335;
    v[624][11] = 2755;
    v[625][11] = 1909;
    v[626][11] = 3603;
    v[627][11] = 2397;
    v[628][11] = 653;
    v[629][11] = 87;
    v[630][11] = 2025;
    v[631][11] = 2617;
    v[632][11] = 3257;
    v[633][11] = 287;
    v[634][11] = 3051;
    v[635][11] = 3809;
    v[636][11] = 897;
    v[637][11] = 2215;
    v[638][11] = 63;
    v[639][11] = 2043;
    v[640][11] = 1757;
    v[641][11] = 3671;
    v[642][11] = 297;
    v[643][11] = 3131;
    v[644][11] = 1305;
    v[645][11] = 293;
    v[646][11] = 3865;
    v[647][11] = 3173;
    v[648][11] = 3397;
    v[649][11] = 2269;
    v[650][11] = 3673;
    v[651][11] = 717;
    v[652][11] = 3041;
    v[653][11] = 3341;
    v[654][11] = 3595;
    v[655][11] = 3819;
    v[656][11] = 2871;
    v[657][11] = 3973;
    v[658][11] = 1129;
    v[659][11] = 513;
    v[660][11] = 871;
    v[661][11] = 1485;
    v[662][11] = 3977;
    v[663][11] = 2473;
    v[664][11] = 1171;
    v[665][11] = 1143;
    v[666][11] = 3063;
    v[667][11] = 3547;
    v[668][11] = 2183;
    v[669][11] = 3993;
    v[670][11] = 133;
    v[671][11] = 2529;
    v[672][11] = 2699;
    v[673][11] = 233;
    v[674][11] = 2355;
    v[675][11] = 231;
    v[676][11] = 3241;
    v[677][11] = 611;
    v[678][11] = 1309;
    v[679][11] = 3829;
    v[680][11] = 1839;
    v[681][11] = 1495;
    v[682][11] = 301;
    v[683][11] = 1169;
    v[684][11] = 1613;
    v[685][11] = 2673;
    v[686][11] = 243;
    v[687][11] = 3601;
    v[688][11] = 3669;
    v[689][11] = 2813;
    v[690][11] = 2671;
    v[691][11] = 2679;
    v[692][11] = 3463;
    v[693][11] = 2477;
    v[694][11] = 1795;
    v[695][11] = 617;
    v[696][11] = 2317;
    v[697][11] = 1855;
    v[698][11] = 1057;
    v[699][11] = 1703;
    v[700][11] = 1761;
    v[701][11] = 2515;
    v[702][11] = 801;
    v[703][11] = 1205;
    v[704][11] = 1311;
    v[705][11] = 473;
    v[706][11] = 3963;
    v[707][11] = 697;
    v[708][11] = 1221;
    v[709][11] = 251;
    v[710][11] = 381;
    v[711][11] = 3887;
    v[712][11] = 1761;
    v[713][11] = 3093;
    v[714][11] = 3721;
    v[715][11] = 2079;
    v[716][11] = 4085;
    v[717][11] = 379;
    v[718][11] = 3601;
    v[719][11] = 3845;
    v[720][11] = 433;
    v[721][11] = 1781;
    v[722][11] = 29;
    v[723][11] = 1897;
    v[724][11] = 1599;
    v[725][11] = 2163;
    v[726][11] = 75;
    v[727][11] = 3475;
    v[728][11] = 3957;
    v[729][11] = 1641;
    v[730][11] = 3911;
    v[731][11] = 2959;
    v[732][11] = 2833;
    v[733][11] = 1279;
    v[734][11] = 1099;
    v[735][11] = 403;
    v[736][11] = 799;
    v[737][11] = 2183;
    v[738][11] = 2699;
    v[739][11] = 1711;
    v[740][11] = 2037;
    v[741][11] = 727;
    v[742][11] = 289;
    v[743][11] = 1785;
    v[744][11] = 1575;
    v[745][11] = 3633;
    v[746][11] = 2367;
    v[747][11] = 1261;
    v[748][11] = 3953;
    v[749][11] = 1735;
    v[750][11] = 171;
    v[751][11] = 1959;
    v[752][11] = 2867;
    v[753][11] = 859;
    v[754][11] = 2951;
    v[755][11] = 3211;
    v[756][11] = 15;
    v[757][11] = 1279;
    v[758][11] = 1323;
    v[759][11] = 599;
    v[760][11] = 1651;
    v[761][11] = 3951;
    v[762][11] = 1011;
    v[763][11] = 315;
    v[764][11] = 3513;
    v[765][11] = 3351;
    v[766][11] = 1725;
    v[767][11] = 3793;
    v[768][11] = 2399;
    v[769][11] = 287;
    v[770][11] = 4017;
    v[771][11] = 3571;
    v[772][11] = 1007;
    v[773][11] = 541;
    v[774][11] = 3115;
    v[775][11] = 429;
    v[776][11] = 1585;
    v[777][11] = 1285;
    v[778][11] = 755;
    v[779][11] = 1211;
    v[780][11] = 3047;
    v[781][11] = 915;
    v[782][11] = 3611;
    v[783][11] = 2697;
    v[784][11] = 2129;
    v[785][11] = 3669;
    v[786][11] = 81;
    v[787][11] = 3939;
    v[788][11] = 2437;
    v[789][11] = 915;
    v[790][11] = 779;
    v[791][11] = 3567;
    v[792][11] = 3701;
    v[793][11] = 2479;
    v[794][11] = 3807;
    v[795][11] = 1893;
    v[796][11] = 3927;
    v[797][11] = 2619;
    v[798][11] = 2543;
    v[799][11] = 3633;
    v[800][11] = 2007;
    v[801][11] = 3857;
    v[802][11] = 3837;
    v[803][11] = 487;
    v[804][11] = 1769;
    v[805][11] = 3759;
    v[806][11] = 3105;
    v[807][11] = 2727;
    v[808][11] = 3155;
    v[809][11] = 2479;
    v[810][11] = 1341;
    v[811][11] = 1657;
    v[812][11] = 2767;
    v[813][11] = 2541;
    v[814][11] = 577;
    v[815][11] = 2105;
    v[816][11] = 799;
    v[817][11] = 17;
    v[818][11] = 2871;
    v[819][11] = 3637;
    v[820][11] = 953;
    v[821][11] = 65;
    v[822][11] = 69;
    v[823][11] = 2897;
    v[824][11] = 3841;
    v[825][11] = 3559;
    v[826][11] = 4067;
    v[827][11] = 2335;
    v[828][11] = 3409;
    v[829][11] = 1087;
    v[830][11] = 425;
    v[831][11] = 2813;
    v[832][11] = 1705;
    v[833][11] = 1701;
    v[834][11] = 1237;
    v[835][11] = 821;
    v[836][11] = 1375;
    v[837][11] = 3673;
    v[838][11] = 2693;
    v[839][11] = 3925;
    v[840][11] = 1541;
    v[841][11] = 1871;
    v[842][11] = 2285;
    v[843][11] = 847;
    v[844][11] = 4035;
    v[845][11] = 1101;
    v[846][11] = 2029;
    v[847][11] = 855;
    v[848][11] = 2733;
    v[849][11] = 2503;
    v[850][11] = 121;
    v[851][11] = 2855;
    v[852][11] = 1069;
    v[853][11] = 3463;
    v[854][11] = 3505;
    v[855][11] = 1539;
    v[856][11] = 607;
    v[857][11] = 1349;
    v[858][11] = 575;
    v[859][11] = 2301;
    v[860][11] = 2321;
    v[861][11] = 1101;
    v[862][11] = 333;
    v[863][11] = 291;
    v[864][11] = 2171;
    v[865][11] = 4085;
    v[866][11] = 2173;
    v[867][11] = 2541;
    v[868][11] = 1195;
    v[869][11] = 925;
    v[870][11] = 4039;
    v[871][11] = 1379;
    v[872][11] = 699;
    v[873][11] = 1979;
    v[874][11] = 275;
    v[875][11] = 953;
    v[876][11] = 1755;
    v[877][11] = 1643;
    v[878][11] = 325;
    v[879][11] = 101;
    v[880][11] = 2263;
    v[881][11] = 3329;
    v[882][11] = 3673;
    v[883][11] = 3413;
    v[884][11] = 1977;
    v[885][11] = 2727;
    v[886][11] = 2313;
    v[887][11] = 1419;
    v[888][11] = 887;
    v[889][11] = 609;
    v[890][11] = 2475;
    v[891][11] = 591;
    v[892][11] = 2613;
    v[893][11] = 2081;
    v[894][11] = 3805;
    v[895][11] = 3435;
    v[896][11] = 2409;
    v[897][11] = 111;
    v[898][11] = 3557;
    v[899][11] = 3607;
    v[900][11] = 903;
    v[901][11] = 231;
    v[902][11] = 3059;
    v[903][11] = 473;
    v[904][11] = 2959;
    v[905][11] = 2925;
    v[906][11] = 3861;
    v[907][11] = 2043;
    v[908][11] = 3887;
    v[909][11] = 351;
    v[910][11] = 2865;
    v[911][11] = 369;
    v[912][11] = 1377;
    v[913][11] = 2639;
    v[914][11] = 1261;
    v[915][11] = 3625;
    v[916][11] = 3279;
    v[917][11] = 2201;
    v[918][11] = 2949;
    v[919][11] = 3049;
    v[920][11] = 449;
    v[921][11] = 1297;
    v[922][11] = 897;
    v[923][11] = 1891;
    v[924][11] = 411;
    v[925][11] = 2773;
    v[926][11] = 749;
    v[927][11] = 2753;
    v[928][11] = 1825;
    v[929][11] = 853;
    v[930][11] = 2775;
    v[931][11] = 3547;
    v[932][11] = 3923;
    v[933][11] = 3923;
    v[934][11] = 987;
    v[935][11] = 3723;
    v[936][11] = 2189;
    v[937][11] = 3877;
    v[938][11] = 3577;
    v[939][11] = 297;
    v[940][11] = 2763;
    v[941][11] = 1845;
    v[942][11] = 3083;
    v[943][11] = 2951;
    v[944][11] = 483;
    v[945][11] = 2169;
    v[946][11] = 3985;
    v[947][11] = 245;
    v[948][11] = 3655;
    v[949][11] = 3441;
    v[950][11] = 1023;
    v[951][11] = 235;
    v[952][11] = 835;
    v[953][11] = 3693;
    v[954][11] = 3585;
    v[955][11] = 327;
    v[956][11] = 1003;
    v[957][11] = 543;
    v[958][11] = 3059;
    v[959][11] = 2637;
    v[960][11] = 2923;
    v[961][11] = 87;
    v[962][11] = 3617;
    v[963][11] = 1031;
    v[964][11] = 1043;
    v[965][11] = 903;
    v[966][11] = 2913;
    v[967][11] = 2177;
    v[968][11] = 2641;
    v[969][11] = 3279;
    v[970][11] = 389;
    v[971][11] = 2009;
    v[972][11] = 525;
    v[973][11] = 4085;
    v[974][11] = 3299;
    v[975][11] = 987;
    v[976][11] = 2409;
    v[977][11] = 813;
    v[978][11] = 2683;
    v[979][11] = 373;
    v[980][11] = 2695;
    v[981][11] = 3775;
    v[982][11] = 2375;
    v[983][11] = 1119;
    v[984][11] = 2791;
    v[985][11] = 223;
    v[986][11] = 325;
    v[987][11] = 587;
    v[988][11] = 1379;
    v[989][11] = 2877;
    v[990][11] = 2867;
    v[991][11] = 3793;
    v[992][11] = 655;
    v[993][11] = 831;
    v[994][11] = 3425;
    v[995][11] = 1663;
    v[996][11] = 1681;
    v[997][11] = 2657;
    v[998][11] = 1865;
    v[999][11] = 3943;
    v[1000][11] = 2977;
    v[1001][11] = 1979;
    v[1002][11] = 2271;
    v[1003][11] = 3247;
    v[1004][11] = 1267;
    v[1005][11] = 1747;
    v[1006][11] = 811;
    v[1007][11] = 159;
    v[1008][11] = 429;
    v[1009][11] = 2001;
    v[1010][11] = 1195;
    v[1011][11] = 3065;
    v[1012][11] = 553;
    v[1013][11] = 1499;
    v[1014][11] = 3529;
    v[1015][11] = 1081;
    v[1016][11] = 2877;
    v[1017][11] = 3077;
    v[1018][11] = 845;
    v[1019][11] = 1793;
    v[1020][11] = 2409;
    v[1021][11] = 3995;
    v[1022][11] = 2559;
    v[1023][11] = 4081;
    v[1024][11] = 1195;
    v[1025][11] = 2955;
    v[1026][11] = 1117;
    v[1027][11] = 1409;
    v[1028][11] = 785;
    v[1029][11] = 287;
    v[1030][11] = 1521;
    v[1031][11] = 1607;
    v[1032][11] = 85;
    v[1033][11] = 3055;
    v[1034][11] = 3123;
    v[1035][11] = 2533;
    v[1036][11] = 2329;
    v[1037][11] = 3477;
    v[1038][11] = 799;
    v[1039][11] = 3683;
    v[1040][11] = 3715;
    v[1041][11] = 337;
    v[1042][11] = 3139;
    v[1043][11] = 3311;
    v[1044][11] = 431;
    v[1045][11] = 3511;
    v[1046][11] = 2299;
    v[1047][11] = 365;
    v[1048][11] = 2941;
    v[1049][11] = 3067;
    v[1050][11] = 1331;
    v[1051][11] = 1081;
    v[1052][11] = 1097;
    v[1053][11] = 2853;
    v[1054][11] = 2299;
    v[1055][11] = 495;
    v[1056][11] = 1745;
    v[1057][11] = 749;
    v[1058][11] = 3819;
    v[1059][11] = 619;
    v[1060][11] = 1059;
    v[1061][11] = 3559;
    v[1062][11] = 183;
    v[1063][11] = 3743;
    v[1064][11] = 723;
    v[1065][11] = 949;
    v[1066][11] = 3501;
    v[1067][11] = 733;
    v[1068][11] = 2599;
    v[1069][11] = 3983;
    v[1070][11] = 3961;
    v[1071][11] = 911;
    v[1072][11] = 1899;
    v[1073][11] = 985;
    v[1074][11] = 2493;
    v[1075][11] = 1795;
    v[1076][11] = 653;
    v[1077][11] = 157;
    v[1078][11] = 433;
    v[1079][11] = 2361;
    v[1080][11] = 3093;
    v[1081][11] = 3119;
    v[1082][11] = 3679;
    v[1083][11] = 2367;
    v[1084][11] = 1701;
    v[1085][11] = 1445;
    v[1086][11] = 1321;
    v[1087][11] = 2397;
    v[1088][11] = 1241;
    v[1089][11] = 3305;
    v[1090][11] = 3985;
    v[1091][11] = 2349;
    v[1092][11] = 4067;
    v[1093][11] = 3805;
    v[1094][11] = 3073;
    v[1095][11] = 2837;
    v[1096][11] = 1567;
    v[1097][11] = 3783;
    v[1098][11] = 451;
    v[1099][11] = 2441;
    v[1100][11] = 1181;
    v[1101][11] = 487;
    v[1102][11] = 543;
    v[1103][11] = 1201;
    v[1104][11] = 3735;
    v[1105][11] = 2517;
    v[1106][11] = 733;
    v[1107][11] = 1535;
    v[1108][11] = 2175;
    v[1109][11] = 3613;
    v[1110][11] = 3019;

    v[481][12] = 2319;
    v[482][12] = 653;
    v[483][12] = 1379;
    v[484][12] = 1675;
    v[485][12] = 1951;
    v[486][12] = 7075;
    v[487][12] = 2087;
    v[488][12] = 7147;
    v[489][12] = 1427;
    v[490][12] = 893;
    v[491][12] = 171;
    v[492][12] = 2019;
    v[493][12] = 7235;
    v[494][12] = 5697;
    v[495][12] = 3615;
    v[496][12] = 1961;
    v[497][12] = 7517;
    v[498][12] = 6849;
    v[499][12] = 2893;
    v[500][12] = 1883;
    v[501][12] = 2863;
    v[502][12] = 2173;
    v[503][12] = 4543;
    v[504][12] = 73;
    v[505][12] = 381;
    v[506][12] = 3893;
    v[507][12] = 6045;
    v[508][12] = 1643;
    v[509][12] = 7669;
    v[510][12] = 1027;
    v[511][12] = 1549;
    v[512][12] = 3983;
    v[513][12] = 1985;
    v[514][12] = 6589;
    v[515][12] = 7497;
    v[516][12] = 2745;
    v[517][12] = 2375;
    v[518][12] = 7047;
    v[519][12] = 1117;
    v[520][12] = 1171;
    v[521][12] = 1975;
    v[522][12] = 5199;
    v[523][12] = 3915;
    v[524][12] = 3695;
    v[525][12] = 8113;
    v[526][12] = 4303;
    v[527][12] = 3773;
    v[528][12] = 7705;
    v[529][12] = 6855;
    v[530][12] = 1675;
    v[531][12] = 2245;
    v[532][12] = 2817;
    v[533][12] = 1719;
    v[534][12] = 569;
    v[535][12] = 1021;
    v[536][12] = 2077;
    v[537][12] = 5945;
    v[538][12] = 1833;
    v[539][12] = 2631;
    v[540][12] = 4851;
    v[541][12] = 6371;
    v[542][12] = 833;
    v[543][12] = 7987;
    v[544][12] = 331;
    v[545][12] = 1899;
    v[546][12] = 8093;
    v[547][12] = 6719;
    v[548][12] = 6903;
    v[549][12] = 5903;
    v[550][12] = 5657;
    v[551][12] = 5007;
    v[552][12] = 2689;
    v[553][12] = 6637;
    v[554][12] = 2675;
    v[555][12] = 1645;
    v[556][12] = 1819;
    v[557][12] = 689;
    v[558][12] = 6709;
    v[559][12] = 7717;
    v[560][12] = 6295;
    v[561][12] = 7013;
    v[562][12] = 7695;
    v[563][12] = 3705;
    v[564][12] = 7069;
    v[565][12] = 2621;
    v[566][12] = 3631;
    v[567][12] = 6571;
    v[568][12] = 6259;
    v[569][12] = 7261;
    v[570][12] = 3397;
    v[571][12] = 7645;
    v[572][12] = 1115;
    v[573][12] = 4753;
    v[574][12] = 2047;
    v[575][12] = 7579;
    v[576][12] = 2271;
    v[577][12] = 5403;
    v[578][12] = 4911;
    v[579][12] = 7629;
    v[580][12] = 4225;
    v[581][12] = 1209;
    v[582][12] = 6955;
    v[583][12] = 6951;
    v[584][12] = 1829;
    v[585][12] = 5579;
    v[586][12] = 5231;
    v[587][12] = 1783;
    v[588][12] = 4285;
    v[589][12] = 7425;
    v[590][12] = 599;
    v[591][12] = 5785;
    v[592][12] = 3275;
    v[593][12] = 5643;
    v[594][12] = 2263;
    v[595][12] = 657;
    v[596][12] = 6769;
    v[597][12] = 6261;
    v[598][12] = 1251;
    v[599][12] = 3249;
    v[600][12] = 4447;
    v[601][12] = 4111;
    v[602][12] = 3991;
    v[603][12] = 1215;
    v[604][12] = 131;
    v[605][12] = 4397;
    v[606][12] = 3487;
    v[607][12] = 7585;
    v[608][12] = 5565;
    v[609][12] = 7199;
    v[610][12] = 3573;
    v[611][12] = 7105;
    v[612][12] = 7409;
    v[613][12] = 1671;
    v[614][12] = 949;
    v[615][12] = 3889;
    v[616][12] = 5971;
    v[617][12] = 3333;
    v[618][12] = 225;
    v[619][12] = 3647;
    v[620][12] = 5403;
    v[621][12] = 3409;
    v[622][12] = 7459;
    v[623][12] = 6879;
    v[624][12] = 5789;
    v[625][12] = 6567;
    v[626][12] = 5581;
    v[627][12] = 4919;
    v[628][12] = 1927;
    v[629][12] = 4407;
    v[630][12] = 8085;
    v[631][12] = 4691;
    v[632][12] = 611;
    v[633][12] = 3005;
    v[634][12] = 591;
    v[635][12] = 753;
    v[636][12] = 589;
    v[637][12] = 171;
    v[638][12] = 5729;
    v[639][12] = 5891;
    v[640][12] = 1033;
    v[641][12] = 3049;
    v[642][12] = 6567;
    v[643][12] = 5257;
    v[644][12] = 8003;
    v[645][12] = 1757;
    v[646][12] = 4489;
    v[647][12] = 4923;
    v[648][12] = 6379;
    v[649][12] = 5171;
    v[650][12] = 1757;
    v[651][12] = 689;
    v[652][12] = 3081;
    v[653][12] = 1389;
    v[654][12] = 4113;
    v[655][12] = 455;
    v[656][12] = 2761;
    v[657][12] = 847;
    v[658][12] = 7575;
    v[659][12] = 5829;
    v[660][12] = 633;
    v[661][12] = 6629;
    v[662][12] = 1103;
    v[663][12] = 7635;
    v[664][12] = 803;
    v[665][12] = 6175;
    v[666][12] = 6587;
    v[667][12] = 2711;
    v[668][12] = 3879;
    v[669][12] = 67;
    v[670][12] = 1179;
    v[671][12] = 4761;
    v[672][12] = 7281;
    v[673][12] = 1557;
    v[674][12] = 3379;
    v[675][12] = 2459;
    v[676][12] = 4273;
    v[677][12] = 4127;
    v[678][12] = 7147;
    v[679][12] = 35;
    v[680][12] = 3549;
    v[681][12] = 395;
    v[682][12] = 3735;
    v[683][12] = 5787;
    v[684][12] = 4179;
    v[685][12] = 5889;
    v[686][12] = 5057;
    v[687][12] = 7473;
    v[688][12] = 4713;
    v[689][12] = 2133;
    v[690][12] = 2897;
    v[691][12] = 1841;
    v[692][12] = 2125;
    v[693][12] = 1029;
    v[694][12] = 1695;
    v[695][12] = 6523;
    v[696][12] = 1143;
    v[697][12] = 5105;
    v[698][12] = 7133;
    v[699][12] = 3351;
    v[700][12] = 2775;
    v[701][12] = 3971;
    v[702][12] = 4503;
    v[703][12] = 7589;
    v[704][12] = 5155;
    v[705][12] = 4305;
    v[706][12] = 1641;
    v[707][12] = 4717;
    v[708][12] = 2427;
    v[709][12] = 5617;
    v[710][12] = 1267;
    v[711][12] = 399;
    v[712][12] = 5831;
    v[713][12] = 4305;
    v[714][12] = 4241;
    v[715][12] = 3395;
    v[716][12] = 3045;
    v[717][12] = 4899;
    v[718][12] = 1713;
    v[719][12] = 171;
    v[720][12] = 411;
    v[721][12] = 7099;
    v[722][12] = 5473;
    v[723][12] = 5209;
    v[724][12] = 1195;
    v[725][12] = 1077;
    v[726][12] = 1309;
    v[727][12] = 2953;
    v[728][12] = 7343;
    v[729][12] = 4887;
    v[730][12] = 3229;
    v[731][12] = 6759;
    v[732][12] = 6721;
    v[733][12] = 6775;
    v[734][12] = 675;
    v[735][12] = 4039;
    v[736][12] = 2493;
    v[737][12] = 7511;
    v[738][12] = 3269;
    v[739][12] = 4199;
    v[740][12] = 6625;
    v[741][12] = 7943;
    v[742][12] = 2013;
    v[743][12] = 4145;
    v[744][12] = 667;
    v[745][12] = 513;
    v[746][12] = 2303;
    v[747][12] = 4591;
    v[748][12] = 7941;
    v[749][12] = 2741;
    v[750][12] = 987;
    v[751][12] = 8061;
    v[752][12] = 3161;
    v[753][12] = 5951;
    v[754][12] = 1431;
    v[755][12] = 831;
    v[756][12] = 5559;
    v[757][12] = 7405;
    v[758][12] = 1357;
    v[759][12] = 4319;
    v[760][12] = 4235;
    v[761][12] = 5421;
    v[762][12] = 2559;
    v[763][12] = 4415;
    v[764][12] = 2439;
    v[765][12] = 823;
    v[766][12] = 1725;
    v[767][12] = 6219;
    v[768][12] = 4903;
    v[769][12] = 6699;
    v[770][12] = 5451;
    v[771][12] = 349;
    v[772][12] = 7703;
    v[773][12] = 2927;
    v[774][12] = 7809;
    v[775][12] = 6179;
    v[776][12] = 1417;
    v[777][12] = 5987;
    v[778][12] = 3017;
    v[779][12] = 4983;
    v[780][12] = 3479;
    v[781][12] = 4525;
    v[782][12] = 4643;
    v[783][12] = 4911;
    v[784][12] = 227;
    v[785][12] = 5475;
    v[786][12] = 2287;
    v[787][12] = 5581;
    v[788][12] = 6817;
    v[789][12] = 1937;
    v[790][12] = 1421;
    v[791][12] = 4415;
    v[792][12] = 7977;
    v[793][12] = 1789;
    v[794][12] = 3907;
    v[795][12] = 6815;
    v[796][12] = 6789;
    v[797][12] = 6003;
    v[798][12] = 5609;
    v[799][12] = 4507;
    v[800][12] = 337;
    v[801][12] = 7427;
    v[802][12] = 7943;
    v[803][12] = 3075;
    v[804][12] = 6427;
    v[805][12] = 1019;
    v[806][12] = 7121;
    v[807][12] = 4763;
    v[808][12] = 81;
    v[809][12] = 3587;
    v[810][12] = 2929;
    v[811][12] = 1795;
    v[812][12] = 8067;
    v[813][12] = 2415;
    v[814][12] = 1265;
    v[815][12] = 4025;
    v[816][12] = 5599;
    v[817][12] = 4771;
    v[818][12] = 3025;
    v[819][12] = 2313;
    v[820][12] = 6129;
    v[821][12] = 7611;
    v[822][12] = 6881;
    v[823][12] = 5253;
    v[824][12] = 4413;
    v[825][12] = 7869;
    v[826][12] = 105;
    v[827][12] = 3173;
    v[828][12] = 1629;
    v[829][12] = 2537;
    v[830][12] = 1023;
    v[831][12] = 4409;
    v[832][12] = 7209;
    v[833][12] = 4413;
    v[834][12] = 7107;
    v[835][12] = 7469;
    v[836][12] = 33;
    v[837][12] = 1955;
    v[838][12] = 2881;
    v[839][12] = 5167;
    v[840][12] = 6451;
    v[841][12] = 4211;
    v[842][12] = 179;
    v[843][12] = 5573;
    v[844][12] = 7879;
    v[845][12] = 3387;
    v[846][12] = 7759;
    v[847][12] = 5455;
    v[848][12] = 7157;
    v[849][12] = 1891;
    v[850][12] = 5683;
    v[851][12] = 5689;
    v[852][12] = 6535;
    v[853][12] = 3109;
    v[854][12] = 6555;
    v[855][12] = 6873;
    v[856][12] = 1249;
    v[857][12] = 4251;
    v[858][12] = 6437;
    v[859][12] = 49;
    v[860][12] = 2745;
    v[861][12] = 1201;
    v[862][12] = 7327;
    v[863][12] = 4179;
    v[864][12] = 6783;
    v[865][12] = 623;
    v[866][12] = 2779;
    v[867][12] = 5963;
    v[868][12] = 2585;
    v[869][12] = 6927;
    v[870][12] = 5333;
    v[871][12] = 4033;
    v[872][12] = 285;
    v[873][12] = 7467;
    v[874][12] = 4443;
    v[875][12] = 4917;
    v[876][12] = 3;
    v[877][12] = 4319;
    v[878][12] = 5517;
    v[879][12] = 3449;
    v[880][12] = 813;
    v[881][12] = 5499;
    v[882][12] = 2515;
    v[883][12] = 5771;
    v[884][12] = 3357;
    v[885][12] = 2073;
    v[886][12] = 4395;
    v[887][12] = 4925;
    v[888][12] = 2643;
    v[889][12] = 7215;
    v[890][12] = 5817;
    v[891][12] = 1199;
    v[892][12] = 1597;
    v[893][12] = 1619;
    v[894][12] = 7535;
    v[895][12] = 4833;
    v[896][12] = 609;
    v[897][12] = 4797;
    v[898][12] = 8171;
    v[899][12] = 6847;
    v[900][12] = 793;
    v[901][12] = 6757;
    v[902][12] = 8165;
    v[903][12] = 3371;
    v[904][12] = 2431;
    v[905][12] = 5235;
    v[906][12] = 4739;
    v[907][12] = 7703;
    v[908][12] = 7223;
    v[909][12] = 6525;
    v[910][12] = 5891;
    v[911][12] = 5605;
    v[912][12] = 4433;
    v[913][12] = 3533;
    v[914][12] = 5267;
    v[915][12] = 5125;
    v[916][12] = 5037;
    v[917][12] = 225;
    v[918][12] = 6717;
    v[919][12] = 1121;
    v[920][12] = 5741;
    v[921][12] = 2013;
    v[922][12] = 4327;
    v[923][12] = 4839;
    v[924][12] = 569;
    v[925][12] = 5227;
    v[926][12] = 7677;
    v[927][12] = 4315;
    v[928][12] = 2391;
    v[929][12] = 5551;
    v[930][12] = 859;
    v[931][12] = 3627;
    v[932][12] = 6377;
    v[933][12] = 3903;
    v[934][12] = 4311;
    v[935][12] = 6527;
    v[936][12] = 7573;
    v[937][12] = 4905;
    v[938][12] = 7731;
    v[939][12] = 1909;
    v[940][12] = 1555;
    v[941][12] = 3279;
    v[942][12] = 1949;
    v[943][12] = 1887;
    v[944][12] = 6675;
    v[945][12] = 5509;
    v[946][12] = 2033;
    v[947][12] = 5473;
    v[948][12] = 3539;
    v[949][12] = 5033;
    v[950][12] = 5935;
    v[951][12] = 6095;
    v[952][12] = 4761;
    v[953][12] = 1771;
    v[954][12] = 1271;
    v[955][12] = 1717;
    v[956][12] = 4415;
    v[957][12] = 5083;
    v[958][12] = 6277;
    v[959][12] = 3147;
    v[960][12] = 7695;
    v[961][12] = 2461;
    v[962][12] = 4783;
    v[963][12] = 4539;
    v[964][12] = 5833;
    v[965][12] = 5583;
    v[966][12] = 651;
    v[967][12] = 1419;
    v[968][12] = 2605;
    v[969][12] = 5511;
    v[970][12] = 3913;
    v[971][12] = 5795;
    v[972][12] = 2333;
    v[973][12] = 2329;
    v[974][12] = 4431;
    v[975][12] = 3725;
    v[976][12] = 6069;
    v[977][12] = 2699;
    v[978][12] = 7055;
    v[979][12] = 6879;
    v[980][12] = 1017;
    v[981][12] = 3121;
    v[982][12] = 2547;
    v[983][12] = 4603;
    v[984][12] = 2385;
    v[985][12] = 6915;
    v[986][12] = 6103;
    v[987][12] = 5669;
    v[988][12] = 7833;
    v[989][12] = 2001;
    v[990][12] = 4287;
    v[991][12] = 6619;
    v[992][12] = 955;
    v[993][12] = 2761;
    v[994][12] = 5711;
    v[995][12] = 6291;
    v[996][12] = 3415;
    v[997][12] = 3909;
    v[998][12] = 2841;
    v[999][12] = 5627;
    v[1000][12] = 4939;
    v[1001][12] = 7671;
    v[1002][12] = 6059;
    v[1003][12] = 6275;
    v[1004][12] = 6517;
    v[1005][12] = 1931;
    v[1006][12] = 4583;
    v[1007][12] = 7301;
    v[1008][12] = 1267;
    v[1009][12] = 7509;
    v[1010][12] = 1435;
    v[1011][12] = 2169;
    v[1012][12] = 6939;
    v[1013][12] = 3515;
    v[1014][12] = 2985;
    v[1015][12] = 2787;
    v[1016][12] = 2123;
    v[1017][12] = 1969;
    v[1018][12] = 3307;
    v[1019][12] = 353;
    v[1020][12] = 4359;
    v[1021][12] = 7059;
    v[1022][12] = 5273;
    v[1023][12] = 5873;
    v[1024][12] = 6657;
    v[1025][12] = 6765;
    v[1026][12] = 6229;
    v[1027][12] = 3179;
    v[1028][12] = 1583;
    v[1029][12] = 6237;
    v[1030][12] = 2155;
    v[1031][12] = 371;
    v[1032][12] = 273;
    v[1033][12] = 7491;
    v[1034][12] = 3309;
    v[1035][12] = 6805;
    v[1036][12] = 3015;
    v[1037][12] = 6831;
    v[1038][12] = 7819;
    v[1039][12] = 713;
    v[1040][12] = 4747;
    v[1041][12] = 3935;
    v[1042][12] = 4109;
    v[1043][12] = 1311;
    v[1044][12] = 709;
    v[1045][12] = 3089;
    v[1046][12] = 7059;
    v[1047][12] = 4247;
    v[1048][12] = 2989;
    v[1049][12] = 1509;
    v[1050][12] = 4919;
    v[1051][12] = 1841;
    v[1052][12] = 3045;
    v[1053][12] = 3821;
    v[1054][12] = 6929;
    v[1055][12] = 4655;
    v[1056][12] = 1333;
    v[1057][12] = 6429;
    v[1058][12] = 6649;
    v[1059][12] = 2131;
    v[1060][12] = 5265;
    v[1061][12] = 1051;
    v[1062][12] = 261;
    v[1063][12] = 8057;
    v[1064][12] = 3379;
    v[1065][12] = 2179;
    v[1066][12] = 1993;
    v[1067][12] = 5655;
    v[1068][12] = 3063;
    v[1069][12] = 6381;
    v[1070][12] = 3587;
    v[1071][12] = 7417;
    v[1072][12] = 1579;
    v[1073][12] = 1541;
    v[1074][12] = 2107;
    v[1075][12] = 5085;
    v[1076][12] = 2873;
    v[1077][12] = 6141;
    v[1078][12] = 955;
    v[1079][12] = 3537;
    v[1080][12] = 2157;
    v[1081][12] = 841;
    v[1082][12] = 1999;
    v[1083][12] = 1465;
    v[1084][12] = 5171;
    v[1085][12] = 5651;
    v[1086][12] = 1535;
    v[1087][12] = 7235;
    v[1088][12] = 4349;
    v[1089][12] = 1263;
    v[1090][12] = 1453;
    v[1091][12] = 1005;
    v[1092][12] = 6893;
    v[1093][12] = 2919;
    v[1094][12] = 1947;
    v[1095][12] = 1635;
    v[1096][12] = 3963;
    v[1097][12] = 397;
    v[1098][12] = 969;
    v[1099][12] = 4569;
    v[1100][12] = 655;
    v[1101][12] = 6737;
    v[1102][12] = 2995;
    v[1103][12] = 7235;
    v[1104][12] = 7713;
    v[1105][12] = 973;
    v[1106][12] = 4821;
    v[1107][12] = 2377;
    v[1108][12] = 1673;
    v[1109][12] = 1;
    v[1110][12] = 6541;
//
//  Check parameters.
//
    if ( dim_num < 1 || DIM_MAX2 < dim_num )
    {
      cout << "\n";
      cout << "I8_SOBOL - Fatal error!\n";
      cout << "  The spatial dimension DIM_NUM should satisfy:\n";
      cout << "    1 <= DIM_NUM <= " << DIM_MAX2 << "\n";
      cout << "  But this input value is DIM_NUM = " << dim_num << "\n";
      exit ( 1 );
    }

    dim_num_save = dim_num;
//
//  Find the number of bits in ATMOST.
//
//  Here, we have short-circuited the computation of MAXCOL from ATMOST, because
//  in some cases, a compiler was complaining that the value of ATMOST could not
//  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
//  so if we know what the answer should be we can try to get it this way!
//  JVB, 24 January 2006.
//
//  maxcol = i8_bit_hi1 ( atmost );
//
    maxcol = 62;
//
//  Initialize row 1 of V.
//
    for ( j = 0; j < maxcol; j++ )
    {
      v[0][j] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; 0 <= k; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
      for ( j = m; j < maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }
        }
        v[i][j] = newv;
      }
    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = maxcol - 2; 0 <= j; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( double ) ( 2 * l ) );
  }

  if ( *seed < 0 )
  {
    *seed = 0;
  }

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == seed_save + 1 )
  {
    l = i8_bit_lo0 ( *seed );
  }
  else if ( *seed <= seed_save )
  {
    seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i8_bit_lo0 ( *seed );
  }
  else if ( seed_save+1 < *seed )
  {
    for ( seed_temp = seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i8_bit_lo0 ( *seed );
  }
//
//  Check that the user is not calling too many times!
//
  if ( maxcol < l )
  {
    cout << "\n";
    cout << "I8_SOBOL - Fatal error!\n";
    cout << "  The value of SEED seems to be too large!\n";
    cout << "  SEED =   " << *seed  << "\n";
    cout << "  MAXCOL = " << maxcol << "\n";
    cout << "  L =      " << l << "\n";
    exit ( 2 );
  }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( double ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef DIM_MAX
# undef DIM_MAX2
# undef LOG_MAX
}

// ================================= i0_bit_lo0 ==================================

int cSamp :: i8_bit_lo0 ( long long int n )

//****************************************************************************80
//
//  Purpose:
//
//    I8_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Lo 0
//    ----    --------  ----
//       0           0     1
//       1           1     2
//       2          10     1
//       3          11     3
//       4         100     1
//       5         101     2
//       6         110     1
//       7         111     4
//       8        1000     1
//       9        1001     2
//      10        1010     1
//      11        1011     3
//      12        1100     1
//      13        1101     2
//      14        1110     1
//      15        1111     5
//      16       10000     1
//      17       10001     2
//    1023  1111111111     1
//    1024 10000000000     1
//    1025 10000000001     1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long long int N, the integer to be measured.
//    N should be nonnegative.
//
//    Output, int I8_BIT_LO0, the position of the low 1 bit.
//
{
  int bit;
  long long int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}

// ================================= CalcSampLHS ==================================

void cSamp :: CalcSampNLHS(int dim_num, int point_num, vector<cVector>& sx)

{
  vector<cVector> sxi;
  int n = 20;
  double s2;
  double s2min = 10e30;
  double maxmin = 0;

  for (int ii = 0; ii < n; ii++)
  {
      int i;
      int j;
      int *perm;
      double *x;
      cVector xs(NumVar);
      int seed;
      for (int jj = 0; jj < 10; jj++)
      {
      seed = get_seed();
      }

      x = r8mat_uniform_01_new ( dim_num, point_num, seed );
//
//  For spatial dimension I,
//    pick a random permutation of 1 to POINT_NUM,
//    force the corresponding I-th components of X to lie in the
//    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
//
      for ( i = 0; i < dim_num; i++ )
      {
        perm = perm_uniform_new ( point_num, seed );

        for ( j = 0; j < point_num; j++ )
        {
          x[i+j*dim_num] = ( ( ( double ) perm[j] ) + (double) x[i+j*dim_num] ) / ( ( double ) point_num );
        }
        delete [] perm;
      }

      for ( int j = 0; j < point_num; j++)
      {
      for ( i = 0; i < dim_num; i++ )
      {
        xs[i] = x[i+j*dim_num];
      }
      sxi.push_back(xs);
      }

      for (int i = 0; i < point_num; i++)
      {
          for (int j = 0; j < i; j++)
          {
              s2 = 0;
              for (int k = 0; k < dim_num; k++)
              {
                 s2 += pow((sxi[i][k] - sxi[j][k]) ,2);
              }
              if (s2 < s2min)
              {
                  s2min = s2;
              }
          }
      }

      if (s2min > maxmin)
      {
          maxmin = s2min;
          sx = sxi;
      }

   }
}

// ================================= CalcSampCVT ==================================

void cSamp :: CalcSampCVT( int dim_num, int n, int batch, int init, int sample, int sample_num, int it_max, int it_fixed, int *seed, double r[],int *it_num, double *it_diff, double *energy, vector<cVector> &sx  )

//****************************************************************************80
//
//  Purpose:
//
//    CVT computes a Centroidal Voronoi Tessellation.
//
//  Discussion:
//
//    This routine initializes the data, and carries out the
//    CVT iteration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int INIT, specifies how the points are to be initialized.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine;
//     4, points are already initialized on input.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Input, int IT_FIXED, the maximum number of iterations to take
//    with a fixed set of sample points.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the approximate CVT points.
//    If INIT = 4 on input, then it is assumed that these values have been
//    initialized.  On output, the CVT iteration has been applied to improve
//    the value of the points.
//
//    Output, int *IT_NUM, the number of iterations taken.  Generally,
//    this will be equal to IT_MAX, unless the iteration tolerance was
//    satisfied early.
//
//    Output, double *IT_DIFF, the L2 norm of the difference
//    between the iterates.
//
//    Output, double *ENERGY,  the discrete "energy", divided
//    by the number of sample points.
//
{
  bool DEBUG = false;
  int i;
  bool initialize;
  int seed_base;
  int seed_init;
  cVector xs(dim_num);

  if ( batch < 1 )
  {
    cout << "\n";
    cout << "CVT - Fatal error!\n";
    cout << "  The input value BATCH < 1.\n";
    exit ( 1 );
  }

  if ( *seed <= 0 )
  {
    cout << "\n";
    cout << "CVT - Fatal error!\n";
    cout << "  The input value SEED <= 0.\n";
    exit ( 1 );
  }

  if ( DEBUG )
  {
    cout << "\n";
    cout << "  Step       SEED          L2-Change        Energy\n";
    cout << "\n";
  }

  *it_num = 0;
  *it_diff = 0.0;
  *energy = 0.0;
  seed_init = *seed;
//
//  Initialize the data, unless the user has already done that.
//
  if ( init != 4 )
  {
    initialize = true;
    cvt_sample ( dim_num, n, n, init, initialize, seed, r );
  }
  if ( DEBUG )
  {
    cout                          << "  "
         << setw(4)  << *it_num   << "  "
         << setw(12) << seed_init << "\n";
  }
//
//  If the initialization and sampling steps use the same random number
//  scheme, then the sampling scheme does not have to be initialized.
//
  if ( init == sample )
  {
    initialize = false;
  }
  else
  {
    initialize = true;
  }
//
//  Carry out the iteration.
//
  while ( *it_num < it_max )
  {
//
//  If it's time to update the seed, save its current value
//  as the starting value for all iterations in this cycle.
//  If it's not time to update the seed, restore it to its initial
//  value for this cycle.
//
    if ( ( (*it_num) % it_fixed ) == 0 )
    {
      seed_base = *seed;
    }
    else
    {
      *seed = seed_base;
    }

    *it_num = *it_num + 1;
    seed_init = *seed;

    cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, seed,
      r, it_diff, energy );

    initialize = false;

    if ( DEBUG )
    {
      cout                          << "  "
           << setw(4)  << *it_num   << "  "
           << setw(12) << seed_init << "  "
           << setw(14) << *it_diff  << "  "
           << setw(14) << *energy    << "\n";
    }
  }

  for ( int j = 0; j < n; j++)
  {
  for ( int i = 0; i < dim_num; i++ )
  {
    xs[i] = r[i+j*dim_num];
  }
  sx.push_back(xs);
  }

  return;
}

// ================================= cvt_iterate ==================================

void cSamp :: cvt_iterate ( int dim_num, int n, int batch, int sample, bool initialize,
  int sample_num, int *seed, double r[], double *it_diff, double *energy )
//****************************************************************************80
//
//  Purpose:
//
//    CVT_ITERATE takes one step of the CVT iteration.
//
//  Discussion:
//
//    The routine is given a set of points, called "generators", which
//    define a tessellation of the region into Voronoi cells.  Each point
//    defines a cell.  Each cell, in turn, has a centroid, but it is
//    unlikely that the centroid and the generator coincide.
//
//    Each time this CVT iteration is carried out, an attempt is made
//    to modify the generators in such a way that they are closer and
//    closer to being the centroids of the Voronoi cells they generate.
//
//    A large number of sample points are generated, and the nearest generator
//    is determined.  A count is kept of how many points were nearest to each
//    generator.  Once the sampling is completed, the location of all the
//    generators is adjusted.  This step should decrease the discrepancy
//    between the generators and the centroids.
//
//    The centroidal Voronoi tessellation minimizes the "energy",
//    defined to be the integral, over the region, of the square of
//    the distance between each point in the region and its nearest generator.
//    The sampling technique supplies a discrete estimate of this
//    energy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, bool INITIALIZE, is TRUE if the SEED must be reset to SEED_INIT
//    before computation.  Also, the pseudorandom process may need to be
//    reinitialized.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the Voronoi
//    cell generators.  On output, these have been modified
//
//    Output, double *IT_DIFF, the L2 norm of the difference
//    between the iterates.
//
//    Output, double *ENERGY,  the discrete "energy", divided
//    by the number of sample points.
//
{
  int *count;
  int get;
  int have;
  int i;
  int j;
  int j2;
  int *nearest;
  double *r2;
  double *s;
  bool success;
  double term;
//
//  Take each generator as the first sample point for its region.
//  This can slightly slow the convergence, but it simplifies the
//  algorithm by guaranteeing that no region is completely missed
//  by the sampling.
//
  *energy = 0.0;
  r2 = new double[dim_num*n];
  count = new int[n];
  nearest = new int[sample_num];
  s = new double[dim_num*sample_num];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r2[i+j*dim_num] = r[i+j*dim_num];
    }
  }
  for ( j = 0; j < n; j++ )
  {
    count[j] = 1;
  }
//
//  Generate the sampling points S.
//
  have = 0;

  while ( have < sample_num )
  {
    get = i4_min ( sample_num - have, batch );

    cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s );

    initialize = false;
    have = have + get;
//
//  Find the index N of the nearest cell generator to each sample point S.
//
    find_closest ( dim_num, n, get, s, r, nearest );
//
//  Add S to the centroid associated with generator N.
//
    for ( j = 0; j < get; j++ )
    {
      j2 = nearest[j];
      for ( i = 0; i < dim_num; i++ )
      {
        r2[i+j2*dim_num] = r2[i+j2*dim_num] + s[i+j*dim_num];
      }
      for ( i = 0; i < dim_num; i++ )
      {
        *energy = *energy + pow ( r[i+j2*dim_num] - s[i+j*dim_num], 2 );
      }
      count[j2] = count[j2] + 1;
    }
  }
//
//  Estimate the centroids.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r2[i+j*dim_num] = r2[i+j*dim_num] / ( double ) ( count[j] );
    }
  }
//
//  Determine the sum of the distances between generators and centroids.
//
  *it_diff = 0.0;

  for ( j = 0; j < n; j++ )
  {
    term = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      term = term + ( r2[i+j*dim_num] - r[i+j*dim_num] )
                  * ( r2[i+j*dim_num] - r[i+j*dim_num] );
    }
    *it_diff = *it_diff + sqrt ( term );
  }
//
//  Replace the generators by the centroids.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r[i+j*dim_num] = r2[i+j*dim_num];
    }
  }
//
//  Normalize the discrete energy estimate.
//
  *energy = *energy / sample_num;

  delete [] count;
  delete [] nearest;
  delete [] r2;
  delete [] s;

  return;
}
// ================================= cvt_sample ==================================

void cSamp :: cvt_sample ( int dim_num, int n, int n_now, int sample, bool initialize, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_SAMPLE returns sample points.
//
//  Discussion:
//
//    N sample points are to be taken from the unit box of dimension DIM_NUM.
//
//    These sample points are usually created by a pseudorandom process
//    for which the points are essentially indexed by a quantity called
//    SEED.  To get N sample points, we generate values with indices
//    SEED through SEED+N-1.
//
//    It may not be practical to generate all the sample points in a
//    single call.  For that reason, the routine allows the user to
//    request a total of N points, but to require that only N_NOW be
//    generated now (on this call).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int N_NOW, the number of sample points to be generated
//    on this call.  N_NOW must be at least 1.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, bool INITIALIZE, is TRUE if the pseudorandom process should be
//    reinitialized.
//
//    Input/output, int *SEED, the random number seed.
//
//    Output, double R[DIM_NUM*N_NOW], the sample points.
//
{
  double exponent;
  static int *halton_base = NULL;
  static int *halton_leap = NULL;
  static int *halton_seed = NULL;
  int halton_step;
  int i;
  int j;
  int k;
  static int ngrid;
  static int rank;
  int rank_max;
  static int *tuple = NULL;

  if ( n_now < 1 )
  {
    cout << "\n";
    cout << "CVT_SAMPLE - Fatal error!\n";
    cout << "  N_NOW < 1.\n";
    exit ( 1 );
  }

  if ( sample == -1 )
  {
    if ( initialize )
    {
      random_initialize ( *seed );
    }

    for ( j = 0; j < n_now; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = (double) rand( ) / ( double ) RAND_MAX;
      }
    }
    *seed = ( *seed ) + n_now * dim_num;
  }
  else if ( sample == 0 )
  {
    r8mat_uniform_01 ( dim_num, n_now, seed, r );
  }
  else if ( sample == 1 )
  {
    halton_seed = new int[dim_num];
    halton_leap = new int[dim_num];
    halton_base = new int[dim_num];

    halton_step = *seed;

    for ( i = 0; i < dim_num; i++ )
    {
      halton_seed[i] = 0;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      halton_leap[i] = 1;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      halton_base[i] = prime ( i + 1 );
    }

    i4_to_halton_sequence ( dim_num, n_now, halton_step, halton_seed,
      halton_leap, halton_base, r );

    delete [] halton_seed;
    delete [] halton_leap;
    delete [] halton_base;

    *seed = *seed + n_now;
  }
  else if ( sample == 2 )
  {
    exponent = 1.0 / ( double ) ( dim_num );
    ngrid = ( int ) pow ( ( double ) n, exponent );
    rank_max = ( int ) pow ( ( double ) ngrid, ( double ) dim_num );
    tuple = new int[dim_num];

    if ( rank_max < n )
    {
      ngrid = ngrid + 1;
      rank_max = ( int ) pow ( ( double ) ngrid, ( double ) dim_num );
    }

    if ( initialize )
    {
      rank = -1;
      tuple_next_fast ( ngrid, dim_num, rank, tuple );
    }

    rank = ( *seed ) % rank_max;

    for ( j = 0; j < n_now; j++ )
    {
      tuple_next_fast ( ngrid, dim_num, rank, tuple );
      rank = rank + 1;
      rank = rank % rank_max;
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = double ( 2 * tuple[i] - 1 ) / double ( 2 * ngrid );
      }
    }
    delete [] tuple;
    *seed = *seed + n_now;
  }
  else if ( sample == 3 )
  {
    user ( dim_num, n_now, seed, r );
  }
  else
  {
    cout << "\n";
    cout << "CVT_SAMPLE - Fatal error!\n";
    cout << "  The value of SAMPLE = " << sample << " is illegal.\n";
    exit ( 1 );
  }

  return;
}

// ================================= i4_min ==================================

int cSamp :: i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1 and I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}

// ================================= find_closest ==================================

void cSamp :: find_closest ( int dim_num, int n, int sample_num, double s[], double r[], int nearest[] )

//****************************************************************************80
//
//  Purpose:
//
//    FIND_CLOSEST finds the nearest R point to each S point.
//
//  Discussion:
//
//    This routine finds the closest Voronoi cell generator by checking every
//    one.  For problems with many cells, this process can take the bulk
//    of the CPU time.  Other approaches, which group the cell generators into
//    bins, can run faster by a large factor.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of cell generators.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, double S[DIM_NUM*SAMPLE_NUM], the points to be checked.
//
//    Input, double R[DIM_NUM*N], the cell generators.
//
//    Output, int NEAREST[SAMPLE_NUM], the (0-based) index of the nearest
//    cell generator.
//
{
  double dist_sq_min;
  double dist_sq;
  int i;
  int jr;
  int js;

  for ( js = 0; js < sample_num; js++ )
  {
    dist_sq_min = r8_huge ( );
    nearest[js] = -1;

    for ( jr = 0; jr < n; jr++ )
    {
      dist_sq = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        dist_sq = dist_sq + ( s[i+js*dim_num] - r[i+jr*dim_num] )
                          * ( s[i+js*dim_num] - r[i+jr*dim_num] );
      }

      if ( jr == 0 || dist_sq < dist_sq_min )
      {
        dist_sq_min = dist_sq;
        nearest[js] = jr;
      }
    }
  }

  return;
}

// ================================= random_initialize ==================================

unsigned long cSamp :: random_initialize ( int seed )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator,
//    it will behave as though it were seeded with value 1.
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRAND (the appropriate routine
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRAND, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
# define DEBUG 0

  unsigned long ul_seed;

  if ( seed != 0 )
  {
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
    }
  }
  else
  {
    seed = get_seed( );
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RAND with arbitrary SEED = " << seed << "\n";
    }
  }
//
//  Now set the seed.
//
  cout << "rodou aqui" << endl;
  ul_seed = ( unsigned long ) seed;
  srand ( ul_seed );

  return ul_seed;
# undef DEBUG
}

// ================================= get_seed ==================================

int cSamp :: get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I_MAX 2147483647
  time_t clock;
  int i;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  clock = time ( &tloc );
  lt = localtime ( &clock );
//
//  Hours is 1, 2, ..., 12.
//
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
//
//  Move Hours to 0, 1, ..., 11
//
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  isec = rand( )%60;

  seed = isec + 60 * ( imin + 60 * ihour );
//
//  We want values in [1,43200], not [0,43199].
//
  seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,IMAX].
//
  seed = ( int ) ( ( ( double ) seed )* ( ( double ) I_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
# undef I_MAX
}

// ================================= r8mat_uniform_01 ==================================

void cSamp :: r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  On output, SEED has
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}

// ================================= i4_to_halton_sequence ==================================

void cSamp :: i4_to_halton_sequence ( int dim_num, int n, int step, int seed[], int leap[],
  int base[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
//
//  Discussion:
//
//    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
//    sequences, each generated by a particular base.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Halton sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
//
//    of the original Halton sequence.
//
//
//    The data to be computed has two dimensions.
//
//    The number of data items is DIM_NUM * N, where DIM_NUM is the spatial dimension
//    of each element of the sequence, and N is the number of elements of the sequence.
//
//    The data is stored in a one dimensional array R.  The first element of the
//    sequence is stored in the first DIM_NUM entries of R, followed by the DIM_NUM entries
//    of the second element, and so on.
//
//    In particular, the J-th element of the sequence is stored in entries
//    0+(J-1)*DIM_NUM through (DIM_NUM-1) + (J-1)*DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
//
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of elements of the sequence.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required
//
//    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
//    to STEP = 0.
//
//    Input, int LEAP[DIM_NUM], the succesive jumps in the Halton sequence.
//
//    Input, int BASE[DIM_NUM], the Halton bases.
//
//    Output, double R[DIM_NUM*N], the next N elements of the
//    leaped Halton subsequence, beginning with element STEP.
//
{
  double base_inv;
  int digit;
  int i;
  int j;
  int *seed2;
//
//  Check the input.
//
  if ( !halham_dim_num_check ( dim_num ) )
  {
    exit ( 1 );
  }

  if ( !halham_n_check ( n ) )
  {
    exit ( 1 );
  }

  if ( !halham_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halham_seed_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !halham_leap_check ( dim_num, leap ) )
  {
    exit ( 1 );
  }

  if ( !halton_base_check ( dim_num, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  seed2 = new int[n];

  for ( i = 0; i < dim_num; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      seed2[j] = seed[i] + ( step + j ) * leap[i];
    }

    for ( j = 0; j < n; j++ )
    {
      r[i+j*dim_num] = 0.0;
    }

    for ( j = 0; j < n; j++ )
    {
      base_inv = 1.0 / ( ( double ) base[i] );

      while ( seed2[j] != 0 )
      {
        digit = seed2[j] % base[i];
        r[i+j*dim_num] = r[i+j*dim_num] + ( ( double ) digit ) * base_inv;
        base_inv = base_inv / ( ( double ) base[i] );
        seed2[j] = seed2[j] / base[i];
      }
    }
  }

  delete [] seed2;

  return;
}

// ================================= tuple_next_fast ==================================

void cSamp :: tuple_next_fast ( int m, int n, int rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between 1 and M.  The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2,
//    M = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank          X
//    ----          ----
//   -1            -1 -1
//
//    0             1  1
//    1             1  2
//    2             1  3
//    3             2  1
//    4             2  2
//    5             2  3
//    6             3  1
//    7             3  2
//    8             3  3
//    9             1  1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the maximum entry in each component.
//    M must be greater than 0.
//
//    Input, int N, the number of components.
//    N must be greater than 0.
//
//    Input, integer RANK, indicates the rank of the tuples.
//    Typically, 0 <= RANK < N**M; values larger than this are legal
//    and meaningful, and are equivalent to the corresponding value
//    MOD N**M.  If RANK < 0, this indicates that this is the first
//    call for the given values of (M,N).  Initialization is done,
//    and X is set to a dummy value.
//
//    Output, int X[N], the next tuple, or a dummy value if initialization
//    is being done.
//
{
  static int *base = NULL;
  int i;
//
  if ( rank < 0 )
  {
    if ( m <= 0 )
    {
      cout << "\n";
      cout << "TUPLE_NEXT_FAST - Fatal error!\n";
      cout << "  The value M <= 0 is not legal.\n";
      cout << "  M = " << m << "\n";
      exit ( 1 );
    }
    if ( n <= 0 )
    {
      cout << "\n";
      cout << "TUPLE_NEXT_FAST - Fatal error!\n";
      cout << "  The value N <= 0 is not legal.\n";
      cout << "  N = " << n << "\n";
      exit ( 1 );
    }

    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; 0 <= i; i-- )
    {
      base[i] = base[i+1] * m;
    }
    for ( i = 0; i < n; i++ )
    {
      x[i] = -1;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( rank / base[i] ) % m ) + 1;
    }
  }
  return;
}

// ================================= tuple_next_fast ==================================

void cSamp :: user ( int dim_num, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    USER samples points in a user-specified region with given density.
//
//  Discussion:
//
//    This routine can be used to
//
//    * specify an interesting initial configuration for the data,
//      by specifing that USER be used for initialization (INIT = 3);
//
//    * specify the shape of the computational region, by specifying
//      that sample points are to be generated by this routine,
//      (SAMPLE = 3) and then returning sample points uniformly at random.
//
//    * specify the distribution or density function, by specifying
//      that sample points are to be generated by this routine,
//      (SAMPLE = 3 ) and then returning sample points according to a
//      given probability density function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension.
//
//    Input, integer N, the number of sample points desired.
//
//    Input/output, int *SEED, the "seed" value.  On output, SEED has
//    been updated.
//
//    Output, double R[DIM_NUM*N], the sample values.
//
{
# define PI 3.141592653589793

  double angle;
  int j;
  double radius;

  for ( j = 0; j < n; j++ )
  {
    angle = 2.0 * PI * (double) rand( ) / (double) RAND_MAX;
    radius = sqrt ( (double) rand( ) / (double) RAND_MAX );
    r[0+j*2] = radius*cos(angle);
    r[1+j*2] = radius*sin(angle);
  }

  return;
# undef PI
}

// ================================= r8_huge ==================================

double cSamp :: r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8.
//
{
  return HUGE_VAL;
}

// ================================= halham_dim_num_check ==================================

bool cSamp :: halham_dim_num_check ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//    DIM_NUM must be positive.
//
//    Output, bool HALHAM_DIM_NUM_CHECK, is true if DIM_NUM is legal.
//
{
  bool value;

  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "HALHAM_DIM_NUM_CHECK - Fatal error!\n";
    cout << "  DIM_NUM < 0.";
    cout << "  DIM_NUM = " << dim_num << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}

// ================================= halham_n_check ==================================

bool cSamp :: halham_n_check ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of the subsequence.
//    N must be positive.
//
//    Output, bool HALHAM_N_CHECK, is true if N is legal.
//
{
  bool value;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "HALHAM_N_CHECK - Fatal error!\n";
    cout << "  N < 0.";
    cout << "  N = " << n << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}

// ================================= halham_step_check ==================================

bool cSamp :: halham_step_check ( int step )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int STEP, the index of the subsequence element.
//    STEP must be 1 or greater.
//
//    Output, bool HALHAM_STEP_CHECK, is true if STEP is legal.
//
{
  int i;
  bool value;

  if ( step < 0 )
  {
    cout << "\n";
    cout << "HALHAM_STEP_CHECK - Fatal error!\n";
    cout << "  STEP < 0.";
    cout << "  STEP = " << step << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}

// ================================= halham_seed_check ==================================

bool cSamp :: halham_seed_check ( int dim_num, int seed[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int SEED[DIM_NUM], the sequence index
//    corresponding to STEP = 0.  Each entry must be 0 or greater.
//
//    Output, bool HALHAM_SEED_CHECK, is true if SEED is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( seed[i] < 0 )
    {
      cout << "\n";
      cout << "HALHAM_SEED_CHECK - Fatal error!\n";
      cout << "  SEED entries must be nonnegative.\n";
      cout << "  seed[" << i << "] = " << seed[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}

// ================================= halham_leap_check ==================================

bool cSamp :: halham_leap_check ( int dim_num, int leap[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEAP[DIM_NUM], the successive jumps in the sequence.
//    Each entry must be greater than 0.
//
//    Output, bool HALHAM_LEAP_CHECK, is true if LEAP is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( leap[i] < 1 )
    {
      cout << "\n";
      cout << "HALHAM_LEAP_CHECK - Fatal error!\n";
      cout << "  Leap entries must be greater than 0.\n";
      cout << "  leap[" << i << "] = " << leap[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}

// ================================= halham_base_check ==================================

bool cSamp :: halton_base_check ( int dim_num, int base[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_BASE_CHECK checks BASE for a Halton sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases.
//
//    Output, bool HALTON_BASE_CHECK, is true if BASE is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( base[i] <= 1 )
    {
      cout << "\n";
      cout << "HALTON_BASE_CHECK - Fatal error!\n";
      cout << "  Bases must be greater than 1.\n";
      cout << "  base[" << i << "] = " << base[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}

// ================================= CalcSampLCVT ==================================

void cSamp :: CalcSampLCVT( int dim_num, int n, int batch, int init, int sample, int sample_num, int it_max, int it_fixed, int *seed, double r[],int *it_num, double *it_diff, double *energy, vector<cVector> &sx  )

//****************************************************************************80
//
//  Purpose:
//
//    CVT computes a Centroidal Voronoi Tessellation.
//
//  Discussion:
//
//    This routine initializes the data, and carries out the
//    CVT iteration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int INIT, specifies how the points are to be initialized.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine;
//     4, points are already initialized on input.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Input, int IT_FIXED, the maximum number of iterations to take
//    with a fixed set of sample points.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the approximate CVT points.
//    If INIT = 4 on input, then it is assumed that these values have been
//    initialized.  On output, the CVT iteration has been applied to improve
//    the value of the points.
//
//    Output, int *IT_NUM, the number of iterations taken.  Generally,
//    this will be equal to IT_MAX, unless the iteration tolerance was
//    satisfied early.
//
//    Output, double *IT_DIFF, the L2 norm of the difference
//    between the iterates.
//
//    Output, double *ENERGY,  the discrete "energy", divided
//    by the number of sample points.
//
{
  bool DEBUG = false;
  int i;
  bool initialize;
  int seed_base;
  int seed_init;
  cVector xs(dim_num);

  if ( batch < 1 )
  {
    cout << "\n";
    cout << "CVT - Fatal error!\n";
    cout << "  The input value BATCH < 1.\n";
    exit ( 1 );
  }

  if ( *seed <= 0 )
  {
    cout << "\n";
    cout << "CVT - Fatal error!\n";
    cout << "  The input value SEED <= 0.\n";
    exit ( 1 );
  }

  if ( DEBUG )
  {
    cout << "\n";
    cout << "  Step       SEED          L2-Change        Energy\n";
    cout << "\n";
  }

  *it_num = 0;
  *it_diff = 0.0;
  *energy = 0.0;
  seed_init = *seed;
//
//  Initialize the data, unless the user has already done that.
//
  if ( init != 4 )
  {
    initialize = true;
    cvt_sample ( dim_num, n, n, init, initialize, seed, r );
  }
  if ( DEBUG )
  {
    cout                          << "  "
         << setw(4)  << *it_num   << "  "
         << setw(12) << seed_init << "\n";
  }
//
//  If the initialization and sampling steps use the same random number
//  scheme, then the sampling scheme does not have to be initialized.
//
  if ( init == sample )
  {
    initialize = false;
  }
  else
  {
    initialize = true;
  }
//
//  Carry out the iteration.
//
  while ( *it_num < it_max )
  {
//
//  If it's time to update the seed, save its current value
//  as the starting value for all iterations in this cycle.
//  If it's not time to update the seed, restore it to its initial
//  value for this cycle.
//
    if ( ( (*it_num) % it_fixed ) == 0 )
    {
      seed_base = *seed;
    }
    else
    {
      *seed = seed_base;
    }

    *it_num = *it_num + 1;
    seed_init = *seed;

    cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, seed,
      r, it_diff, energy );

    initialize = false;

    if ( DEBUG )
    {
      cout                          << "  "
           << setw(4)  << *it_num   << "  "
           << setw(12) << seed_init << "  "
           << setw(14) << *it_diff  << "  "
           << setw(14) << *energy    << "\n";
    }
  }
//
//  If the initialization and sampling steps use the same random number
//  scheme, then the sampling scheme does not have to be initialized.
//
  if ( init == sample )
  {
    initialize = false;
  }
  else
  {
    initialize = true;
  }
//
//  Carry out the iteration.
//
  while ( *it_num < it_max )
  {
//
//  If it's time to update the seed, save its current value
//  as the starting value for all iterations in this cycle.
//  If it's not time to update the seed, restore it to its initial
//  value for this cycle.
//
    if ( ( (*it_num) % it_fixed ) == 0 )
    {
      seed_base = *seed;
    }
    else
    {
      *seed = seed_base;
    }

    *it_num = *it_num + 1;
    seed_init = *seed;

    cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, seed, r, it_diff, energy );

    initialize = false;

    if ( DEBUG )
    {
      cout                          << "  "
           << setw(4)  << *it_num   << "  "
           << setw(12) << seed_init << "  "
           << setw(14) << *it_diff  << "  "
           << setw(14) << *energy    << "\n";
    }
  }


  r8mat_latinize(dim_num, n, r);

  for ( int j = 0; j < n; j++)
  {
  for ( int i = 0; i < dim_num; i++ )
  {
    xs[i] = r[i+j*dim_num];
  }
  sx.push_back(xs);
  }

  return;
}

// ================================= r8mat_latinize ==================================

void cSamp :: r8mat_latinize ( int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_LATINIZE "Latinizes" an R8MAT.
//
//  Discussion:
//
//    It is assumed, though not necessary, that the input dataset
//    has points that lie in the unit hypercube.
//
//    In any case, the output dataset will have this property.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of cells.
//
//    Input/output, double TABLE[M*N].  On input, the dataset to
//    be "Latinized".  On output, the Latinized dataset.
//
{
  double *column;
  int i;
  int *indx;
  int j;

  column = new double[n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      column[j] = table[i+j*m];
    }
    indx = r8vec_sort_heap_index_a ( n, column );

    for ( j = 0; j < n; j++ )
    {
      table[i+indx[j]*m] = ( double ) ( 2 * j + 1 ) / ( double ) ( 2 * n );
    }

    delete [] indx;
  }

  delete [] column;

  return;
}

// ================================= r8vec_sort_heap_index_a ==================================

int* cSamp :: r8vec_sort_heap_index_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(INDX(I)), I = 1 to N is sorted,
//
//    after which A(I), I = 1 to N is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], an array to be index-sorted.
//
//    Output, int R8VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  indx = new int[n];

  for ( i = 1; i <= n; i++ )
  {
    indx[i-1] = i;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        for ( i = 0; i < n; i++ )
        {
          indx[i] = indx[i] - 1;
        }
        return indx;
      }

    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]-1] < a[indx[j]-1] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]-1] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
}

// ======================================================= End of file =====


