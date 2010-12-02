// This file is derived from the code  at 
// http://www.jcornwall.me.uk/2009/04/mersenne-twisters-in-cuda/
// which in turn is derived from the NVIDIA CUDA SDK example 'MersenneTwister'.

/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

// Some parts contain code from Makoto Matsumoto and Takuji Nishimura's dci.h

/* Copyright (C) 2001-2006 Makoto Matsumoto and Takuji Nishimura.  */
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */

//#include <cassert>
//#include <cstdio>
//#include <vector>

#define MT_MM     9
#define MT_NN     19
#define MT_WMASK  0xFFFFFFFFU
#define MT_UMASK  0xFFFFFFFEU
#define MT_LMASK  0x1U
#define MT_RNG_COUNT 32768
#define MT_SHIFT0 12
#define MT_SHIFTB 7
#define MT_SHIFTC 15
#define MT_SHIFT1 18

// Record format for MersenneTwister.dat, created by spawnTwisters.c
// size = 16 bytes
struct __align__(16) mt_struct_stripped {
  unsigned int matrix_a;
  unsigned int mask_b;
  unsigned int mask_c;
  unsigned int seed;
};

// Per-thread state object for a single twister.
// size = 84 bytes but aligned = 96 bytes
struct __align__(16) MersenneTwisterState {
  unsigned int mt[MT_NN];
  int iState;
  unsigned int mti1;
};

// Preloaded, offline-generated seed data structure.
__device__ static mt_struct_stripped MT[MT_RNG_COUNT];

// Hold the current states of the twisters
__device__ static MersenneTwisterState MTS[MT_RNG_COUNT];


__device__ void MersenneTwisterInitialise(MersenneTwisterState *state, unsigned int threadID) {
  state->mt[0] = MT[threadID].seed;
  for(int i = 1; i < MT_NN; ++ i) {
    state->mt[i] = (1812433253U * (state->mt[i - 1] ^ (state->mt[i - 1] >> 30)) + i) & MT_WMASK;
  }
  
  state->iState = 0;
  state->mti1 = state->mt[0];
}

__device__ unsigned int MersenneTwisterGenerate(MersenneTwisterState *state, unsigned int threadID) {
  int iState1 = state->iState + 1;
  int iStateM = state->iState + MT_MM;

  if(iState1 >= MT_NN) iState1 -= MT_NN;
  if(iStateM >= MT_NN) iStateM -= MT_NN;

  unsigned int mti = state->mti1;
  state->mti1 = state->mt[iState1];
  unsigned int mtiM = state->mt[iStateM];

  unsigned int x = (mti & MT_UMASK) | (state->mti1 & MT_LMASK);
  x = mtiM ^ (x >> 1) ^ ((x & 1) ? MT[threadID].matrix_a : 0);
  state->mt[state->iState] = x;
  state->iState = iState1;

  // Tempering transformation.
  x ^= (x >> MT_SHIFT0);
  x ^= (x << MT_SHIFTB) & MT[threadID].mask_b;
  x ^= (x << MT_SHIFTC) & MT[threadID].mask_c;
  x ^= (x >> MT_SHIFT1);
  
  return x;
}

__global__ void InitialiseAllMersenneTwisters() {

  for(int i=0; i<MT_RNG_COUNT; ++i){
    MersenneTwisterInitialise(&(MTS[i]),i);
  }
}

