//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file RandomGenerator.h
 * @brief Declare a global Random Number Generator
 *
 * Selected among
 * - boost::random
 * - sprng
 * - math::random
 * qmcplusplus::Random() returns a random number [0,1)
 * For OpenMP is enabled, it is important to use thread-safe boost::random. Each
 * thread uses its own random number generator with a distinct seed. This prevents
 * a use of global lock which will slow down the applications significantly.
 */
#ifndef OHMMS_RANDOMGENERATOR
#define OHMMS_RANDOMGENERATOR
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cmath>
#include <ctime>

#include <stdint.h>

inline uint32_t make_seed(int i, int n)
{
  return static_cast<uint32_t>(std::time(0))%10474949+(i+1)*n+i;
}

// The definition of the fake RNG should always be available for unit testing
#include "Utilities/FakeRandom.h"
#ifdef USE_FAKE_RNG
namespace qmcplusplus
{
typedef FakeRandom RandomGenerator_t;
extern RandomGenerator_t Random;
}
#else

#ifdef HAVE_LIBBOOST

#include "Utilities/BoostRandom.h"
namespace qmcplusplus
{
typedef BoostRandom<OHMMS_PRECISION_FULL> RandomGenerator_t;
extern RandomGenerator_t Random;
}
#else

#ifdef USE_SPRNG
#include "Utilities/SprngRandom.h"
namespace qmcplusplus
{
typedef SprngRandom<0> RandomGenerator_t;
extern RandomGenerator_t Random;
}
#else
#include "Utilities/SimpleRandom.h"
namespace qmcplusplus
{
typedef SimpleRandom<MTRand> RandomGenerator_t;
extern RandomGenerator_t Random;
}
#endif
#endif
#endif
#endif

