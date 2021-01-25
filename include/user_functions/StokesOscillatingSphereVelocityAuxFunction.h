// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef StokesOscillatingSphereVelocityAuxFunction_h
#define StokesOscillatingSphereVelocityAuxFunction_h

#include <AuxFunction.h>

#include <complex>
#include <vector>

namespace sierra{
namespace nalu{

class StokesOscillatingSphereVelocityAuxFunction : public AuxFunction
{
public:

  StokesOscillatingSphereVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos);

  virtual ~StokesOscillatingSphereVelocityAuxFunction() {}
  
  using AuxFunction::do_evaluate;
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  double omg_{0.0};
  double pi_{0.0};
  double r0_{0.0};
  double visc_{0.0};
  double eps_{1e-8};

  std::complex<double> lam_{0.0,0.0};
  std::complex<double> B_{0.0,0.0};
  std::complex<double> Q_{0.0,0.0};

  std::vector<double> cen_{0.0,0.0,0.0};
  std::vector<double> vCoeff_{0.0,0.0,0.0};
};

} // namespace nalu
} // namespace Sierra

#endif
