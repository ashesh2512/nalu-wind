// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//



#include <user_functions/StokesOscillatingSpherePressureAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>

namespace sierra{
namespace nalu{

StokesOscillatingSpherePressureAuxFunction::StokesOscillatingSpherePressureAuxFunction() :
  AuxFunction(0,1),
  pi_(std::acos(-1.0)),
  r0_(0.5),
  visc_(1.0)
{
  vCoeff_[0] = 1.0;

  omg_ = 2*pi_;
  lam_ = std::complex<double>{1.,-1.}*r0_*std::sqrt(0.5*std::abs(omg_)/visc_);
  B_   = 6.*pi_*visc_*r0_*(1.+lam_+std::pow(lam_,2.)/3.);
}


void
StokesOscillatingSpherePressureAuxFunction::do_evaluate(
  const double *coords,
  const double t,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  // center as a function of time
  double cenX = cen_[0]-std::real(exp(std::complex<double>{0.,-1.}*omg_*t)/omg_*vCoeff_[0]);
  double cenY = cen_[1]-std::real(exp(std::complex<double>{0.,-1.}*omg_*t)/omg_*vCoeff_[1]);
  double cenZ = cen_[2]-std::real(exp(std::complex<double>{0.,-1.}*omg_*t)/omg_*vCoeff_[2]);

  // premultipliers
  std::complex<double> PCoeff = std::complex<double>{0.,1.} *
                                exp(std::complex<double>{0.,-1.}*omg_*t)*B_/(4.*pi_);

  for(unsigned p=0; p < numPoints; ++p) {
    // get realtive position of point w.r.t. sphere center
    double x = coords[0] - cenX;
    double y = coords[1] - cenY;
    double z = coords[2] - cenZ;
    double dist = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));

    double vDotx = vCoeff_[0]*x+vCoeff_[1]*y+vCoeff_[2]*z;

    fieldPtr[0] = std::real(PCoeff*vDotx/std::pow(dist,3.));

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

StokesOscillatingSpherePressureGradAuxFunction::StokesOscillatingSpherePressureGradAuxFunction(
  const unsigned beginPos,
  const unsigned endPos) :
  AuxFunction(beginPos, endPos),
  pi_(std::acos(-1.0)),
  r0_(0.5),
  visc_(1.0)
{
  vCoeff_[0] = 1.0;

  omg_ = 2*pi_;
  lam_ = std::complex<double>{1.,-1.}*r0_*std::sqrt(0.5*std::abs(omg_)/visc_);
  B_   = 6.*pi_*visc_*r0_*(1.+lam_+std::pow(lam_,2.)/3.);
}


void
StokesOscillatingSpherePressureGradAuxFunction::do_evaluate(
  const double *coords,
  const double t,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  // center as a function of time
  double cenX = cen_[0]-std::real(exp(std::complex<double>{0.,-1.}*omg_*t)/omg_*vCoeff_[0]);
  double cenY = cen_[1]-std::real(exp(std::complex<double>{0.,-1.}*omg_*t)/omg_*vCoeff_[1]);
  double cenZ = cen_[2]-std::real(exp(std::complex<double>{0.,-1.}*omg_*t)/omg_*vCoeff_[2]);

  // premultipliers
  std::complex<double> PCoeff = std::complex<double>{0.,1.} *
                                exp(std::complex<double>{0.,-1.}*omg_*t)*B_/(4.*pi_);

  for(unsigned p=0; p < numPoints; ++p) {
    // get realtive position of point w.r.t. sphere center
    double x = coords[0] - cenX;
    double y = coords[1] - cenY;
    double z = coords[2] - cenZ;
    double dist = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));

    double vDotx = vCoeff_[0]*x+vCoeff_[1]*y+vCoeff_[2]*z;

    fieldPtr[0] =  std::real(PCoeff*(vCoeff_[0]/std::pow(dist,3.) - 3.*vDotx*x/std::pow(dist,5.)));
    fieldPtr[1] =  std::real(PCoeff*(vCoeff_[1]/std::pow(dist,3.) - 3.*vDotx*y/std::pow(dist,5.)));
    fieldPtr[2] =  std::real(PCoeff*(vCoeff_[2]/std::pow(dist,3.) - 3.*vDotx*z/std::pow(dist,5.)));

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
