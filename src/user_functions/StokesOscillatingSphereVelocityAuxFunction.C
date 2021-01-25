// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//



#include <user_functions/StokesOscillatingSphereVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>

namespace sierra{
namespace nalu{

StokesOscillatingSphereVelocityAuxFunction::StokesOscillatingSphereVelocityAuxFunction(
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
  Q_   = -6.*pi_*std::pow(r0_,3.)*(exp(lam_)-1.-lam_-std::pow(lam_,2.)/3.)/lam_/lam_;
}

void
StokesOscillatingSphereVelocityAuxFunction::do_evaluate(
  const double *coords,
  const double t,
  const unsigned /*spatialDimension*/,
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
  std::complex<double> BCoeff = std::complex<double>{0.,1.} *
                                exp(std::complex<double>{0.,-1.}*omg_*t)*B_/(8.*pi_*visc_);
  std::complex<double> QCoeff = std::complex<double>{0.,1.} *
                                exp(std::complex<double>{0.,-1.}*omg_*t)*Q_/(4.*pi_);

  for(unsigned p=0; p < numPoints; ++p) {
    // get realtive position of point w.r.t. sphere center
    double x = coords[0] - cenX;
    double y = coords[1] - cenY;
    double z = coords[2] - cenZ;
    double dist = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
    
    // solution does not exist inside sphere
    if(dist < r0_-eps_) {
      fieldPtr[0] = 0.0;
      fieldPtr[1] = 0.0;
      fieldPtr[2] = 0.0;
      
      fieldPtr += fieldSize;
      coords += fieldSize;
      continue;
    }

    // dimensionless distance metrics
    std::complex<double> R = lam_*dist/r0_;
    std::complex<double> R2 = std::pow(R,2.);
    std::complex<double> expR = exp(-R);

    double vDotx = vCoeff_[0]*x+vCoeff_[1]*y+vCoeff_[2]*z;

    // different components of 1st term
    std::vector<std::complex<double>> B1 {(2.*expR*(1.+1./R+1./R2) - 2./R2)*vCoeff_[0]/(dist+eps_),
                                          (2.*expR*(1.+1./R+1./R2) - 2./R2)*vCoeff_[1]/(dist+eps_),
                                          (2.*expR*(1.+1./R+1./R2) - 2./R2)*vCoeff_[2]/(dist+eps_)};
    std::vector<std::complex<double>> B2 {(6./R2 - 2.*expR*(1.+3./R+3./R2))*vDotx*x/(std::pow(dist,3.)+eps_),
                                          (6./R2 - 2.*expR*(1.+3./R+3./R2))*vDotx*y/(std::pow(dist,3.)+eps_),
                                          (6./R2 - 2.*expR*(1.+3./R+3./R2))*vDotx*z/(std::pow(dist,3.)+eps_)};
    std::vector<std::complex<double>> T1 {BCoeff*(B1[0] + B2[0]),
                                          BCoeff*(B1[1] + B2[1]),
                                          BCoeff*(B1[2] + B2[2])};

    // different components of 2nd term
    std::vector<std::complex<double>> Q1 {-expR*(1.+R+R2)*vCoeff_[0]/(std::pow(dist,3.)+eps_),
                                          -expR*(1.+R+R2)*vCoeff_[1]/(std::pow(dist,3.)+eps_),
                                          -expR*(1.+R+R2)*vCoeff_[2]/(std::pow(dist,3.)+eps_)};
    std::vector<std::complex<double>> Q2 {3.*expR*(1.+R+R2/3.)*vDotx*x/(std::pow(dist,5.)+eps_),
                                          3.*expR*(1.+R+R2/3.)*vDotx*y/(std::pow(dist,5.)+eps_),
                                          3.*expR*(1.+R+R2/3.)*vDotx*z/(std::pow(dist,5.)+eps_)};
    std::vector<std::complex<double>> T2 {QCoeff*(Q1[0] + Q2[0]),
                                          QCoeff*(Q1[1] + Q2[1]),
                                          QCoeff*(Q1[2] + Q2[2])};

    // extract real parts
    fieldPtr[0] = std::real(T1[0]) + std::real(T2[0]);
    fieldPtr[1] = std::real(T1[1]) + std::real(T2[1]);
    fieldPtr[2] = std::real(T1[2]) + std::real(T2[2]);

    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
