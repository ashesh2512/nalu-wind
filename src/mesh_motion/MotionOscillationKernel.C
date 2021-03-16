#include "mesh_motion/MotionOscillationKernel.h"

#include<FieldTypeDef.h>
#include <NaluParsing.h>

// stk_mesh/base/fem
#include <stk_mesh/base/FieldBLAS.hpp>

namespace sierra{
namespace nalu{

MotionOscillationKernel::MotionOscillationKernel(
  stk::mesh::MetaData& meta,
  const YAML::Node& node)
  : NgpMotionKernel<MotionOscillationKernel>()
{
  load(node);
}

void MotionOscillationKernel::load(const YAML::Node& node)
{
  // perturb start and end times with a small value for
  // accurate comparison with floats
  get_if_present(node, "start_time", startTime_, startTime_);
  startTime_ = startTime_-DBL_EPSILON;

  get_if_present(node, "end_time", endTime_, endTime_);
  endTime_ = endTime_+DBL_EPSILON;

  // get origin based on if it was defined
  if( node["centroid"] ) {
    for (int d=0; d < nalu_ngp::NDimMax; ++d)
      origin_[d] = node["centroid"][d].as<double>();
  }

  // get frequency if it was defined
  if( node["frequency"] ) {
    for (int d=0; d < nalu_ngp::NDimMax; ++d)
      frequency_[d] = node["frequency"][d].as<double>();
  }

  // get amplitude if it was defined
  if( node["amplitude"] ) {
    for (int d=0; d < nalu_ngp::NDimMax; ++d)
      amplitude_[d] = node["amplitude"][d].as<double>();
  }
}

mm::TransMatType MotionOscillationKernel::build_transformation(
  const double& time,
  const mm::ThreeDVecType& /* xyz */)
{
  mm::TransMatType transMat;

  if(time < (startTime_)) return transMat;
  double currTime = (time < endTime_)? time : endTime_;

  // add appropriate translation
  transMat[0*mm::matSize+3] = origin_[0] - amplitude_[0] * stk::math::cos(frequency_[0]*currTime);
  transMat[1*mm::matSize+3] = origin_[1] - amplitude_[1] * stk::math::cos(frequency_[1]*currTime);
  transMat[2*mm::matSize+3] = origin_[2] - amplitude_[2] * stk::math::cos(frequency_[2]*currTime);

  // composite addition of motions
  return transMat;
}

mm::ThreeDVecType MotionOscillationKernel::compute_velocity(
  const double& time,
  const mm::TransMatType&  /* compTrans */,
  const mm::ThreeDVecType& /* mxyz */,
  const mm::ThreeDVecType& /* cxyz */)
{
  mm::ThreeDVecType vel;

  if((time < startTime_) || (time > endTime_))
    return vel;

  for (int d=0; d < nalu_ngp::NDimMax; d++)
    vel[d] = amplitude_[d]*frequency_[d] * stk::math::sin(frequency_[d]*time);

  return vel;
}

} // nalu
} // sierra
