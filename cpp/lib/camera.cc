//
// Created by daeyun on 12/20/17.
//

#include "camera.h"

//#include "file_io.h"

namespace mesh_to_depth {

FrustumParams MakePerspectiveFrustumParams(double hw_ratio, double x_fov, double near, double far) {
  FrustumParams ret;
  ret.right = std::abs(std::tan(x_fov) * near);
  ret.top = ret.right * hw_ratio;

  ret.left = -ret.right;
  ret.bottom = -ret.top;

  ret.near = near;
  ret.far = far;
  return ret;
}

FrustumParams ForceFixedAspectRatio(double hw_ratio, const FrustumParams &frustum) {
  double left = frustum.left;
  double right = frustum.right;
  double top = frustum.top;
  double bottom = frustum.bottom;

  double lr = std::abs(right - left);
  double bt = std::abs(top - bottom);
  double box_hw_ratio = bt / lr;

  if (box_hw_ratio < hw_ratio) {
    // image is squeezed horizontally.
    double padding = (hw_ratio * lr - bt) * 0.5;
    top += padding;
    bottom -= padding;
  } else {
    double padding = (bt - hw_ratio * lr) / hw_ratio * 0.5;
    right += padding;
    left -= padding;
  }

  FrustumParams ret = frustum;
  ret.left = left;
  ret.right = right;
  ret.top = top;
  ret.bottom = bottom;

  return ret;
}
}
