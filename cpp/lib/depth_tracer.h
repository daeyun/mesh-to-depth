//
// Created by Daeyun Shin on 12/21/20.
//

#pragma once

#include "ray_mesh_intersection.h"
#include "camera.h"

namespace mesh_to_depth {

class DepthTracer {
 public:

  // `ray_tracer` is managed externally.
  DepthTracer(const RayTracer *ray_tracer,
              const Camera *camera,
              size_t width,
              size_t height)
      : ray_tracer_(ray_tracer), camera_(camera), width_(width), height_(height) {
    if (camera_->is_perspective()) {
      double xf, yf;
      // This whole block of code is from an older version. It makes sure aspect stretching does not happen. There's probably a better way to do this.
      camera->fov(&xf, &yf);
      // Distance to the image plane according to the x fov.
      double xl = 0.5 * width_ / std::tan(xf);
      // Distance to the image plane according to the y fov.
      double yl = 0.5 * height_ / std::tan(yf);

      // For now, we assume the aspect ratio is always 1.0. So the distance to image plane should end up being the same according to both x and y.
      // Otherwise the image size or focal length is wrong. This can also happen because of precision error.
      // 0.01 is an arbitrary threshold.
      if (std::abs(xl - yl) > 0.01) {
        LOGGER->warn("xf: {}, yf: {}, width: {}, height: {}, xl: {}, yl: {}", xf, yf, width_, height_, xl, yl);
        throw std::runtime_error("Inconsistent distance to image plane.");
      }
      // Compute the average of the two distances. There are probably other, better ways to do this.
      image_focal_length_ = (xl + yl) * 0.5;
    }

    image_optical_center_ = Vec3{width_ * 0.5, height_ * 0.5, 0};

  }

  const RayTracer *ray_tracer() const {
    return ray_tracer_;
  }

  Vec3 RayOrigin(int x, int y) const {
    if (camera_->is_perspective()) {
      return camera_->position();
    } else {
      double im_x = static_cast<double>(x) + 0.5;
      double im_y = height_ - (static_cast<double>(y) + 0.5);
      const auto &frustum = camera_->frustum();
      double cam_x = im_x / width_ * (frustum.right - frustum.left) + frustum.left;
      double cam_y = im_y / height_ * (frustum.top - frustum.bottom) + frustum.bottom;
      Vec3 ret;
      camera_->CamToWorld(Vec3{cam_x, cam_y, 0}, &ret);
      return ret;
    }
  }

  Vec3 RayDirection(int x, int y) const {
    Vec3 ray_direction;
    if (camera_->is_perspective()) {
      Vec3 image_plane_coord{static_cast<double>(x) + 0.5, height_ - (static_cast<double>(y) + 0.5), -image_focal_length()};
      Vec3 cam_ray_direction = (image_plane_coord - image_optical_center()).normalized();
      camera_->CamToWorldNormal(cam_ray_direction, &ray_direction);
    } else {
      ray_direction = camera_->viewing_direction();
    }
    ray_direction.normalize();
    return ray_direction;
  }

  // Implementation specific. If `is_depth` is false, ray displacement (t) value will be assigned to `out_depth` instead.
  virtual bool DepthValue(int x, int y, bool is_depth, float *out_depth, uint32_t *out_prim_id) const = 0;
//  virtual int DepthValues(int x, int y, vector<float> *out, vector<uint32_t> *prim_ids) const = 0;
//  virtual int ObjectCenteredRayDisplacement(int x, int y, vector<float> *out, vector<uint32_t> *prim_ids) const = 0;

  size_t width() const {
    return width_;
  }
  size_t height() const {
    return height_;
  }

 protected:
  const RayTracer *ray_tracer_;
  const Camera *camera_;  // Managed externally
  size_t width_;
  size_t height_;

 private:
  double image_focal_length() const {
    if (camera_->is_perspective()) {
      return image_focal_length_;
    } else {
      throw std::runtime_error("image focal length should not be used in orthographic projection.");
    }
  }

  const Vec3 &image_optical_center() const {
    if (camera_->is_perspective()) {
      return image_optical_center_;
    } else {
      throw std::runtime_error("image optical center should not be used in orthographic projection.");
    }
  }

  double image_focal_length_;
  Vec3 image_optical_center_;
};

class SimpleMultiLayerDepthRenderer : public DepthTracer {
 public:
  SimpleMultiLayerDepthRenderer(const RayTracer *ray_tracer,
                                const Camera *camera,
                                size_t width,
                                size_t height)
      : DepthTracer(ray_tracer, camera, width, height) {}

  virtual bool DepthValue(int x, int y, bool is_depth, float *out_value, uint32_t *out_prim_id) const override {
    const Vec3 ray_direction = this->RayDirection(x, y);
    bool is_hit = false;

    // Depth values are collected in the callback function, in the order traversed.
    ray_tracer_->Traverse(this->RayOrigin(x, y), ray_direction, [&](float t, float u, float v, unsigned int prim_id) -> bool {
      if (t > camera_->frustum().far) {
        return false;
      }
      // TODO(daeyun): this is not efficient. better to move the ray origin.
      if (t < camera_->frustum().near) {
        return true;
      }

      *out_value = t;
      *out_prim_id = prim_id;
      is_hit = true;
      return false;  // Stop at first hit.
    });

    if (is_depth) {
      // Convert ray displacement to depth.
      const float z = camera_->viewing_direction().dot(ray_direction);
      *out_value *= z;
    }

    return is_hit;
  }

/* TODO(daeyun): WIP
  virtual int DepthValues(int x, int y, vector<float> *out_values, vector<uint32_t> *prim_ids) const override {
    Vec3 ray_direction = this->RayDirection(x, y);

    int count = 0;

    // Depth values are collected in the callback function, in the order traversed.
    ray_tracer_->Traverse(this->RayOrigin(x, y), ray_direction, [&](float t, float u, float v, unsigned int prim_id) -> bool {
      out_values->push_back(t);
      prim_ids->push_back(prim_id);
      ++count;
    });

    // Convert ray displacement to depth.
    const double z = camera_->viewing_direction().dot(ray_direction);
    for (auto &t : *out_values) {
      t *= z;
    }

    return count;
  }
*/

};
}
