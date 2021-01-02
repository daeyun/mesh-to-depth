//
// Created by daeyun on 4/11/17.
//

#pragma once

#include <fstream>
#include <iomanip>

#include "common.h"

namespace mesh_to_depth {
struct FrustumParams {
  double left = -1;
  double right = 1;
  double bottom = -1;
  double top = 1;
  double near = 0.001;
  double far = 50;
};

// `hw_ratio` is height/width. e.g. 0.75
// `x_fov` must be half-angle, in radians.
FrustumParams MakePerspectiveFrustumParams(double hw_ratio, double x_fov, double near, double far);

struct Plane {
  Vec3 ref_point;
  Vec3 up_normal;

  double Displacement(const Vec3 &point) const {
    return up_normal.dot(point - ref_point);
  }

  double Displacement(const array<float, 3> &point) const {
    return up_normal.dot(Vec3{point[0], point[1], point[2]} - ref_point);
  }

  bool IsAbove(const Vec3 &point) const {
    return Displacement(point) >= 0;
  }

  bool IsAbove(const array<float, 3> &point) const {
    return Displacement(point) >= 0;
  }

  bool IntersectRay(const Vec3 &origin, const Vec3 &dir, double *t) const {
    // Assuming vectors are all normalized.
    double denom = up_normal.dot(dir);
    if (std::abs(denom) > 1e-6) {
      *t = up_normal.dot(ref_point - origin) / denom;
      return *t >= 0;
    }
    return false;
  }

  // Used for visualization.
  bool ToTriangle(float size, vector<array<unsigned int, 3>> *faces, vector<array<float, 3>> *vertices) const {
    Vec3 a = up_normal.cross(Vec3::Random()).normalized();
    Vec3 b = up_normal.cross(a).normalized();
    Vec3 v1 = ref_point + a * size;
    Vec3 v2 = ref_point + b * size;
    unsigned int offset = vertices->size();
    vertices->push_back(array < float, 3 > {static_cast<float>(ref_point[0]), static_cast<float>(ref_point[1]), static_cast<float>(ref_point[2])});
    vertices->push_back(array < float, 3 > {static_cast<float>(v1[0]), static_cast<float>(v1[1]), static_cast<float>(v1[2])});
    vertices->push_back(array < float, 3 > {static_cast<float>(v2[0]), static_cast<float>(v2[1]), static_cast<float>(v2[2])});
    faces->push_back(array < unsigned int, 3 > {offset, offset + 1, offset + 2});
  }
};

class Camera {
 public:
  Camera(const Vec3 &camera_position,
         const Vec3 &lookat_position,
         const Vec3 &up,
         const FrustumParams &frustum)
      : position_(camera_position),
        lookat_position_(lookat_position),
        frustum_(frustum) {
    Vec3 viewing_direction = (lookat_position - camera_position).normalized();
    Vec3 right = viewing_direction.cross(up).normalized();
    Vec3 up_vector = right.cross(viewing_direction).normalized();

    viewing_direction_ = viewing_direction;
    up_ = up_vector;  // Up vector can change if the initial value is not orthogonal.

    view_mat_(0, 0) = right[0];
    view_mat_(0, 1) = right[1];
    view_mat_(0, 2) = right[2];
    view_mat_(1, 0) = up_vector[0];
    view_mat_(1, 1) = up_vector[1];
    view_mat_(1, 2) = up_vector[2];
    view_mat_(2, 0) = -viewing_direction[0];
    view_mat_(2, 1) = -viewing_direction[1];
    view_mat_(2, 2) = -viewing_direction[2];

    Vec3 translation = -(view_mat_.topLeftCorner<3, 3>() * camera_position);
    view_mat_(0, 3) = translation[0];
    view_mat_(1, 3) = translation[1];
    view_mat_(2, 3) = translation[2];

    view_mat_(3, 0) = 0;
    view_mat_(3, 1) = 0;
    view_mat_(3, 2) = 0;
    view_mat_(3, 3) = 1;

    view_mat_inv_.topLeftCorner<3, 3>() =
        view_mat_.topLeftCorner<3, 3>().transpose();
    view_mat_inv_.block<3, 1>(0, 3) = camera_position;
  }

  void WorldToCam(const Vec3 &xyz, Vec3 *out) const {
    *out = view_mat_.topRows<3>() * xyz.homogeneous();
  }

  void WorldToCam(const Points3d &xyz, Points3d *out) const {
    *out = view_mat_.topRows<3>() * xyz.colwise().homogeneous();
  }

  void CamToWorld(const Vec3 &xyz, Vec3 *out) const {
    *out = view_mat_inv_.topRows<3>() * xyz.homogeneous();
  }

  void CamToWorld(const Points3d &xyz, Points3d *out) const {
    *out = view_mat_inv_.topRows<3>() * xyz.colwise().homogeneous();
  }

  void WorldToCamNormal(const Vec3 &xyz, Vec3 *out) const {
    *out = view_mat_.topLeftCorner<3, 3>() * xyz;
  }

  void WorldToCamNormal(const Points3d &xyz, Points3d *out) const {
    *out = view_mat_.topLeftCorner<3, 3>() * xyz;
  }

  void CamToWorldNormal(const Vec3 &xyz, Vec3 *out) const {
    *out = view_mat_inv_.topLeftCorner<3, 3>() * xyz;
  }

  void CamToWorldNormal(const Points3d &xyz, Points3d *out) const {
    *out = view_mat_inv_.topLeftCorner<3, 3>() * xyz;
  }

  // Camera coordinates to NDC.
  void CamToFrustum(const Vec3 &xyz, Vec3 *out) const {
    *out = (projection_mat() * xyz.homogeneous()).hnormalized();
  }

  void CamToFrustum(const Points3d &xyz, Points3d *out) const {
    *out = (projection_mat() * xyz.colwise().homogeneous()).colwise().hnormalized();
  }

  // NDC to Camera coordinates.
  void FrustumToCam(const Vec3 &xyz, Vec3 *out) const {
    *out = (projection_mat_inv() * xyz.homogeneous()).hnormalized();
  }

  void FrustumToCam(const Points3d &xyz, Points3d *out) const {
    *out = (projection_mat_inv() * xyz.colwise().homogeneous()).colwise().hnormalized();
  }

  // World coordinates to NDC.
  void WorldToFrustum(const Vec3 &xyz, Vec3 *out) const {
    *out = ((projection_mat() * view_mat_) * xyz.homogeneous()).hnormalized();
  }

  void WorldToFrustum(const Points3d &xyz, Points3d *out) const {
    *out = ((projection_mat() * view_mat_) * xyz.colwise().homogeneous()).colwise().hnormalized();
  }

  // NDC to World coordinates.
  void FrustumToWorld(const Vec3 &xyz, Vec3 *out) const {
    // For some reason (view_mat_inv_ * projection_mat_inv()) is unstable. They are mathematically the same.
    Mat44 pv_inv = (projection_mat() * view_mat_).inverse();  // TODO(daeyun): Avoid re-computing this.
    *out = (pv_inv * xyz.homogeneous()).hnormalized();
  }

  void FrustumToWorld(const Points3d &xyz, Points3d *out) const {
    Mat44 pv_inv = (projection_mat() * view_mat_).inverse();
    *out = (pv_inv * xyz.colwise().homogeneous()).colwise().hnormalized();
  }

  // `cam_depth_value` is optional.
  void CamToImage(const Vec3 &xyz, unsigned int height, unsigned int width, Vec2i *image_xy, double *cam_depth_value = nullptr) const {
    Vec3 ndc;
    CamToFrustum(xyz, &ndc);
    *image_xy = ((ndc.topRows<2>().array() + 1.0).colwise() * (Vec2{width, height} * 0.5).array()).cast<int>();
    image_xy->row(1) = height - image_xy->row(1).array() - 1;  // flip y.

    if (cam_depth_value) {
      // World-to-cam on z only.
      *cam_depth_value = -xyz[2];
    }
  }

  // `cam_depth_value` is optional.
  void CamToImage(const Points3d &xyz, unsigned int height, unsigned int width, Points2i *image_xy, Points1d *cam_depth_value = nullptr) const {
    Points3d ndc;
    CamToFrustum(xyz, &ndc);
    // https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glViewport.xhtml
    *image_xy = ((ndc.topRows<2>().array() + 1.0).colwise() * (Vec2{width, height} * 0.5).array()).cast<int>();
    image_xy->row(1) = height - image_xy->row(1).array() - 1;  // flip y.

    if (cam_depth_value) {
      // World-to-cam on z only.
      *cam_depth_value = -xyz.row(2);
    }
  }

  // `cam_depth_value` is optional.
  void WorldToImage(const Vec3 &xyz, unsigned int height, unsigned int width, Vec2i *image_xy, double *cam_depth_value = nullptr) const {
    Vec3 ndc;
    WorldToFrustum(xyz, &ndc);
    *image_xy = ((ndc.topRows<2>().array() + 1.0).colwise() * (Vec2{width, height} * 0.5).array()).cast<int>();
    image_xy->row(1) = height - image_xy->row(1).array() - 1;  // flip y.

    if (cam_depth_value) {
      // World-to-cam on z only.
      *cam_depth_value = -view_mat_.row(2) * xyz.homogeneous();
    }
  }

  // `cam_depth_value` is optional.
  void WorldToImage(const Points3d &xyz, unsigned int height, unsigned int width, Points2i *image_xy, Points1d *cam_depth_value = nullptr) const {
    Points3d ndc;
    WorldToFrustum(xyz, &ndc);
    *image_xy = ((ndc.topRows<2>().array() + 1.0).colwise() * (Vec2{width, height} * 0.5).array()).cast<int>();
    image_xy->row(1) = height - image_xy->row(1).array() - 1;  // flip y.

    if (cam_depth_value) {
      // World-to-cam on z only.
      *cam_depth_value = -view_mat_.row(2) * xyz.colwise().homogeneous();
    }
  }

  void ImageToWorld(const Points2i &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points3d *out) const {
    Points3d cam;
    ImageToCam(xy, cam_depth, height, width, &cam);
    CamToWorld(cam, out);
  }

  void ImageToCam(const Points2i &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points3d *out) const {
    Points2d xy_double = xy.cast<double>().array() + 0.5;
    ImageToCamTopLeftCorner(xy_double, cam_depth, height, width, out);
  }

  // Output is in the order of left, right, bottom, top, near, far.
  // Same ordering as http://docs.gl/gl3/glFrustum
  // Plane up direction is inward.
  void CamFrustumPlanes(array<Plane, 6> *out) const {
    Vec3 topleft = Vec3{frustum().left, frustum().top, -frustum().near}.normalized();
    Vec3 topright = Vec3{frustum().right, frustum().top, -frustum().near}.normalized();
    Vec3 bottomleft = Vec3{frustum().left, frustum().bottom, -frustum().near}.normalized();
    Vec3 bottomright = Vec3{frustum().right, frustum().bottom, -frustum().near}.normalized();

    Vec3 topleft2bottom = bottomleft - topleft;
    Plane left{
        .ref_point = topleft,
        .up_normal = topleft2bottom.cross(topleft).normalized(),
    };

    Vec3 topleft2right = topright - topleft;
    Plane top{
        .ref_point = topright,
        .up_normal = topright.cross(topleft2right).normalized(),
    };

    Vec3 bottomright2top = topright - bottomright;
    Plane right{
        .ref_point = bottomright,
        .up_normal = bottomright2top.cross(bottomright).normalized(),
    };

    Vec3 bottomright2left = bottomleft - bottomright;
    Plane bottom{
        .ref_point = bottomleft,
        .up_normal = bottomleft.cross(bottomright2left).normalized(),
    };

    Plane near{
        .ref_point = Vec3{0, 0, -frustum().near},
        .up_normal = Vec3{0, 0, -1},
    };

    Plane far{
        .ref_point = Vec3{0, 0, -frustum().far},
        .up_normal = Vec3{0, 0, 1},
    };

    *out = array < Plane, 6 > {
        left, right, bottom, top, near, far
    };
  }

  void WorldFrustumPlanes(array<Plane, 6> *out) const {
    CamFrustumPlanes(out);

    for (size_t i = 0; i < 6; ++i) {
      Vec3 world_ref_point;
      CamToWorld(out->at(i).ref_point, &world_ref_point);
      out->at(i).ref_point = world_ref_point;

      Vec3 world_normal;
      CamToWorldNormal(out->at(i).up_normal, &world_normal);
      out->at(i).up_normal = world_normal.normalized();
    }
  }

  const Mat44 &view_mat() const {
    return view_mat_;
  }

  const Mat44 &view_mat_inv() const {
    return view_mat_inv_;
  }

  const FrustumParams &frustum() const {
    return frustum_;
  }

  const Vec3 &position() const {
    return position_;
  }

  const Vec3 &lookat_position() const {
    return lookat_position_;
  }

  const Vec3 &viewing_direction() const {
    return viewing_direction_;
  }

  const Vec3 &up() const {
    return up_;
  }

  virtual void ImageToCamTopLeftCorner(const Points2d &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points3d *out) const = 0;
  virtual void PixelFootprintX(const Points2i &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points1d *out) const = 0;

  virtual const Mat44 &projection_mat() const = 0;
  virtual const Mat44 &projection_mat_inv() const = 0;
  virtual bool is_perspective() const = 0;

  // Returns false if this is orthographic. `x_fov` and `y_fov` are half of end-to-end angles, in radians.
  virtual bool fov(double *x_fov, double *y_fov) const = 0;

 private:
  Vec3 position_;
  Vec3 lookat_position_;
  Vec3 up_;
  Vec3 viewing_direction_;
  Mat44 view_mat_;
  Mat44 view_mat_inv_;
  FrustumParams frustum_;
};

class OrthographicCamera : public Camera {
 public:
  OrthographicCamera(const Vec3 &camera_position,
                     const Vec3 &lookat_position,
                     const Vec3 &up,
                     const FrustumParams &frustum_params)
      : Camera(camera_position, lookat_position, up, frustum_params) {
    auto rl = frustum().right - frustum().left;
    auto tb = frustum().top - frustum().bottom;
    auto fn = frustum().far - frustum().near;
    projection_mat_.setIdentity();
    projection_mat_(0, 0) = 2.0 / rl;
    projection_mat_(1, 1) = 2.0 / tb;
    projection_mat_(2, 2) = -2.0 / fn;
    projection_mat_(0, 3) = -(frustum().right + frustum().left) / rl;
    projection_mat_(1, 3) = -(frustum().top + frustum().bottom) / tb;
    projection_mat_(2, 3) = -(frustum().far + frustum().near) / fn;

    projection_mat_inv_ = projection_mat_.inverse();
  }

  // Skip NDC and directly find camera coordinates.
  void ImageToCamTopLeftCorner(const Points2d &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points3d *out) const override {
    Vec2 xy_scale{(frustum().right - frustum().left) / width, -(frustum().top - frustum().bottom) / height};
    Vec2 xy_offset{frustum().left, frustum().top};
    out->resize(3, xy.cols());

    out->topRows<2>() = (xy.array().colwise() * xy_scale.array()).array().colwise() + xy_offset.array();
    out->bottomRows<1>() = -cam_depth;
  }

  void PixelFootprintX(const Points2i &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points1d *out) const override {
    // Constant for a given orthographic frustum.
    double x_footprint = (frustum().right - frustum().left) / static_cast<double>(width);

    out->resizeLike(cam_depth);
    out->fill(x_footprint);
  }

  const Mat44 &projection_mat() const override {
    return projection_mat_;
  }

  const Mat44 &projection_mat_inv() const override {
    return projection_mat_inv_;
  }

  bool is_perspective() const override {
    return false;
  }

  bool fov(double *x_fov, double *y_fov) const override {
    return false;
  }

 private:
  Mat44 projection_mat_;
  Mat44 projection_mat_inv_;
};

class PerspectiveCamera : public Camera {
 public:
  PerspectiveCamera(const Vec3 &camera_position,
                    const Vec3 &lookat_position,
                    const Vec3 &up,
                    const FrustumParams &frustum_params)
      : Camera(camera_position, lookat_position, up, frustum_params) {
    // Expects symmetric frustum.
    Expects(std::abs(frustum().bottom + frustum().top) < 1e-7);  // Expects a symmetric frustum for perspective camera.
    Expects(std::abs(frustum().left + frustum().right) < 1e-7);

    auto rl = frustum().right - frustum().left;
    auto tb = frustum().top - frustum().bottom;
    auto fn = frustum().far - frustum().near;
    projection_mat_.setZero();
    projection_mat_(0, 0) = 2.0 * frustum().near / rl;
    projection_mat_(1, 1) = 2.0 * frustum().near / tb;
    projection_mat_(2, 2) = -(frustum().far + frustum().near) / fn;
    projection_mat_(0, 2) = (frustum().right + frustum().left) / rl;
    projection_mat_(1, 2) = (frustum().top + frustum().bottom) / tb;
    projection_mat_(2, 3) = -2.0 * frustum().far * frustum().near / fn;
    projection_mat_(3, 2) = -1.0;
    projection_mat_inv_ = projection_mat_.inverse();
  }

  // Skip NDC and directly find camera coordinates.
  void ImageToCamTopLeftCorner(const Points2d &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points3d *out) const override {
    Vec2 xy_scale{(frustum().right - frustum().left) / width, -(frustum().top - frustum().bottom) / height};
    Vec2 xy_offset{frustum().left, frustum().top};
    Points1d z = cam_depth.array() / frustum().near;
    out->resize(3, xy.cols());

    out->topRows<2>() = ((xy.array().colwise() * xy_scale.array()).array().colwise() + xy_offset.array()).array().rowwise() * z.array();
    out->bottomRows<1>() = -cam_depth;
  }

  void PixelFootprintX(const Points2i &xy, const Points1d &cam_depth, unsigned int height, unsigned int width, Points1d *out) const override {
    Points2d xy_double = xy.cast<double>();
    xy_double.row(1).array() += 0.5;
    Points3d out_left;
    ImageToCamTopLeftCorner(xy_double, cam_depth, height, width, &out_left);

    xy_double.row(0).array() += 1;
    Points3d out_right;
    ImageToCamTopLeftCorner(xy_double, cam_depth, height, width, &out_right);

    *out = (out_left - out_right).colwise().norm();
  }

  const Mat44 &projection_mat() const override {
    return projection_mat_;
  }

  const Mat44 &projection_mat_inv() const override {
    return projection_mat_inv_;
  }

  bool is_perspective() const override {
    return true;
  }

  bool fov(double *x_fov, double *y_fov) const override {
    // Frustum is expected to be symmetric.
    if (x_fov != nullptr) {
      *x_fov = std::abs(std::atan2(frustum().right, frustum().near));
    }
    if (y_fov != nullptr) {
      *y_fov = std::abs(std::atan2(frustum().top, frustum().near));
    }
    return true;
  }

 private:
  Mat44 projection_mat_;
  Mat44 projection_mat_inv_;
};

// `hw_ratio` is height/width.
FrustumParams ForceFixedAspectRatio(double hw_ratio, const FrustumParams &frustum);

}
