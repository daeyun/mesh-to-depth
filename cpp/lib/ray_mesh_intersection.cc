//
// Created by daeyun on 4/27/18.
//

#include "ray_mesh_intersection.h"

namespace mesh_to_depth {
RayTracer::RayTracer(const std::vector<std::array<unsigned int, 3>> &faces, const std::vector<std::array<float, 3>> &vertices)
    : faces_(faces), vertices_(vertices), triangle_intersector_(vertices_.data()->data(), faces_.data()->data(), sizeof(float) * 3) {
  nanort::TriangleMesh<float> triangle_mesh(vertices_.data()->data(), faces_.data()->data(), sizeof(float) * 3);
  nanort::TriangleSAHPred<float> triangle_pred(vertices_.data()->data(), faces_.data()->data(), sizeof(float) * 3);
  nanort::BVHBuildOptions<float> build_options;
  build_options.cache_bbox = true;  // ~1 second difference when loading suncg house models.

  bool build_ok = accel_.Build(static_cast<const unsigned int>(faces_.size()), triangle_mesh, triangle_pred, build_options);
  Ensures(build_ok);

}
}