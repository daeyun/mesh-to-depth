//
// Created by daeyun on 4/27/18.
//

#pragma once

#include <functional>

#include "lib/common.h"

#include "nanort.h"

namespace mesh_to_depth {
class RayTracer {
 public:
  RayTracer(const std::vector<std::array<unsigned int, 3>> &faces,
            const std::vector<std::array<float, 3>> &vertices);

  void Traverse(const Vec3 &ray_origin, const Vec3 &ray_direction, const std::function<bool(float t, float u, float v, unsigned int prim_id)>& callback) const {
    nanort::Ray<float> ray;

    ray.org[0] = static_cast<float>(ray_origin[0]);
    ray.org[1] = static_cast<float>(ray_origin[1]);
    ray.org[2] = static_cast<float>(ray_origin[2]);

    ray.dir[0] = static_cast<float>(ray_direction[0]);
    ray.dir[1] = static_cast<float>(ray_direction[1]);
    ray.dir[2] = static_cast<float>(ray_direction[2]);

    std::map<unsigned int, float> ignored_prim_ids;
    nanort::BVHTraceOptions trace_options;
    trace_options.ignored_prim_ids = &ignored_prim_ids;  // Managed externally

    while (true) {
      nanort::TriangleIntersection<> isect{};
      bool hit = accel_.Traverse(ray, triangle_intersector_, &isect, trace_options);

      if (hit) {
        bool keep_going = callback(isect.t, isect.u, isect.v, isect.prim_id);

        if (!keep_going) {
          break;
        }

        ray.min_t = static_cast<float>(isect.t - kEps);

        if (!ignored_prim_ids.empty()) {
          for (auto it = ignored_prim_ids.cbegin(); it != ignored_prim_ids.cend();) {
            if (isect.t - it->second > kEps) {
              ignored_prim_ids.erase(it++);
            } else {
              ++it;
            }
          }
        }
        if (isect.t - trace_options.skip_prim_id_t < kEps) {
          ignored_prim_ids[trace_options.skip_prim_id] = trace_options.skip_prim_id_t;
        }
        trace_options.skip_prim_id = isect.prim_id;
        trace_options.skip_prim_id_t = isect.t;

      } else {
        break;
      }
    }
  }

  void PrintStats() const {
    nanort::BVHBuildStatistics stats = accel_.GetStatistics();

    printf("  BVH statistics:\n");
    printf("    # of leaf   nodes: %d\n", stats.num_leaf_nodes);
    printf("    # of branch nodes: %d\n", stats.num_branch_nodes);
    printf("  Max tree depth     : %d\n", stats.max_tree_depth);
    float bmin[3], bmax[3];
    accel_.BoundingBox(bmin, bmax);
    printf("  Bmin               : %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    printf("  Bmax               : %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  }

  Vec3 vertex(size_t i) const {
    return Vec3{vertices_[i][0], vertices_[i][1], vertices_[i][2]};
  }

  std::array<unsigned int, 3> face(size_t i) const {
    return {faces_[i][0], faces_[i][1], faces_[i][2]};
  }

  Vec3 normal(size_t face_i) const {
    auto f = this->face(face_i);
    Vec3 p1 = this->vertex(f[0]);
    Vec3 p2 = this->vertex(f[1]);
    Vec3 p3 = this->vertex(f[2]);
    return (p2 - p1).cross(p3 - p1).normalized();
  }

 private:
  std::vector<std::array<unsigned int, 3>> faces_;
  std::vector<std::array<float, 3>> vertices_;

  nanort::BVHAccel<float> accel_;
  nanort::TriangleIntersector<> triangle_intersector_;
  const float kEps = 0.0001;
};
}
