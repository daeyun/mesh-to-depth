//
// Created by Daeyun Shin on 12/21/20.
//

#include <string>
#include "lib/camera.h"
#include "lib/ray_mesh_intersection.h"
#include "lib/depth_tracer.h"
//#include <omp.h>
#include "lib/common.h"

namespace mesh_to_depth {
/**
 * Generates multi-view, multi-layer depth representations of a mesh, potentially representing a scene containing multiple objects.
 *
 * @param vertices An array of shape (num_vertices, 3) and type float32.
 * @param num_vertices Number of vertices.
 * @param faces An array of shape (num_faces, 3) and type uint32.
 * @param num_faces Number of triangle mesh faces.
 * @param face_instances An array of shape (num_faces,) and type int32. Optional (can be null).
 *     jth value is the instance number of triangle face j. Negative value means impassible face. Negative values are only meaningful when doing multi-hit traversal for multi-layer depth map generation.
 *     A common use case is when each instance (aka group) represents an object.
 *     face_instances provides a many-to-one mapping face->instance.
 * @param camera_params An array of shape (num_cameras, 16) and type double.
 *     Each row represents a perspective or orthographic camera:
 *         cam_x, cam_y, cam_z,
 *         lookat_x, lookat_y, lookat_z,
 *         up_x, up_y, up_z,
 *         is_perspective,
 *         left, right, bottom, top, near, far
 *
 *     The first 9 parameters define the extrinsics of a camera. Same convention as [gluLookAt](https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluLookAt.xml).
 *     is_perspective: Determines whether the frustum is perspective or orthographic. A non-zero value means perspective.
 *     left, right, bottom, top, near, far: Define the intrinsics of a camera. Same convention as [glFrustum](https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glFrustum.xml) and [glOrtho](https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml).
 * @param num_cameras Number of cameras.
 * @param render_params An array of shape (num_cameras, 3) and type double.
 *     Row i contains the rendering parameters for the ith camera, casted to integer as needed:
 *         image_height, image_width, max_layers
 * @param out_ml_depth_values A pointer array of shape (num_cameras,) and type float32*.
 *     `out_ml_depth_values[i]` stores the values of a multi-layer depth map array of shape (h, w, \f$\max_{j} n_{ij}\f$) as seen from camera view i, in row-major order ignoring invalid values.
 *     Multi-layer depth map i is stored as a dynamically allocated 1d array of size \f$ \sum_{j=1}^{h_i*w_i} n_{ij} \f$  and type float32,
 *     where \f$n_ij\f$ is the number of values at pixel j of image i.
 *     Each depth map i is a multi-layer depth map when \f$ \exists n_{ij} in N_i: n_{ij} > 1 \f$.
 *     All dynamically allocated arrays must be freed by the caller.
 * @param out_pixel_num_layers A pointer array of shape (num_cameras,) and type uint32*.
 *     `out_pixel_num_layers[i]` is an array of shape (h_i, w_i) and type uint32, preallocated based on the image resolutions stored in `camera_params`.
 *     `out_pixel_num_layers[i]` is \f$N_i\f$.
 *     `out_pixel_num_layers[i][j]` is \f$n_{ij}\f$, i.e. the number of values at pixel j of image i, stored in row-major order.
 *     \f$n_{ij}\f$ can be zero if there is no valid value at that pixel.
 *     See \p out_ml_depth_values.
 * @param out_face_maps A pointer array of shape (num_cameras,) and type uint32*.
 *     Same spec as \p out_ml_depth_values except it stores the indices of the corresponding triangle face in \p face.
 *     All dynamically allocated arrays must be freed by the caller, same as \p out_ml_depth_values.
 *     Ignored if null.
 * @param out_instance_maps A pointer array of shape (num_cameras,) and type uint32*.
 *     Same as \p out_face_maps, except each value is transformed by the mapping face_instances[i_face]. See \p face_instances.
 *     All dynamically allocated arrays must be freed by the caller, same as \p out_ml_depth_values.
 *     Ignored if null.
 */
void MeshToDepth(
    const float *vertices,
    size_t num_vertices,
    const uint32_t *faces,
    size_t num_faces,
    const int32_t *face_instances,  // Optional input.
    const double *camera_params,
    size_t num_cameras,
    const double *render_params,
    float **out_ml_depth_values,  // Dynamically allocated.
    uint32_t **out_pixel_num_layers,
    uint32_t **out_face_maps,  // Optional output. Dynamically allocated.
    int32_t **out_instance_maps  // Optional output. Dynamically allocated.
) {
  // TODO(daeyun)
}

/**
 * Generates depth maps from a mesh, representing visible surfaces.
 *
 * @param vertices An array of shape (num_vertices, 3) and type float32.
 * @param num_vertices Number of vertices.
 * @param faces An array of shape (num_faces, 3) and type uint32.
 * @param num_faces Number of triangle mesh faces.
 * @param camera_params An array of shape (num_cameras, 16) and type double.
 *     Each row represents a perspective or orthographic camera:
 *         cam_x, cam_y, cam_z,
 *         lookat_x, lookat_y, lookat_z,
 *         up_x, up_y, up_z,
 *         is_perspective,
 *         left, right, bottom, top, near, far
 *
 *     The first 9 parameters define the extrinsics of a camera. Same convention as [gluLookAt](https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluLookAt.xml).
 *     is_perspective: Determines whether the frustum is perspective or orthographic. A non-zero value means perspective.
 *     left, right, bottom, top, near, far: Define the intrinsics of a camera. Same convention as [glFrustum](https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glFrustum.xml) and [glOrtho](https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml).
 * @param num_cameras Number of cameras.
 * @param render_params An array of shape (num_cameras, 3) and type double.
 *     Row i contains the rendering parameters for the ith camera, casted to integer as needed:
 *         image_height, image_width, is_depth
 *     If is_depth is false, ray displacement (t) values will be used instead.
 * @param out_depth_values A pointer array of shape (num_cameras,) and type float32*.
 */
void MeshToDepth_simple(
    const float *vertices,
    size_t num_vertices,
    const uint32_t *faces,
    size_t num_faces,
    const double *camera_params,
    size_t num_cameras,
    const double *render_params,
    float **out_depth_values  // Dynamically allocated.
) {
  std::vector<std::array<float, 3>> vertices_vector(num_vertices);
  for (int i = 0; i < num_vertices; ++i) {
    vertices_vector[i][0] = vertices[i * 3];
    vertices_vector[i][1] = vertices[i * 3 + 1];
    vertices_vector[i][2] = vertices[i * 3 + 2];
  }

  std::vector<std::array<unsigned int, 3>> faces_vector(num_faces);
  for (int i = 0; i < num_faces; ++i) {
    faces_vector[i][0] = faces[i * 3];
    faces_vector[i][1] = faces[i * 3 + 1];
    faces_vector[i][2] = faces[i * 3 + 2];
  }
  RayTracer ray_tracer(faces_vector, vertices_vector);

  std::vector<std::unique_ptr<Camera>> cameras;
  std::vector<std::array<size_t, 2>> height_width_vector;
  std::vector<bool> is_depth_vector;
  for (int i = 0; i < num_cameras; ++i) {
    bool is_perspective = camera_params[i * 16 + 9] != 0;

    // For now, only accept perspective cameras.
    Expects(is_perspective);

    unique_ptr<PerspectiveCamera> camera = std::make_unique<PerspectiveCamera>(
        Vec3{camera_params[i * 16 + 0], camera_params[i * 16 + 1], camera_params[i * 16 + 2]},
        Vec3{camera_params[i * 16 + 3], camera_params[i * 16 + 4], camera_params[i * 16 + 5]},
        Vec3{camera_params[i * 16 + 6], camera_params[i * 16 + 7], camera_params[i * 16 + 8]},
        FrustumParams{
            .left = camera_params[i * 16 + 10],
            .right = camera_params[i * 16 + 11],
            .bottom = camera_params[i * 16 + 12],
            .top = camera_params[i * 16 + 13],
            .near = camera_params[i * 16 + 14],
            .far = camera_params[i * 16 + 15],
        }
    );

    cameras.push_back(move(camera));
    height_width_vector.push_back({static_cast<size_t>(render_params[i * 3]), static_cast<size_t>(render_params[i * 3 + 1])});
    is_depth_vector.push_back(render_params[i * 3 + 2] != 0);
  }

  for (int i = 0; i < num_cameras; ++i) {
    const auto height = height_width_vector[i][0];
    const auto width = height_width_vector[i][1];
    const auto is_depth = is_depth_vector[i];

    // We want to be able to call free() in `free_arrays()`.
    out_depth_values[i] = static_cast<float *>(malloc(sizeof(float) * height * width));

    SimpleMultiLayerDepthRenderer depth_renderer(&ray_tracer, cameras[i].get(), width, height);

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        float out_value;
        uint32_t out_prim_id;
        bool is_hit = depth_renderer.DepthValue(x, y, is_depth, &out_value, &out_prim_id);
        if (is_hit) {
          out_depth_values[i][y * width + x] = out_value;
        } else {
          out_depth_values[i][y * width + x] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }

  }

}

/**
 * Used for development purposes, to make sure the python binding works, memory order is correct, and there's no memory leak.
 * Implements `meshgrid` in Numpy and Matlab.
 * @param num_images Number of output images.
 * @param resolutions An array of shape (num_images, 2), containing (h_i, w_i) resolutions.
 * @param out_y Array of y coordinates of grid i.
 * @param out_x Array of x coordinates of grid i.
 */
void MeshgridDebug(
    size_t num_images,
    const uint32_t *resolutions,
    uint32_t **out_y,
    uint32_t **out_x
) {
  for (int i = 0; i < num_images; ++i) {
    uint32_t h = resolutions[i * 2];
    uint32_t w = resolutions[i * 2 + 1];

    // We want to be able to call free() in `free_arrays()`.
    out_y[i] = static_cast<uint32_t *>(malloc(sizeof(uint32_t) * h * w));
    out_x[i] = static_cast<uint32_t *>(malloc(sizeof(uint32_t) * h * w));

    for (int j = 0; j < h; ++j) {
      for (int k = 0; k < w; ++k) {
        out_y[i][j * w + k] = j;
        out_x[i][j * w + k] = k;
      }
    }
  }
}
}

extern "C" {

void mesh2depth(
    const float *vertices,
    size_t num_vertices,
    const uint32_t *faces,
    size_t num_faces,
    const int32_t *face_instances,  // Optional input.
    const double *camera_params,
    size_t num_cameras,
    const double *render_params,
    float **out_ml_depth_values,  // Dynamically allocated.
    uint32_t **out_pixel_num_layers,
    uint32_t **out_face_maps,  // Optional output. Dynamically allocated.
    int32_t **out_instance_maps  // Optional output. Dynamically allocated.
) {
  // TODO(daeyun): Not implemented yet.
  mesh_to_depth::MeshToDepth(vertices, num_vertices, faces, num_faces, face_instances, camera_params, num_cameras, render_params, out_ml_depth_values, out_pixel_num_layers, out_face_maps, out_instance_maps);
}

void mesh2depth_simple(
    const float *vertices,
    size_t num_vertices,
    const uint32_t *faces,
    size_t num_faces,
    const double *camera_params,
    size_t num_cameras,
    const double *render_params,
    float **out_depth_values  // Dynamically allocated.
) {
  if (spdlog::get("console") == nullptr) {
    spdlog::stdout_color_mt("console");
  }
  mesh_to_depth::MeshToDepth_simple(vertices, num_vertices, faces, num_faces, camera_params, num_cameras, render_params, out_depth_values);
}

void _meshgrid_debug(
    size_t num_images,
    const uint32_t *resolutions,
    uint32_t **out_y,
    uint32_t **out_x
) {
  mesh_to_depth::MeshgridDebug(num_images, resolutions, out_y, out_x);
}

/**
 * Used to free the dynamically allocated output values of \p mesh2depth.
 * @param arr_ptr Array of malloc pointers. `arr_ptr` itself is not dynamically allocated.
 * @param size Size of the array provided by \p arr_ptr.
 */
void free_arrays(void **arr_ptr, size_t size) {
  if (arr_ptr) {
    for (int i = 0; i < size; ++i) {
      free(arr_ptr[i]);
    }
  }
}

}
